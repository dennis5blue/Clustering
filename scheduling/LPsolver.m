function sol = LPsolver(iToFix,nToFix)

    map = dlmread('../sourceData/24cam_r500_map.out');
    idtEntropy = dlmread('../sourceData/24cam_r500_idt.out');
    corrEntropy = dlmread('../sourceData/24cam_r500_corr.out');
    CSRead = dlmread('../sourceData/CS_test_CSA_24.out');
    clusterHead = CSRead( 1,find(CSRead(1,:)>0) );
    clusterStructure = CSRead(2:length(CSRead(:,1)),:);

    % Parameter settings
    numNodes = map(1,1);
    numMembers = numNodes - length(clusterHead);
    radius = map(1,2);
    BSx = 0;
    BSy = 0;
    N0 = 1e-16;
    tau = 2000; % Tx time per slot
    tier2NumSlot = 5;
    bandwidthhz = 180000; %kHz
    Gamma = 1;
    Phi = 0;
    maxPower = 1;
    clusterMembers = [1:numNodes];
    for i=1:length(clusterHead)
        m_head = clusterHead(i);
        nodeRm = find(clusterMembers==m_head);
        clusterMembers(nodeRm) = [];
    end

    % Calculate channel gain
    vecGi0 = zeros(1,numNodes);
    matGij = zeros(numNodes,numNodes);
    for i = 1:numNodes
        vecGi0(i) = CalChannelGain(map(i+1,1),map(i+1,2),BSx,BSy);
        for j = 1:numNodes
            matGij(i,j) = CalChannelGain(map(i+1,1),map(i+1,2),map(j+1,1),map(j+1,2));
        end
    end

    % Fill vecFi = 2^( H(V_i)/tau*W ) -1
    vecFi = zeros(1,numNodes);
    for i = 1:numNodes
        vecFi(i) = 2^( idtEntropy(i)/(tau*bandwidthhz) - 1.0 );
    end

    % Initial the constraints matrix and bounds of selection variables
    numIC = (numNodes - length(clusterHead))*tier2NumSlot; % number of interference constraints
    numSC = length(clusterHead)*tier2NumSlot + (numNodes - length(clusterHead)); % number of selection constraints
    numSV = 2*numNodes*tier2NumSlot; % number of selection variables
    matConstraints = zeros((numIC+numSC),numSV);
    vecConstraints = zeros((numIC+numSC),1);
    Aeq = zeros( (numMembers+numMembers*tier2NumSlot),numSV );
    beq = zeros( (numMembers+numMembers*tier2NumSlot),1 );
    lb = zeros(numSV,1);
    ub = inf*ones(numSV,1);
    objective = zeros(1,numSV);

    %
    % y1,1
    % y1,2
    % y1,3
    % y2,1
    % y2,2
    % y2,3
    % q1,1
    % q1,2
    % q1,3
    % q2,1
    % q2,2
    % q2,3
    %

    % interference constraint: \Psi Y - A^T Q \leq \Phi-B
    conIndex = 1; % constarint index
    for i=1:length(clusterMembers)
        m_node = clusterMembers(i);
        m_vecInterNode = [];
        for j=1:length(clusterHead)
            if (find(clusterStructure(j,:)==m_node))
                m_head = clusterHead(j);
            else
                m_vecInterNode = [m_vecInterNode clusterStructure( j,2:length(find(clusterStructure(j,:)>0)))];
            end
        end
    
        for n=1:tier2NumSlot
            matConstraints(conIndex, (numSV/2 + (m_node-1)*tier2NumSlot +n) ) = -matGij(m_node,m_head); %q_in
            for k=1:length(m_vecInterNode)
                m_interNode = m_vecInterNode(k);
                matConstraints(conIndex, (numSV/2 + (m_interNode-1)*tier2NumSlot +n) ) = ...
                    vecFi(m_node)*Gamma*matGij(m_interNode,m_head); % q_kn
            end
            matConstraints(conIndex, ( (m_node-1)*tier2NumSlot + n )) = Phi; % y_in
            vecConstraints(conIndex) = Phi - vecFi(m_node)*Gamma*bandwidthhz*N0;
            conIndex = conIndex + 1;
        end
    end

    % selection constraint 1 (for cluster)
    for j=1:length(clusterHead)
        m_head = clusterHead(j);
        m_vecNode = clusterStructure(j,2:length(clusterStructure(j,:)));
        for n=1:tier2NumSlot
            for i=1:length(m_vecNode)
                m_node = m_vecNode(i);
                matConstraints(conIndex,( (m_node-1)*tier2NumSlot + n ))=1;
            end
            vecConstraints(conIndex) = 1;
            conIndex = conIndex + 1;
        end
    end
    % selection constraint 2 (for time slot)
    for i=1:length(clusterMembers)
        m_node = clusterMembers(i);
        for n=1:tier2NumSlot
            matConstraints(conIndex, ( (m_node-1)*tier2NumSlot + n ))=1;
        end
        vecConstraints(conIndex) = 1;
        conIndex = conIndex + 1;
    end

    % set the equality constraint (make sure all the nodes are supported)
    % mouda does not have this equaltiy constraint
    eqIndex = numMembers*tier2NumSlot+1;
    Aeq( ((iToFix-1)*tier2NumSlot + nToFix),((iToFix-1)*tier2NumSlot + nToFix) ) = 1;
    beq((iToFix-1)*tier2NumSlot + nToFix) = 1;
    for i=1:length(clusterMembers)
        m_node = clusterMembers(i);
        for n=1:tier2NumSlot
            Aeq(eqIndex, ( (m_node-1)*tier2NumSlot + n )) = 1;
        end
        beq(eqIndex) = 1;
        eqIndex = eqIndex + 1;
    end

    % set the bounds of selection variables
    for i=1:(numSV/2)
        ub(i) = 1;
    end

    % set the objective function
    for i=1:length(clusterMembers)
        m_node = clusterMembers(i);
        for n=1:tier2NumSlot
            objective( numSV/2 + (m_node-1)*tier2NumSlot +n ) = 1;
        end
    end

    sol = linprog(objective,matConstraints,vecConstraints,Aeq,beq,lb,ub);
end