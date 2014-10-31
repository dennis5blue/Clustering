% Consider the interference SINR constraints
function [Q] = LPsolver(Y,map,idtEntropy,corrEntropy,clusterHead,clusterStructure,tier2NumSlot,powerMax)
    % Parameter settings
    numNodes = map(1,1);
    numMembers = numNodes - length(clusterHead);
    numHeads = length(clusterHead);
    radius = map(1,2);
    BSx = 0;
    BSy = 0;
    N0 = 1e-15;
    tau = 2000; % Tx time per slot (ms)
    bandwidthhz = 180000; %Hz
    Gamma = 1;
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
        vecFi(i) = Gamma*2^( idtEntropy(i)/(tau*bandwidthhz) - 1.0 );
    end
    
    % Write Interference constraints based on selection matrix Y
    vec_objectiveWeight = ones(numNodes*tier2NumSlot,1);
    mat_constraintWeight = zeros(numHeads*tier2NumSlot,numNodes*tier2NumSlot);
    vec_constraintValue = zeros(numHeads*tier2NumSlot,1);
    Aeq = zeros(numHeads*tier2NumSlot,numNodes*tier2NumSlot);
    beq = zeros(numHeads*tier2NumSlot,1);
    lb = zeros(numNodes*tier2NumSlot,1);
    ub = powerMax*ones(numNodes*tier2NumSlot,1);
    m_clusterStructure = zeros(numHeads,tier2NumSlot+1); % determined by Y
    m_clusterStructure(:,1) = clusterStructure(:,1);
    for i=1:tier2NumSlot
        schedNodes = find(Y(:,i)==1);
        for j=1:length(schedNodes)
            m_node = schedNodes(j);
            m_whichCluster = mod(find(clusterStructure == m_node),numHeads);
            if m_whichCluster==0
                m_whichCluster = numHeads;
            end
            m_clusterStructure(m_whichCluster,i+1) = m_node;
        end
    end
    
    % Start LP formulation
    for n=1:tier2NumSlot
        for j=1:numHeads
            constraintIndex = numHeads*(n-1)+j;
            m_head = m_clusterStructure(j,1);
            m_Tx = m_clusterStructure(j,n+1);
            m_vecInter = m_clusterStructure(:,n+1);
            m_vecInter(j)=[];
            for k=1:length(m_vecInter)
                m_interNode = m_vecInter(k);
                mat_constraintWeight( constraintIndex,((m_interNode-1)*tier2NumSlot + n) ) = vecFi(m_Tx)*matGij(m_interNode,m_head);
            end
            mat_constraintWeight( constraintIndex,((m_Tx-1)*tier2NumSlot + n) ) = -matGij(m_Tx,m_head);
            vec_constraintValue(constraintIndex) = -vecFi(m_Tx)*bandwidthhz*N0;
        end
    end
    
    % Transform sol to matrix Q
    sol = linprog(vec_objectiveWeight,mat_constraintWeight,vec_constraintValue,Aeq,beq,lb,ub);
    Q = zeros(numNodes,tier2NumSlot);
    for i=1:length(sol)
        if sol(i) ~= 0
            if mod(i,tier2NumSlot) ~= 0
                Q( floor(i/tier2NumSlot)+1,mod(i,tier2NumSlot) ) = sol(i); 
            else
                Q( floor(i/tier2NumSlot),tier2NumSlot ) = sol(i);
            end
        end
    end
    
end