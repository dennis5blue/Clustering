% Consider the interference SINR constraints
function [Q, Y] = LPsolver_relaxed(Y,map,idtEntropy,corrEntropy,clusterHead,clusterStructure,tier2NumSlot,powerMax,Phi,magicValue)
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
    %Phi = 1;
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
    constraintList = find(Y>0);
    equalityList0 = find(Y==0);
    equalityList1 = find(Y==1);
    numEquality = length(equalityList0)+length(equalityList1);
    numInterConstraint = length(constraintList);
    numScheConstraint = numNodes + numHeads*tier2NumSlot;
    vec_objectiveWeight = [ones(numNodes*tier2NumSlot,1) zeros(numNodes*tier2NumSlot,1)]; % First qin then yin
    mat_constraintWeight = zeros(numInterConstraint+numScheConstraint,numNodes*tier2NumSlot*2);
    vec_constraintValue = zeros(numInterConstraint+numScheConstraint,1);
    Aeq = zeros(numEquality,numNodes*tier2NumSlot*2);
    beq = zeros(numEquality,1);
    lb = zeros(numNodes*tier2NumSlot*2,1);
    ub = [powerMax*ones(numNodes*tier2NumSlot,1);ones(numNodes*tier2NumSlot,1)];
    m_clusterStructure = clusterStructure;
    
    % Fill in the interference constraint
    for c=1:numInterConstraint
        if mod(constraintList(c),numNodes)~=0
            m_Tx = mod(constraintList(c),numNodes);
        else
            m_Tx = numNodes;
        end
        m_slot = ceil(constraintList(c)/numNodes);
        if mod(find(m_clusterStructure==m_Tx),numHeads)~=0
            m_head = m_clusterStructure(mod(find(m_clusterStructure==m_Tx),numHeads),1);
        else
            m_head = m_clusterStructure(numHeads,1);
        end
        m_vecInter = m_clusterStructure(:,ceil(find(m_clusterStructure==m_Tx)/numHeads));
        m_vecInter(find(m_vecInter==m_Tx)) = [];
        
        for k=1:length(m_vecInter)
            m_interNode = m_vecInter(k);
            mat_constraintWeight( c,((m_interNode-1)*tier2NumSlot + m_slot) ) = vecFi(m_Tx)*matGij(m_interNode,m_head);
            if length(find(m_clusterStructure(:,1)==m_interNode))~=0
                display('Error: yin of cluster head must be assigned to 0');
            end
        end
        mat_constraintWeight( c,((m_Tx-1)*tier2NumSlot + m_slot) ) = -matGij(m_Tx,m_head);
        mat_constraintWeight(  c,(numNodes*tier2NumSlot + (m_Tx-1)*tier2NumSlot + m_slot) ) = Phi;
        vec_constraintValue(c) = -vecFi(m_Tx)*bandwidthhz*N0 + Phi;
    end
    
    % Fill in the scheduling constraint 1
    m_Tx = 1;
    for c=(numInterConstraint+1):(numInterConstraint+numNodes)
        for n=1:tier2NumSlot
            mat_constraintWeight( c,(numNodes*tier2NumSlot + (m_Tx-1)*tier2NumSlot + n) ) = magicValue;
        end
        vec_constraintValue(c) = magicValue;
        m_Tx = m_Tx + 1;
    end
    % Fill in the scheduling constraint 2
    c = numInterConstraint + numNodes + 1;
    for n=1:tier2NumSlot
        for j=1:numHeads
            m_TxList = m_clusterStructure(j,2:length(m_clusterStructure(j,:)));
            for i=1:length(m_TxList)
                m_Tx = m_TxList(i);
                mat_constraintWeight( c,(numNodes*tier2NumSlot + (m_Tx-1)*tier2NumSlot + n) ) = magicValue;
            end
            vec_constraintValue(c) = magicValue;
            c = c + 1;
        end
    end
    %mat_constraintWeight
    
    % Write down equality constarint (for fixing yin)
    for e=1:length(equalityList1)
        if mod(equalityList1(e),numNodes)~=0
            m_Tx = mod(equalityList1(e),numNodes);
        else
            m_Tx = numNodes;
        end
        m_slot = ceil(equalityList1(e)/numNodes);
        lb(numNodes*tier2NumSlot + (m_Tx-1)*tier2NumSlot + m_slot) = 1;
    end
    for e=1:length(equalityList0)
        if mod(equalityList0(e),numNodes)~=0
            m_Tx = mod(equalityList0(e),numNodes);
        else
            m_Tx = numNodes;
        end
        m_slot = ceil(equalityList0(e)/numNodes);
        ub(numNodes*tier2NumSlot + (m_Tx-1)*tier2NumSlot + m_slot) = 0;
    end
    
    
    % Transform sol to matrix Q
    sol = linprog(vec_objectiveWeight,mat_constraintWeight,vec_constraintValue,Aeq,beq,lb,ub);
    Q = zeros(numNodes,tier2NumSlot);
    Y = zeros(numNodes,tier2NumSlot);
    QderSol = sol(1:numNodes*tier2NumSlot);
    YderSol = sol((numNodes*tier2NumSlot+1):length(sol));
    for i=1:length(QderSol)
        if QderSol(i) ~= 0
            if mod(i,tier2NumSlot) ~= 0
                Q( floor(i/tier2NumSlot)+1,mod(i,tier2NumSlot) ) = QderSol(i);
                Y( floor(i/tier2NumSlot)+1,mod(i,tier2NumSlot) ) = YderSol(i);
            else
                Q( floor(i/tier2NumSlot),tier2NumSlot ) = QderSol(i);
                Y( floor(i/tier2NumSlot),tier2NumSlot ) = YderSol(i);
            end
        end
    end
    
end