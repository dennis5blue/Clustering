clear;
clc;
algorithmFlag = 3;
map = dlmread('../sourceData/20cam_r500_map.out');
idtEntropy = dlmread('../sourceData/20cam_r500_idt.out');
corrEntropy = dlmread('../sourceData/20cam_r500_corr.out');

% Parameter settings
numNodes = 20;
numHeads = 4;
numMembers = numNodes - numHeads;
radius = map(1,2);
BSx = 0;
BSy = 0;
N0 = 1e-14;
totalTime = 10000; % 10000ms
tier1Time = 9000; % 5000ms~9000ms
tier2NumSlot = 4;
tau = (totalTime - tier1Time)/tier2NumSlot; % Tier-2 reosurce (2000ms)
bandwidthhz = 180; %kHz
Gamma = 1;
%Phi = 10000;
maxPower = 1; % 1
epsilon = 0.01; % stop criteria for power update

% This is CSA (use v4)
if algorithmFlag == 1
    display('Clustering Algorithm: CSA');
    map = dlmread('../sourceData/20cam_r500_map.out');
    idtEntropy = dlmread('../sourceData/20cam_r500_idt.out');
    corrEntropy = dlmread('../sourceData/20cam_r500_corr.out');
    CSRead = dlmread('../sourceData/CS_test_CSA_20_v4.out');
    clusterHead = CSRead( 1,find(CSRead(1,:)>0) );
    clusterStructure = CSRead(2:length(CSRead(:,1)),:);

    % Parameter settings
    numNodes = map(1,1);
    numMembers = numNodes - length(clusterHead);
    tier2NumSlot = 4;
    bandwidthhz = 180; %kHz
    Gamma = 1;
    %Phi = 10000;
    maxPower = 1; % 1
    clusterMembers = [1:numNodes];
    for i=1:length(clusterHead)
        m_head = clusterHead(i);
        nodeRm = find(clusterMembers==m_head);
        clusterMembers(nodeRm) = [];
    end
    imageSolution = [0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 1 0 0 0; ...
                     0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 1 1; ...
                     0 0 1 1 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0; ...
                     1 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0];
    
% This is Kmeans (use v6)         
elseif algorithmFlag == 2
    display('Clustering Algorithm: Kmeans');
    clusterHead = [8 2 14 12];
    clusterStructure = [8  3  6  9  17 0  0  0  0; ...
                        2  1  5  10 11 15 16 18 19; ...
                        14 4  13 0  0  0  0  0  0; ...
                        7  12 20 0  0  0  0  0  0]; % The order of head should match clusterHead
            
% This is DMCP
elseif algorithmFlag == 3
    display('Clustering Algorithm: DMCP');
    clusterHead = [4 15 13 12];
    clusterStructure = [4  14 0  0  0  0  0  0  0;...
                        15 1  2  5  10 11 16 18 19;...
                        13 7  0  0  0  0  0  0  0;...
                        12 3  6  8  9  17 20 0  0]; % The order of head should match clusterHead
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

if algorithmFlag == 1 % Solve optimization problem to get tier-2 power; otherwise we need to use power update
    clusterMembers = [1:numNodes];
    for i=1:length(clusterHead)
        m_head = clusterHead(i);
        nodeRm = find(clusterMembers==m_head);
        clusterMembers(nodeRm) = [];
    end

    % Fill vecFi = 2^( H(V_i)/tau*W ) -1
    vecFi = zeros(1,numNodes);
    for i = 1:numNodes
        vecFi(i) = 2^( idtEntropy(i)/(tau*bandwidthhz) - 1.0 );
    end

    %
    %       ------------
    %       |  q1,1    |
    %       |  q1,2    |
    %       |    .     |
    %       |    .     |
    %       |    .     |
    %       |  q1,N    |
    % sol = |  q2,1    |
    %       |  q2,2    |
    %       |    .     |
    %       |    .     |
    %       |    .     |
    %       | q|S|,N   |
    %       ------------

    numConstraint = length(find(imageSolution==1));
    matConstraint = zeros(numConstraint,numNodes*tier2NumSlot);
    vecConstraint = zeros(numConstraint,1);
    conIndex = 1;
    for n=1:tier2NumSlot
        activeNode = find(imageSolution(n,:)==1);
        for i=1:length(activeNode)
            m_node = activeNode(i);
            temp = mod(find(clusterStructure==m_node),tier2NumSlot);
            if temp == 0
                m_head = clusterStructure(tier2NumSlot,1);
            else
                m_head = clusterStructure(temp,1);
            end
            matConstraint(conIndex,(m_node-1)*tier2NumSlot+n) = -matGij(m_node,m_head);
            m_vecInterNode = [];
            for k=1:length(activeNode)
                if activeNode(k) ~= m_node
                    m_interNode = activeNode(k);
                    matConstraint(conIndex,(m_interNode-1)*tier2NumSlot+n) = vecFi(m_node)*Gamma*matGij(m_interNode,m_head);
                end
            end
            vecConstraint(conIndex) = -vecFi(m_node)*Gamma*bandwidthhz*N0;
            conIndex = conIndex + 1;
        end
    end

    objective = ones(1,numNodes*tier2NumSlot);
    Aeq = zeros(1,numNodes*tier2NumSlot);
    beq = 0;
    lb = zeros(numNodes*tier2NumSlot,1);
    ub = inf*ones(numNodes*tier2NumSlot,1);

    [sol, fval, exitflag] = linprog(objective,matConstraint,vecConstraint,Aeq,beq,lb,ub);

    Q = zeros(numNodes,tier2NumSlot);
    Y = imageSolution';
    for i=1:numNodes
        for n=1:tier2NumSlot
            Q(i,n) = sol( (i-1)*tier2NumSlot + n ); % This is qin
        end
    end

    vec_power = [];
    for i=1:numNodes
        vec_power = [vec_power sum(Q(i,:).*Y(i,:))];
    end
    totalPower = sum(vec_power)
    
else % we use power update for these cases
    % Calculate vector C
    vecC = zeros(1,numNodes);
    for i=1:numNodes
        vecC(i) = Gamma*(2^(idtEntropy(i)/(tau*bandwidthhz) - 1.0));
    end
    
    iter = 0;
    updateStep = inf;
    vec_power = 0.00000000001*ones(1,numNodes);
    vec_deltaP = zeros(1,numNodes);
    
    % Initial state
    %for i=1:numNodes
    %    vec_power(i) = (bandwidthhz*N0+ max(vec_power.*m_InterfereGain) )/(matGij(m_node,m_head)/vecC(m_node))
    %end    
    
    while updateStep > epsilon
        for i=1:numNodes
            m_node = i;
            if mod(find(clusterStructure==i),numHeads) > 0
                m_head = clusterStructure( mod(find(clusterStructure==i),numHeads),1 );
            elseif mod(find(clusterStructure==i),numHeads) == 0
                m_head = clusterStructure( numHeads,1 );
            end
            m_interNodes = [];
            for h=1:numHeads
                if clusterStructure(h,1) ~= m_head
                    m_interNodes = [m_interNodes clusterStructure(h,2:length(find(clusterStructure(h,:)>0)))];
                end
            end

            % if a node is not in the set m_interNodes, then we set the channel
            % gain to 0, since it will not interfere traget node i
            m_InterfereGain = zeros(1,numNodes);
            for k=1:length(m_interNodes)
                node_k = m_interNodes(k);
                m_InterfereGain(node_k) = matGij(node_k,m_head);
            end
            oldvec_power = vec_power;
            oldvec_deltaP = vec_deltaP;
            vec_power(i) = (bandwidthhz*N0+ max(vec_power.*m_InterfereGain) )/(matGij(m_node,m_head)/vecC(m_node));
            vec_deltaP = vec_power - oldvec_power;
            vec_updateRatio = zeros(1,numNodes);
            for u=1:numNodes
                if vec_deltaP(u) > 0 & oldvec_deltaP > 0
                    vec_updateRatio(u) = vec_deltaP(u)/oldvec_deltaP(u);
                else
                    vec_updateRatio(u) = 0;
                end
            end
            vec_updateRatio;
            updateStep = sum(vec_updateRatio)/numNodes;
        end        
    end
    totalPower = sum(vec_power)
end