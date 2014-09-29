clear;
clc;
algorithmFlag = 4;

if algorithmFlag == 1
    display('Scheduling Algorithm: branch-and-bound');
    % This is branch-and-bound
                   % 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0
    imageSolution = [0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 1 0 0 0; ...
                     0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 1 1; ...
                     0 0 1 1 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0; ...
                     1 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0];
elseif algorithmFlag == 2
    display('Scheduling Algorithm: max SNR');
    % This is maxSNR 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0
    imageSolution = [1 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0; ... 
                     0 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 1 0 0 0; ...
                     0 0 0 0 1 0 0 0 1 0 0 1 0 1 0 0 0 0 0 0; ...   
                     0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1];
elseif algorithmFlag == 3
    display('Scheduling Algorithm: max sum rate');
    % This is maxSumate
    imageSolution = [1 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0; ... 
                     1 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0; ...
                     1 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0; ...   
                     1 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0];
elseif algorithmFlag == 4
    display('Scheduling Algorithm: greedy physical');
    % This is greedyPhysical
    %                1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0
    imageSolution = [0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 1; ... 
                     0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 1 0; ...
                     0 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 1 0 0 0; ...   
                     1 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0];
end

map = dlmread('../sourceData/20cam_r500_map.out');
idtEntropy = dlmread('../sourceData/20cam_r500_idt.out');
corrEntropy = dlmread('../sourceData/20cam_r500_corr.out');
CSRead = dlmread('../sourceData/CS_test_CSA_20_v4.out');
clusterHead = CSRead( 1,find(CSRead(1,:)>0) );
clusterStructure = CSRead(2:length(CSRead(:,1)),:);

% Parameter settings
numNodes = map(1,1);
numMembers = numNodes - length(clusterHead);
radius = map(1,2);
BSx = 0;
BSy = 0;
N0 = 1e-14;
tau = 5000; % Tx time per slot (2000ms)
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