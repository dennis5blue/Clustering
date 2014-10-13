clear;
clc;
algorithmFlag = 1;
tau = 100; % Tx time per slot (100ms)
%tau = (1000 - 900)/4;
slotStart = 1;
slotEnd = 4;

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
    imageSolution = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1; ... 
                     0 0 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0; ...
                     0 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0 1 0 0 0; ...   
                     1 0 0 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0];
elseif algorithmFlag == 3
    display('Scheduling Algorithm: max idtEntropy');
    % This is maxIdtEntropy
    %                1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0
    imageSolution = [0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 1 1; ... 
                     0 0 0 0 1 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0; ...
                     1 0 1 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0; ...   
                     0 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 1 0 0 0];
elseif algorithmFlag == 4
    display('Scheduling Algorithm: greedy physical');
    % This is greedyPhysical
    %                1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0
    imageSolution = [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 1; ... % remove 4 (must add it back)
                     0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 1 0; ...
                     0 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 1 0 0 0; ...   
                     1 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0];
end

map = dlmread('../sourceData/20cam_r500_map.out');
idtEntropy = dlmread('../sourceData/20cam_r500_idt.out');
corrEntropy = dlmread('../sourceData/20cam_r500_corr.out');
CSRead = dlmread('../sourceData/CS_test_CSA_20_v4.out');
clusterHead = CSRead( 1,find(CSRead(1,:)>0) );
numClusters = length(clusterHead);
clusterStructure = CSRead(2:length(CSRead(:,1)),:);

% Parameter settings
numNodes = map(1,1);
numMembers = numNodes - length(clusterHead);
radius = map(1,2);
BSx = 0;
BSy = 0;
N0 = 1e-16; %1e-16
bandwidthhz = 180000; %kHz
Gamma = 1;
clusterMembers = [1:numNodes];
for i=1:length(clusterHead)
    m_head = clusterHead(i);
    nodeRm = find(clusterMembers==m_head);
    clusterMembers(nodeRm) = [];
end

% head for each nodes
vecNodeDerHead = zeros(1,numNodes);
for i=1:numNodes
    m_row = mod(find(clusterStructure == i),numClusters);
    if m_row == 0
        m_row = numClusters;
    end
    vecNodeDerHead(i) = clusterStructure(m_row,1);
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

powerEachSlot = [];

for focusSlot=slotStart:slotEnd
    activeNodes = find(imageSolution(focusSlot,:)==1);

    % Initial solution (does not consider interference)
    initialSol = zeros(1,length(activeNodes));
    for i=1:length(activeNodes)
        m_node = activeNodes(i);
        initialSol(i) = ( (bandwidthhz*N0)/matGij(m_node,vecNodeDerHead(m_node)) )*vecFi(m_node);
    end

    powerSol = initialSol;
    for iter=1:1000
        interference = powerSol;
        for i=1:length(activeNodes)
            m_node = activeNodes(i);
            m_interference = interference;
            m_channelGain = zeros(1,length(activeNodes));
            for j=1:length(activeNodes)
                m_interNode = activeNodes(j);
                m_channelGain(j) = matGij(m_interNode,vecNodeDerHead(m_node));
            end

            m_interference(i) = 0;
            m_channelGain(i) = 0;
            m_totalInterference = sum(m_interference.*m_channelGain);
            powerSol(i) = ( (bandwidthhz*N0+m_totalInterference)/matGij(m_node,vecNodeDerHead(m_node)) )*vecFi(m_node);
        end
        interference;
    end
    powerEachSlot = [powerEachSlot sum(powerSol)];
end

sum(powerEachSlot)