clear;
clc;

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
N0 = 1e-18;
tau = 2000; % Tx time per slot (2000ms)
tier2NumSlot = 4;
bandwidthhz = 180000; %kHz
Gamma = 1;
Phi = 10000;
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

% max SNR scheduling
imageSolution = zeros(tier2NumSlot, numNodes);
m_vecTemp = zeros(1,numNodes);
for i=1:numNodes
    if length(find(clusterHead==i)) == 1
        m_vecTemp(i) = -1;
    else
        if mod(find(clusterStructure==i),tier2NumSlot) == 0
            m_head = clusterStructure(tier2NumSlot,1);
        else
            m_head = clusterStructure(mod(find(clusterStructure==i),tier2NumSlot),1);
        end
        m_vecTemp(i) = length( find(matGij(:,m_head) > matGij(i,m_head)) );
    end
end
[interV interI] = sort(m_vecTemp,'descend')


