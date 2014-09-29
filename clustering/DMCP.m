clear;
clc;

map = dlmread('../sourceData/20cam_r500_map.out');
idtEntropy = dlmread('../sourceData/20cam_r500_idt.out');
corrEntropy = dlmread('../sourceData/20cam_r500_corr.out');

% Parameter settings
numNodes = map(1,1);
numHeads = 4;
BSx = 0;
BSy = 0;
magicNum = 2;

% Calculate channel gain
vecGi0 = zeros(1,numNodes);
matGij = zeros(numNodes,numNodes);
for i = 1:numNodes
    vecGi0(i) = CalChannelGain(map(i+1,1),map(i+1,2),BSx,BSy);
    for j = 1:numNodes
        matGij(i,j) = CalChannelGain(map(i+1,1),map(i+1,2),map(j+1,1),map(j+1,2));
    end
end

% Calculate the nodes covered by each node
coverTopo = zeros(numNodes,numNodes);
for i=1:numNodes
    for j=1:numNodes
        if vecGi0(i) > matGij(i,j) % means i covers j
            coverTopo(i,j) = 1;
        end
    end
end
constCoverTopo = coverTopo;

% start DMCP
recordCoverTimes = zeros(1,numNodes);
metric = zeros(1,numNodes);
clusterHeads = zeros(1,numHeads);
m_numHeads = 0;

% Determine cluster heads
while m_numHeads < numHeads
    m_numHeads = m_numHeads + 1;
    for i=1:numNodes
        m_coverdNodesTemp = find(coverTopo(i,:)==1);
        m_coverdNodes = [];
        for k=1:length(m_coverdNodesTemp)
            node_k = m_coverdNodesTemp(k);
            if recordCoverTimes(node_k) < magicNum
                m_coverdNodes = [m_coverdNodes node_k];
            end
        end
        m_coverdNum = length(m_coverdNodes);
        if m_coverdNum > 0
            metric(i) = CalJointEmtropy(m_coverdNodes,idtEntropy,corrEntropy)/m_coverdNum;
        else
            metric(i) = 0;
        end
    end
    [m_val m_index] = sort(metric,'descend');
    clusterHeads(m_numHeads) = m_index(1);
    coverTopo(:,m_index(1)) = zeros(numNodes,1);
    coverTopo(m_index(1),:) = zeros(1,numNodes);
    
    % Update recordCoverTimes
    recordCoverTimes = recordCoverTimes + coverTopo(m_index(1),:);
end

% Determine cluster members
attachCondition = zeros(2,numNodes); % 1st row is node index, 2nd row is its cluster head
for i=1:numNodes
    node_i = i;
    channelToHead = [];
    for h=1:numHeads
        channelToHead = [channelToHead matGij(node_i,clusterHeads(h))];
    end
    [m_val m_index] = sort(channelToHead,'descend');
    attachCondition(1,node_i) = node_i;
    attachCondition(2,node_i) = clusterHeads(m_index(1));
end
attachCondition
