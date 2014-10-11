clear;
clc;

map = dlmread('../sourceData/20cam_r500_map.out');
idtEntropy = dlmread('../sourceData/20cam_r500_idt.out');
corrEntropy = dlmread('../sourceData/20cam_r500_corr.out');
algorithmFlag = 2;

% Parameter settings
numNodes = 20;
numClusters = 4;
radius = 900;
BSx = 0;
BSy = 0;
N0 = 1e-16; %1e-14
totalTime = 1000; % 10000ms
tier1Time = 900/numClusters; % 5000ms~9000ms
tier2NumSlot = 4;
tau = (totalTime - tier1Time)/tier2NumSlot; % Tier-2 time slot (2000ms)
bandwidthhz = 180000; %kHz
Gamma = 1;
%Phi = 10000;
maxPower = 1; % 1

% This is CSA (use v4)
if algorithmFlag == 1
    display('Clustering Algorithm: CSA');
    clusterHead = [6 2 18 7];
    clusterStructure = [6  8  17 3  0; ...
                        2  5  11 15 14; ...
                        18 1  10 19 16; ...
                        7  9  12 13 4]; % The order of head should match clusterHead
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
% This is direct access
elseif algorithmFlag == 4
    display('Direct access');
    clusterHead = [1:20];
    clusterStructure = [1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20];
    tier1Time = totalTime;
end

numHeads = length(clusterHead);

% Calculate channel gain
vecGi0 = zeros(1,numNodes);
matGij = zeros(numNodes,numNodes);
for i = 1:numNodes
    vecGi0(i) = CalChannelGain(map(i+1,1),map(i+1,2),BSx,BSy);
    for j = 1:numNodes
        matGij(i,j) = CalChannelGain(map(i+1,1),map(i+1,2),map(j+1,1),map(j+1,2));
    end
end

% Calculate joint entropy
jointEntropy = zeros(1,numHeads);
channelGainToBS = zeros(1,numHeads);
for i=1:numHeads
    m_nodes = clusterStructure(i,find(clusterStructure(i,:)>0));
    m_numNodes = length(m_nodes);
    m_corrMatrix = inf*ones(m_numNodes,m_numNodes);
    % Create partial correlation matrix
    for j=1:m_numNodes
        node_j = m_nodes(j);
        for k=1:m_numNodes
            node_k = m_nodes(k);
            if node_j == node_k
                m_corrMatrix(j,k) = idtEntropy(node_j);
            else
                m_corrMatrix(j,k) = corrEntropy(node_j,node_k);
            end
        end
    end
    % Start calculate joint entropy (greedy method)
    m_jointEntropy = 0;
    m_candidateList = [];
    for c=1:m_numNodes
        m_candidateList = [m_candidateList m_corrMatrix(c,c)];
    end
    [m_idtEntropy m_index] = sort(m_candidateList, 'ascend');
    m_jointEntropy = m_jointEntropy + m_idtEntropy(1);
    m_corrMatrix(m_index(1),:) = inf;
    for p=2:m_numNodes
        m_candidateList = m_corrMatrix(:,m_index(1))';
        [m_condEntropy m_index] = sort(m_candidateList, 'ascend');
        m_jointEntropy = m_jointEntropy + m_condEntropy(1);
        m_corrMatrix(m_index(1),:) = inf;
    end
    % Record each cluster's joint entropy
    jointEntropy(i) = m_jointEntropy;
    channelGainToBS(i) = vecGi0(clusterHead(i));
end

% Calculate r_j^*
rj = zeros(1,numHeads);
for i=1:numHeads
    m_head = clusterHead(i);
    m_otherHeads = clusterHead;
    m_otherHeads(i) = [];
    m_temp = (jointEntropy./channelGainToBS).^(1/2);
    rj(i) = ((channelGainToBS(i)*jointEntropy(i))^(1/2)/tier1Time)*(sum(m_temp));
end

% Calculate tier-1 power
pj = zeros(1,numHeads);
for i=1:numHeads
    pj(i) = (bandwidthhz*N0*(2^(rj(i)/bandwidthhz)-1))/channelGainToBS(i);    
end

totalTier1Power = sum(pj)