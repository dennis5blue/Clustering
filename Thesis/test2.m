% use largest eigenvalue to show SIR level
clc;
clear;
map = dlmread('../sourceData/10cam_r500_map.out');
idtEntropy = dlmread('../sourceData/10cam_r500_idt.out');
corrEntropy = dlmread('../sourceData/10cam_r500_corr.out');
CSRead = dlmread('../sourceData/CS_test_CSA_10.out');
clusterHead = CSRead( 1,find(CSRead(1,:)>0) );
CS = CSRead(2:length(CSRead(:,1)),:); % remember to pust cluster head at this 1st position

numNodes = map(1,1);
numMembers = numNodes - length(clusterHead);
numHeads = length(clusterHead);
tier2NumSlot = 4;
powerMax = 1000; %mW
Phi = 0;
N0 = 1e-16; %-16
Gamma = 1;
tau = 2; % Tx time per slot (s)
bandwidthhz = 180000; %Hz
    
matGij = zeros(numNodes,numNodes);
for i = 1:numNodes
    for j = 1:numNodes
        matGij(i,j) = CalChannelGain(map(i+1,1),map(i+1,2),map(j+1,1),map(j+1,2));
    end
end

vecFi = zeros(1,numNodes);
for i = 1:numNodes
    vecFi(i) = Gamma*2^( 8*idtEntropy(i)/(tau*bandwidthhz) - 1.0 );
end

%[AA,BB] = powerSolver([3,6],[1,7],matGij,vecFi,bandwidthhz,N0,powerMax)
%obj = groupBranchAndRelax([3,6;2,4],[1,7;1,7],[5,8,9,10],[7,1,7,1],matGij,vecFi,bandwidthhz,N0,powerMax)
%RemoveSched(CS,[3,6;2,4])
%FindHead(CS,RemoveSched(CS,[3,6;2,4]))

m_Tx = [2,5;10,4;3,9;8,6]
m_Rx = [1,7;1,7;1,7;1,7]
m_others = RemoveSched(CS,m_Tx)
%obj = groupBranchAndRelax(m_Tx,m_Rx,m_others,FindHead(CS,m_others),matGij,vecFi,bandwidthhz,N0,powerMax)