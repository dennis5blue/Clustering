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

group1 = [2,3,8,10];
group2 = [4,5,6,9];
heads = [1,7];
i = 4;
for j=1:length(group2)
    transmitters = [group1(i),group2(j)]
    [lambdaMax,vec_power] = FindSinrBound(matGij,[group1(i),group2(j)],[1,7]);
    [1/(lambdaMax-1) sum(vec_power)]
end