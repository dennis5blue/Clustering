clc;
clear;
map = dlmread('../sourceData/10cam_r500_map.out');
idtEntropy = dlmread('../sourceData/10cam_r500_idt.out');
corrEntropy = dlmread('../sourceData/10cam_r500_corr.out');
CSRead = dlmread('../sourceData/CS_test_CSA_10.out');
clusterHead = CSRead( 1,find(CSRead(1,:)>0) );
clusterStructure = CSRead(2:length(CSRead(:,1)),:); % remember to pust cluster head at this 1st position

% These parameters must be changed inside LPsolver too!!!
numNodes = map(1,1);
numMembers = numNodes - length(clusterHead);
numHeads = length(clusterHead);
tier2NumSlot = 4;
powerMax = 1000; %mW

%Initial solution
Y = zeros(numNodes,tier2NumSlot);
for i=1:numHeads
    for j=1:tier2NumSlot
        Y(clusterStructure(i,j+1),j) = 1;
    end
end

[Q] = LPsolver(Y,map,idtEntropy,corrEntropy,clusterHead,clusterStructure,tier2NumSlot,powerMax)
sum(sum(Q))