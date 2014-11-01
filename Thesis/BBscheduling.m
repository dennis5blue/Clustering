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
Phi = 0;
magicValue = 1e-07;

%Initial solution
Y = zeros(numNodes,tier2NumSlot);
for i=1:numHeads
    for j=1:tier2NumSlot
        Y(clusterStructure(i,j+1),j) = 1;
    end
end

Y = 0.5*ones(numNodes,tier2NumSlot);
Y(2,:) = [1 0 0 0];
Y(1,:) = [0 0 0 0];
Y(7,:) = [0 0 0 0];

%[Q] = LPsolver(Y,map,idtEntropy,corrEntropy,clusterHead,clusterStructure,tier2NumSlot,powerMax)
%minValue = sum(sum(Q))
[Q, Y] = LPsolver_relaxed(Y,map,idtEntropy,corrEntropy,clusterHead,clusterStructure,tier2NumSlot,powerMax,Phi,magicValue);
minValue = sum(sum(Q))
