map = dlmread('../sourceData/24cam_r500_map.out');
idtEntropy = dlmread('../sourceData/24cam_r500_idt.out');
corrEntropy = dlmread('../sourceData/24cam_r500_corr.out');
clusterStructure = dlmread('../sourceData/CS_test_CSA_24.out');

% Parameter settings
numNodes = map(1,1);
radius = map(1,2);
BSx = 0;
BSy = 0;
tau = 2000; % Tx time per slot
tier2NumSlot = 5;
bandwidthKhz = 180; %kHz
Gamma = 1;
Phi = 10000;
maxPower = 1;

% Calculate channel gain
vecGi0 = zeros(1,numNodes);
matGij = zeros(numNodes,numNodes);
for i = 1:numNodes
    vecGi0(i) = CalChannelGain(map(i+1,1),map(i+1,2),BSx,BSy);
    for j = 1:numNodes
        matGij(i,j) = CalChannelGain(map(i+1,1),map(i+1,2),map(j+1,1),map(j+1,2));
    end
end

% Fill vecFi
vecFi = zeros(1,numNodes);
for i = 1:numNodes
    vecFi(i) = 2^( idtEntropy(i)/(tau*bandwidthKhz) - 1.0 );
end

% interference constraint
matIConstraint = zeros(numNodes,numNodes);
for i = 1:numNodes
    for j = 1:numNodes
        if (i==j)
            matIConstraint(i,j) = matGij(i,j);
        else
            matIConstraint(i,j) = -Fi(i)*Gamma*matGij(i,j);
        end
    end
end