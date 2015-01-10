% Use largest eigenvalue to find SIR level and construct enumeration tree
% After that, go through enumeration tree and bound some branches
% Use the idea of SIC while calculating lower bound (treat former one as interferers)
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

vecSirLevel = zeros(1,numNodes);
for i = 1:numNodes
    vecSirLevel(i) = Gamma*2^( 8*idtEntropy(i)/(tau*bandwidthhz) - 1.0 );
end


%group1 = [2,3,8,10];
group1 = [10,2,3,8];
group2 = [4,5,6,9];
heads = [1,7];

% Slot 1
t = tree(group1(1));
branchIndex = 1;
vecBranchIndex = [];
for k = 1:length(group2)
    gammaMax = FindSinrBound(matGij,[group1(1) group2(k)],[heads(1) heads(2)]);
    if gammaMax > vecSirLevel(group1(1)) && gammaMax > vecSirLevel(group2(k)) % consider as feasilbe branch
        residualNodes = RemoveSched(group2,group2(k));
        eval(sprintf('[ t n_%i_%i_%i ] = t.addnode(1, [group2(%i) residualNodes]);',1,2,branchIndex,k));
        % node name namerule: n_(which slot)_(which cluster)_(which branch)
        % what info. we put in a node [schednode residualNodes ...] 
        branchIndex = branchIndex + 1;
    end
end
vecBranchIndex = [vecBranchIndex branchIndex-1];

for slot = 2:tier2NumSlot
    branchIndex = 1;
    for b = 1:vecBranchIndex(length(vecBranchIndex))
        eval(sprintf('[ t n_%i_%i_%i ] = t.addnode(n_%i_%i_%i, group1(%i));',slot,1,branchIndex,slot-1,2,b,slot));
        branchIndex = branchIndex + 1;
    end
    vecBranchIndex = [vecBranchIndex branchIndex-1];

    branchIndex = 1;
    for b = 1:vecBranchIndex(length(vecBranchIndex))
        eval(sprintf('nodesList = t.get(t.getparent(n_%i_%i_%i));',slot,1,b));
        if length(nodesList) >= 2 % since nodesList(1) is the scheduled machine
            for k = 2:length(nodesList)
                gammaMax = FindSinrBound(matGij,[group1(slot) nodesList(k)],[heads(1) heads(2)]);
                if gammaMax > vecSirLevel(group1(slot)) && gammaMax > vecSirLevel(nodesList(k))
                    residualNodes = RemoveSched(nodesList(2:length(nodesList)),nodesList(k));
                    eval(sprintf('[ t n_%i_%i_%i ] = t.addnode(n_%i_%i_%i, [nodesList(%i) residualNodes]);',slot,2,branchIndex,slot,1,b,k));
                    branchIndex = branchIndex + 1;
                end
            end
        end
    end
    vecBranchIndex = [vecBranchIndex branchIndex-1];
end

%disp(t.tostring)

% Start traverse tree and bound branches
objectiveUb = inf;
% Slot 1
solTree = tree(group1(1));
lbTree = tree(CalPower(group1(1),heads(1),[],[],matGij,vecSirLevel,bandwidthhz,N0));
traverseIndex = 1;
vecTraverseIndex = [];
childrenIndex = t.getchildren(1);
for i = 1:length(childrenIndex)
    Tx = t.get(childrenIndex(i));
    Tx = Tx(1);
    m_power = CalPower(Tx,heads(2),group1(1),lbTree.get(1),matGij,vecSirLevel,bandwidthhz,N0);
    eval(sprintf('[ solTree soln_%i_%i_%i ] = solTree.addnode(1, Tx);',1,2,traverseIndex));
    eval(sprintf('[ lbTree lbn_%i_%i_%i ] = lbTree.addnode(1, lbTree.get(1)+m_power);',1,2,traverseIndex));
    traverseIndex = traverseIndex + 1;
end
vecTraverseIndex = [vecTraverseIndex traverseIndex-1];

for slot = 2:tier2NumSlot
    traverseIndex = 1;
    m_power = CalPower(group1(slot),heads(1),[],[],matGij,vecSirLevel,bandwidthhz,N0);
    for b = 1:vecTraverseIndex(length(vecTraverseIndex))
        eval(sprintf('[ solTree soln_%i_%i_%i ] = solTree.addnode(soln_%i_%i_%i, group1(%i));',slot,1,traverseIndex,slot-1,2,b,slot));
        eval(sprintf('[ lbTree lbn_%i_%i_%i ] = lbTree.addnode(lbn_%i_%i_%i, lbTree.get(lbn_%i_%i_%i)+m_power);',slot,1,traverseIndex,slot-1,2,b,slot-1,2,b));
        traverseIndex = traverseIndex + 1;
    end
    vecTraverseIndex = [vecTraverseIndex traverseIndex-1];

    traverseIndex = 1;
    for b = 1:vecTraverseIndex(length(vecTraverseIndex))
        eval(sprintf('childrenIndex = t.getchildren(n_%i_%i_%i);',slot,1,b));
        for i = 1:length(childrenIndex)
            Tx = t.get(childrenIndex(i));
            Tx = Tx(1);
            eval(sprintf('m_power = CalPower(Tx,heads(2),solTree.get(soln_%i_%i_%i),lbTree.get(lbn_%i_%i_%i),matGij,vecSirLevel,bandwidthhz,N0);',slot,1,b,slot,1,b));
            eval(sprintf('[ solTree soln_%i_%i_%i ] = solTree.addnode(soln_%i_%i_%i, Tx);',slot,2,traverseIndex,slot,1,b));
            eval(sprintf('[ lbTree lbn_%i_%i_%i ] = lbTree.addnode(lbn_%i_%i_%i, lbTree.get(lbn_%i_%i_%i)+m_power);',slot,2,traverseIndex,slot,1,b,slot,1,b));
            traverseIndex = traverseIndex + 1;
        end
    end
    vecTraverseIndex = [vecTraverseIndex traverseIndex-1];
end

disp(solTree.tostring)
disp(lbTree.tostring)

% Start figures plotting
plotTree = lbTree;
nodeLeft = ones(1,length(plotTree.Node));
globalub = [1.6817]; % largest one
globallb = [cell2mat(plotTree.Node(1))];
globalIndex = [1];
nodeLeft(1) = 0;
treeIndex = 1;
while(sum(nodeLeft)>0)
    globalIndex = [globalIndex treeIndex];
    nodeLeft(treeIndex) = 0;
    if (plotTree.isleaf(treeIndex) == 1)
        globallb = [globallb globallb(length(globallb))];
        maxValue = cell2mat(plotTree.Node(treeIndex));
        [maxValue globalub(length(globalub))]
        if (maxValue < globalub(length(globalub)) && maxValue>=1.2942) % smallest one
            globalub = [globalub maxValue];
        else
            globalub = [globalub globalub(length(globalub))];
        end
        plotTree = plotTree.set(treeIndex,inf);
        treeIndex = plotTree.getparent(treeIndex);
        subTree = plotTree.subtree(treeIndex);
        mat_subTree = cell2mat(subTree.Node);
        while (length(plotTree.getchildren(treeIndex)) == 1 || all(mat_subTree==inf)==1 )
            treeIndex = plotTree.getparent(treeIndex);
            if treeIndex <= 0
                break;
            end
            subTree = plotTree.subtree(treeIndex);
            mat_subTree = cell2mat(subTree.Node);
        end        
    else
        globalub = [globalub globalub(length(globalub))];
        nextRunCandidate = plotTree.getchildren(treeIndex);
        nextRunCandidateValue = cell2mat(plotTree.Node(nextRunCandidate));
        minValue = min(nextRunCandidateValue);
        if (minValue > globallb(length(globallb)) && minValue < 1.2942) % smallest one
            globallb = [globallb minValue];
        else
            globallb = [globallb globallb(length(globallb))];
        end
        plotTree = plotTree.set(treeIndex,inf);
        treeIndex = nextRunCandidate(find(nextRunCandidateValue==minValue));
    end    
end

xaxis = [1:length(globallb)];
str_lb = sprintf('lower bound');
plot(xaxis, globallb, '-','LineWidth',2,'Color', 'g', 'DisplayName',str_lb); hold on;
str_ub = sprintf('upper bound');
plot(xaxis, globalub, '-','LineWidth',2,'Color', 'r', 'DisplayName',str_ub);
xlabel('Iteration');
ylabel('Tier-2 Power (Watt)');
legend('show','location','Best');
grid on;
hold off;