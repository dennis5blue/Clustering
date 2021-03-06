% deep first search
clc;
clear;
% Use matlan bintprog: we can get tier-2 power = 1.9903e+04
numNodes = 20;
tier2NumSlot = 4;
numClusters = 4;

matBranch = -1*ones(numNodes,tier2NumSlot); % use for fixing selection variables during the branching process 
convergeFlag = 0; % convergance happens when all elements in Y is binary
feasibleFlag = 0;
tempMin = inf; % record the minimum objective value obtained so far
iter = 0;
maxIter = 20000;
maxPath = 100;

%                --           --T
%                |   i    ...  |
% recordPath =   |   n    ...  |
%                |  0or1  ...  |
%                --           --

record_ObjValue = [inf; 0];

m_matBranch = matBranch;
m_matBranch(1) = 0;
k = 0;
while(iter<maxIter)
    iter = iter + 1;
    [tier2Power, convergeFlag, feasibleFlag] = LPsolver3(m_matBranch);
    [k convergeFlag]
    m_matBranch
    if k < numNodes*tier2NumSlot
        if convergeFlag == 0
            m_matBranch(k+1) = 0;
            k = k + 1;
            continue;
        elseif convergeFlag == 1
            while m_matBranch(k) == 1
                m_matBranch(k) = -1;
                k = k - 1;
            end
            m_matBranch(k) = 1;
            continue;
        end
    else
        while m_matBranch(k) == 1
            m_matBranch(k) = -1;
            k = k - 1;
        end
        m_matBranch(k) = 1;
        continue;
    end
end