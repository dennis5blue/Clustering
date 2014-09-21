% Not used, please run BBsolver2.m
numNodes = 24;
tier2NumSlot = 5;
numClusters = 4;

matBranch = -1*ones(24,5); % use for fixing selection variables during the branching process 
iter = 0; % use for recording iterations to check for convergence
convergeFlag = 0; % convergance happens when all elements in Y is binary
tempMin = inf; % record the minimum objective value obtained so far
maxIter = 2000;

%                --           --T
%                |   i    ...  |
% recordList =   |   n    ...  |
%                |  0or1  ...  |
%                --           --

% Fix to 1 first
fixVal = 0;
candidate = [];
for i=1:numNodes
    for n=1:tier2NumSlot
        % check if the maximum iteration time reaches
        if iter>maxIter
            disp('Maximum iteration reaches')
            break;
        end
        iter = iter + 1;
        
        % solve the relaxed LP
        m_matBranch = matBranch;
        m_matBranch(i,n) = fixVal;
        [tier2Power, convergeFlag, feasibleFlag] = LPsolver(m_matBranch);
        
        if ( feasibleFlag == 1 && convergeFlag ~= 1)
            tempMin = tier2Power;
            eval(['recordList_' num2str(iter) ' = [' num2str(i) ';' num2str(n) ';' num2str(fixVal) '];']);
            candidate = [candidate iter];
            disp('Branch')
            continue;
        end
        
        % infeasible
        if feasibleFlag ~= 1
            matBranch(i,n) = -1;
            disp('Bound: no more feasible solution');
            break;
        end
        
        %if tier2Power > tempMin
        %    matBranch(i,n) = -1;
        %    disp('Bound: no more better objective value');
        %    break;
        %end
        
        % convergance
        if convergeFlag == 1
            matBranch(i,n) = -1;
            tempMin = tier2Power;
            disp('Bound: convergance reached');
            break;
        end
        
    end
end
%