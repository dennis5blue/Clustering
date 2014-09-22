clc;
clear;
% Use matlan bintprog: we can get tier-2 power = 1.9903e+04
numNodes = 24;
tier2NumSlot = 5;
numClusters = 4;

matBranch = -1*ones(24,5); % use for fixing selection variables during the branching process 
convergeFlag = 0; % convergance happens when all elements in Y is binary
tempMin = inf; % record the minimum objective value obtained so far
tempNotConvNum = numNodes*tier2NumSlot;
iter = 0;
maxIter = 2000;
iterMax = 5000;

%                --           --T
%                |   i    ...  |
% recordList =   |   n    ...  |
%                |  0or1  ...  |
%                --           --

record_ObjValue = [inf; 0];
for k=1:15 % k means which element to fix (e.g. k=1 means fix y11, k=2 means fix y21 ...)
    if iter>iterMax
        break;
    end
    k
    tempMin
    % Initial branch
    if k==1
        for i=1:2^k % i odd fix to 0, i even fix to 1 (note i grows in the order of 2)
            iter = iter + 1;
            m_matBranch = matBranch;
            if mod(i,2) == 1 % odd i fix to 0
                m_matBranch(k) = 0;
                [tier2Power, convergeFlag, feasibleFlag, numNotConvVariables] = LPsolver2(m_matBranch);
                if tier2Power < tempMin
                    tempMin = tier2Power;
                    record_ObjValue = [record_ObjValue [tempMin; iter]];
                end
                eval(['recordList_' num2str(k) '_' num2str(i) ' = [' num2str(k) '; 0];']);
            else % even i fix to 1
                m_matBranch(k) = 1;
                [tier2Power, convergeFlag, feasibleFlag, numNotConvVariables] = LPsolver2(m_matBranch);
                if tier2Power < tempMin
                    tempMin = tier2Power;
                    record_ObjValue = [record_ObjValue [tempMin; iter]];
                end
                eval(['recordList_' num2str(k) '_' num2str(i) ' = [' num2str(k) '; 1];']);
            end
        end
        
    % Branch based on the previous branch    
    else
        for i=1:2^k
            iter = iter + 1;
            % get the previous branched value
            m_matBranch = matBranch;
            prev_i = ceil(i/2);
            tempString =  ['recordList_' num2str(k-1) '_' num2str(prev_i)];
            if exist(tempString) == 1
                eval(['tempList = recordList_' num2str(k-1) '_' num2str(prev_i) ';']);
                for t=1:length(tempList(1,:))
                    m_matBranch(tempList(1,t)) = tempList(2,t);
                end
                        
                % start this run's branch
                if mod(i,2) == 1 % odd i fix to 0
                    m_matBranch(k) = 0;
                    [tier2Power, convergeFlag, feasibleFlag, numNotConvVariables] = LPsolver2(m_matBranch);
                    if feasibleFlag == 1
                        if tier2Power < tempMin
                            eval(['recordList_' num2str(k) '_' num2str(i) ' = [' num2str(k) '; 0];']);
                            eval(['recordList_' num2str(k) '_' num2str(i) ' = [tempList recordList_' num2str(k) '_' num2str(i) '];']);
                            tempMin = tier2Power;
                            record_ObjValue = [record_ObjValue [tempMin; iter]];
                        elseif tier2Power >= tempMin 
                            eval(['recordList_' num2str(k) '_' num2str(i) ' = [' num2str(k) '; 0];']);
                            eval(['recordList_' num2str(k) '_' num2str(i) ' = [tempList recordList_' num2str(k) '_' num2str(i) '];']);
                            tempNotConvNum = numNotConvVariables;
                        end
                    elseif feasibleFlag ~= 1
                        eval(['recordList_' num2str(k) '_' num2str(i) ' = [' num2str(k) '; 0];']);
                        eval(['recordList_' num2str(k) '_' num2str(i) ' = [tempList recordList_' num2str(k) '_' num2str(i) '];']);
                    end
                else % even i fix to 1
                    m_matBranch(k) = 1;
                    [tier2Power, convergeFlag, feasibleFlag, numNotConvVariables] = LPsolver2(m_matBranch);
                    if feasibleFlag == 1
                        if tier2Power < tempMin
                            eval(['recordList_' num2str(k) '_' num2str(i) ' = [' num2str(k) '; 1];']);
                            eval(['recordList_' num2str(k) '_' num2str(i) ' = [tempList recordList_' num2str(k) '_' num2str(i) '];']);
                            tempMin = tier2Power;
                            record_ObjValue = [record_ObjValue [tempMin; iter]];
                        elseif tier2Power >= tempMin
                            eval(['recordList_' num2str(k) '_' num2str(i) ' = [' num2str(k) '; 1];']);
                            eval(['recordList_' num2str(k) '_' num2str(i) ' = [tempList recordList_' num2str(k) '_' num2str(i) '];']);
                            tempNotConvNum = numNotConvVariables;
                        end
                    elseif feasibleFlag ~= 1
                        eval(['recordList_' num2str(k) '_' num2str(i) ' = [' num2str(k) '; 1];']);
                        eval(['recordList_' num2str(k) '_' num2str(i) ' = [tempList recordList_' num2str(k) '_' num2str(i) '];']);
                            
                        eval(['clear recordList_' num2str(k-1) '_' num2str(prev_i)]); % remember to free memory after two branches
                        continue;
                    end
                    eval(['clear recordList_' num2str(k-1) '_' num2str(prev_i)]); % remember to free memory after two branches
                end
            end
        end
    end
end

        