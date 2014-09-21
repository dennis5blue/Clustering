% Not used, please run BBsolver2.m
numNodes = 24;
tier2NumSlot = 5;
numClusters = 4;
obj1 = inf;
obj0 = inf;

matBranch = -1*ones(24,5); % use for fixing selection variables during the branching process 
iter = 0; % use for recording iterations to check for convergence
convergeFlag = 0; % convergance happens when all elements in Y is binary
tempMin = inf; % record the minimum objective value obtained so far
maxIter = 10000000;

obj1=BBsolver(iter,1,1,matBranch,1,tempMin,maxIter);
obj0=BBsolver(iter,1,1,matBranch,0,tempMin,maxIter);
min(obj1,obj0)
