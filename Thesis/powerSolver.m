% solve transmission power given a set of transmitters
function [feasibility, vec_power] = powerSolver(vec_Tx,vec_Rx,mat_Gij,vec_Fi,bandwidth,N0,powerMax)
    numVariables = length(vec_Tx);
    numConstraints = length(vec_Tx);
    
    mat_constraintWeight = zeros(numVariables,numConstraints);
    vec_constraintValue = zeros(numConstraints,1);
    vec_objectiveWeight = ones(numVariables,1);
    lb = zeros(numVariables,1);
    for i=1:numConstraints
        m_Tx = vec_Tx(i);
        m_Rx = vec_Rx(i);
        m_interfer = vec_Tx;
        m_interfer(find(m_interfer==m_Tx)) = [];
        %m_Tx
        %m_Rx
        %m_interfer
        mat_constraintWeight(i,find(vec_Tx==m_Tx)) = -mat_Gij(m_Tx,m_Rx);
        for k=1:length(m_interfer)
            mat_constraintWeight(i,find(vec_Tx==m_interfer(k))) = vec_Fi(m_Tx)*mat_Gij(m_interfer(k),m_Rx);
        end
        vec_constraintValue(i) = -vec_Fi(m_Tx)*bandwidth*N0;
        %mat_constraintWeight
        %vec_constraintValue
    end
    [x,fval,exitflag] = linprog(vec_objectiveWeight,mat_constraintWeight,vec_constraintValue,[],[],lb);
    feasibility = exitflag;
    if length(find(x>powerMax)) > 0
        feasibility = 0;
    end
    vec_power = x';    
end