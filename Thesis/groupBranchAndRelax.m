function objLowerBound = groupBranchAndRelax(mat_Tx,mat_Rx,vec_othersTx,vec_othersRx,mat_Gij,vec_Fi,bandwidth,N0,powerMax)
    objLowerBound = 0;
    numRuns = length(mat_Tx(:,1));
    vec_power = zeros(1,numRuns);
    for run=1:numRuns
        vec_Tx = mat_Tx(run,:);
        vec_Rx = mat_Rx(run,:);
        [feasibility, m_power] = powerSolver(vec_Tx,vec_Rx,mat_Gij,vec_Fi,bandwidth,N0,powerMax);
        if feasibility ~= 1
            objLowerBound = inf;
            break;
        else
            vec_power(run) = sum(m_power);
        end
    end
    %vec_power
    if feasibility == 1
        othersPower = 0;
        for i=1:length(vec_othersTx)
            m_Tx = vec_othersTx(i);
            m_Rx = vec_othersRx(i);
            othersPower = othersPower + ( vec_Fi(m_Tx)*bandwidth*N0 )/mat_Gij(m_Tx,m_Rx);
        end
        %othersPower
        objLowerBound = sum(vec_power) + othersPower;
    end
end