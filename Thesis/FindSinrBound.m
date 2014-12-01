function [lambdaMax, vec_power] = FindSinrBound(mat_Gij,vec_Tx,vec_Rx)
    lambdaMax = 0;
    vec_power = [];
    mat_normalizedGain = zeros(length(vec_Tx),length(vec_Tx));
    for i=1:length(vec_Tx)
        for j=1:length(vec_Tx)
            mat_normalizedGain(i,j) = mat_Gij(vec_Tx(j),vec_Rx(i))/mat_Gij(vec_Tx(i),vec_Rx(i));
        end
    end
    [V,D] = eig(mat_normalizedGain); % A*V=V*D,
    for i=1:length(D(1,:))
        if D(i,i)>lambdaMax
            lambdaMax = D(i,i);
            vec_power = V(:,i);
        end
    end
end