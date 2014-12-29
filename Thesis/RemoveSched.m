function [vec_residual] = RemoveSched(nodeList,nodeRm)
    vec_residual = nodeList;
    for i=1:length(vec_residual)
        if vec_residual(i) == nodeRm
            vec_residual(i) = -1;
        end
    end
    vec_residual(find(vec_residual == -1)) = [];
end