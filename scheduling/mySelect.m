function selection = mySelect(P,X)
    times = zeros(1,length(P));
    for i=1:length(P)
        times(i) = binornd(10,P(i));
    end
    [val index] = sort(times,'descend');
    selection = X(index(1));
end