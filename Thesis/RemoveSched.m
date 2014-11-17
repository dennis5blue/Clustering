function [vec_others] = RemoveSched(CS,mat_sched)
    vec_others = [];
    for i=1:length(CS(:,1))
        vec_others = [vec_others CS(i,2:length(CS(i,:)))];
    end
    for i=1:length(mat_sched(:,1))
        for j=1:length(mat_sched(i,:))
            node = mat_sched(i,j);
            if length(find(vec_others==node))>0
                vec_others(find(vec_others==node)) = [];
            end
        end
    end
    vec_others(find(vec_others==0)) = [];
end