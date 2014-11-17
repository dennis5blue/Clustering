function [vec_heads] = FindHead(CS,vec_members)
    vec_heads = zeros(1,length(vec_members));
    for i=1:length(vec_members)
        node = vec_members(i);
        head_row = mod(find(CS==node),length(CS(:,1)));
        if head_row == 0
            vec_heads(i) = CS(length(CS(:,1)),1);
        else
            vec_heads(i) = CS(head_row,1);
        end
    end
end