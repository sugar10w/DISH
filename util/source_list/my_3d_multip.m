function [source_list_new] = my_3d_multip(source_list, aaa)
% easy multiplication

proj_cnt = numel(source_list);

source_list_new = cell(1, proj_cnt);

for i = 1:proj_cnt
    source_list_new{i} = source_list{i} * aaa;
end

end