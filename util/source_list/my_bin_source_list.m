function [source_list_new] = my_bin_source_list(source_list, th)
% �򵥵ض� source_list ����ֵ��
%  Binarize source_list in a simple way

proj_cnt = numel(source_list);
source_list_new = cell(1, proj_cnt);

for proj_i = 1:numel(source_list)
    source_list_new{proj_i} = source_list{proj_i} > th;
end


end

