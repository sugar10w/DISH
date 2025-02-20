function [source_list_new] = my_binarize_source_list_dynamic_th(source_list, proj_ratio, max_intensity)
%MY_BINARIZE_SOURCE_LIST_DYNAMIC_TH 使用动态阈值进行二值化


if nargin<3
    max_intensity = my_get_source_list_max(source_list);
end
if nargin<2
    proj_ratio = 3;
end


proj_cnt = numel(source_list);

% 生成阈值列表
th_list = zeros(1, proj_cnt);
for i = 1:proj_cnt
    ii = mod(i-1, proj_ratio);
    th_list(i) = (ii + 0.5) / proj_ratio;
end
th_list = th_list * max_intensity;

% 生成二值图像
source_list_new = cell(1, proj_cnt);
for i = 1:proj_cnt
    source_list_new{i} = source_list{i} > th_list(i);
end

end

