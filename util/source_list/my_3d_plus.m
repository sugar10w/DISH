function [source_list_new] = my_3d_plus(source_list_1, source_list_2)
% my_3d_plus 将不同的source_list 加起来
%  add up different source_lists

proj_cnt = numel(source_list_1);

if numel(source_list_1) ~= numel(source_list_2)
    fprintf('[Fatal Error] numel(source_list_1) %d ~= numel(source_list_2) %d \n', numel(source_list_1), numel(source_list_2));
end

source_list_new = cell(1, proj_cnt);
for i = 1:proj_cnt
    img1 = source_list_1{i};
    img2 = source_list_2{i};
    
    s1 = max( size(img1, 1), size(img2, 1) );
    s2 = max( size(img1, 2), size(img2, 2) );
    
    img1 = my_reshape_img(img1, [s1, s2]);
    img2 = my_reshape_img(img2, [s1, s2]);
    
    source_list_new{i} = img1 + img2;
end


end

