function [source_list_new] = my_duplicate_rotate_source_list( source_list, proj_ratio )
%MY_DUPLCATE_SOURCE_LIST 复制source_list，同时进行旋转补偿

if nargin<2
    proj_ratio = 3;   % 将proj_cnt扩大proj_ratio倍；magnification
end

proj_cnt = numel(source_list);


%%
relative_idx = zeros(1, proj_ratio);
for i = 1:proj_ratio
    relative_idx(i) = (  -mod(i,2)*2+1  ) * ( floor(i/2) );
end

%%

source_list_new = cell(1, proj_cnt*proj_ratio);

for proj_i = 1:proj_cnt
    proj_i_new_c = (proj_i-1)*proj_ratio + 1;
    for gray_i = 1:proj_ratio
        proj_i_new = mod( proj_i_new_c + relative_idx(gray_i) - 1, proj_cnt*proj_ratio ) + 1;
        source_list_new{ proj_i_new } = source_list{proj_i};
        
        % 额外旋转补偿
        the_angle = - relative_idx(gray_i) / (proj_cnt*proj_ratio) * 360;
        source_list_new{ proj_i_new } = imrotate(source_list_new{ proj_i_new }, the_angle, 'bilinear', 'crop');
    end
end

end

