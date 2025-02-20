function [source_list_new] = my_binarize_multiply_source_list_motionblur( source_list, proj_ratio, max_intensity )
% 将source_list中的每一张灰度图拆成多个二值图。考虑运动模糊
%   Split each grayscale image in source_list into multiple binary images

if nargin<3
    max_intensity = my_get_source_list_max(source_list);
end

if nargin<2
    proj_ratio = 3;   % 将proj_cnt扩大proj_ratio倍；magnification
end

proj_cnt = numel(source_list);


%% 重新改变大小
for i = 1:proj_cnt
    source_list{i} = source_list{i} / max_intensity;
end

%% 生成小序列，即从中心向两周展开的序列 Generates small sequences, i.e. sequences spread out from the center to two wings

relative_idx = zeros(1, proj_ratio);
for i = 1:proj_ratio
    relative_idx(i) = (  -mod(i,2)*2+1  ) * ( floor(i/2) );
end

%%
source_list_new = cell(1, proj_cnt*proj_ratio);
source_list_sub = cell(1, proj_ratio);

for proj_i = 1:proj_cnt
    if mod(proj_i, 10)==0
        fprintf('%d / %d\n', proj_i, proj_cnt);
    end

    proj_i_new_c = (proj_i-1)*proj_ratio + 1;

    % 先将图案移动到对称中心，以方便使用imrotate
    [img_moved, v] = move_img_to_center(source_list{proj_i});
    
    % 计算二值图案
    for gray_i = 1:proj_ratio
        the_motion_blur_angle = relative_idx(gray_i) / proj_ratio / proj_cnt * 360;
        img = my_imtranslate( imrotate(img_moved, -the_motion_blur_angle, 'bilinear', 'crop'), v);
        source_list_sub{gray_i} = img > ((gray_i-0.5) / proj_ratio ) ;
    end

    % 移动到指定位置
    for gray_i = 1:proj_ratio
        proj_i_new = mod( proj_i_new_c + relative_idx(gray_i) - 1, proj_cnt*proj_ratio ) + 1;
        source_list_new{ proj_i_new } = source_list_sub{gray_i};
    end
end


end



function [img_moved,v] = move_img_to_center(img)

img_mask = imbinarize(img);

v = [0,0];
[row, col] = find(img_mask == 1);
v(1) = (min(row) + max(row)) /2 - size(img,1)/2;
v(2) = (min(col) + max(col)) /2 - size(img,2)/2;

img_moved = my_imtranslate(img, -v);


end


