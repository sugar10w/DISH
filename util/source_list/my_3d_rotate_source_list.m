function [source_list] = my_3d_rotate_source_list(source_list_raw, the_angle_raw)
% 用一次投影得到的图像，直接产生模型旋转后的系列图像
%  rotate the projected model by shifting source_list

proj_cnt = numel(source_list_raw);
angle_delta = 360 / proj_cnt;

the_shift = round(the_angle_raw / angle_delta);

the_angle = the_shift * angle_delta;

source_list = cell(1, proj_cnt);
for i = 1:proj_cnt
    curr_img = source_list_raw{ mod(i + the_shift -1, proj_cnt) + 1 };
    source_list{i} = imrotate(curr_img, the_angle, 'crop');
end
    
    
end

