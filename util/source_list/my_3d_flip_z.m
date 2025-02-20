function [source_list_new] = my_3d_flip_z(source_list)
% my_3d_flip_z 将每个投影图像进行镜像翻转

proj_cnt = numel(source_list);

source_list_new = cell(1, proj_cnt);

for i = 1:proj_cnt
    the_angle = (i-1)/proj_cnt * 2*pi;
    the_angle_deg = the_angle / pi * 180;

    img = source_list{i};

    img = imrotate(img, the_angle_deg, 'crop'); 
    img = flipud(img);
    img = imrotate(img, -the_angle_deg, 'crop');

    source_list_new{i} = img;
end

end

