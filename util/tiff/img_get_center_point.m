function [center_point] = img_get_center_point(img)
% 获得图像的重心 get the center of gravity of the image

center_point = [0,0];

img_a = sum(img, 1);
center_point(1) = sum( img_a .* (1:numel(img_a)) ) / sum(img_a);

img_b = sum(img, 2)';
center_point(2) = sum( img_b .* (1:numel(img_b)) ) / sum(img_b);

end

