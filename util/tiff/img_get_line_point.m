function [ p1, p2, the_best_angle ] = img_get_line_point(img, angle_list)
% 给出图像中最亮的线段的坐标 give the coordinates of the brightest line in the image

% 缺省参数设置 default value
if nargin<2
    angle_list = 0:179;
end

% 如果不是正方形图像，那裁减一下 crop it if neccessary
if size(img,1) ~= size(img,2)
    fprintf('[Error][img_get_line_point] img is not square-shaped \n');
    img_size = min( size(img) );
    img = my_reshape_img(img, [img_size, img_size] );
end

% 
p1 = [0, 0];
p2 = [0, 0];

the_best_angle = img_get_line_angle(img, angle_list);

img_rot   = imrotate(img, -the_best_angle, 'crop');

% 判断线段的位置 determine the position of the line
img_a     = sum(img_rot, 2);
img_a_bw = img_a > max(img_a)*0.1;
img_a_c = mean(find(img_a_bw==1));

% 判断线段的长度 determine the length of the line
img_b = sum(img_rot, 1);   % img_b用来判断线段的长度 determine the length with img_b
img_b_bw = img_b > max(img_b)*0.05;
img_b_1 = find(img_b_bw==1, 1);
img_b_2 = find(img_b_bw==1, 1, 'last');

% 
img_c = (size(img,1)+1)/2; % 图像中心 image center
img_a_c = img_a_c - img_c;
img_b_1 = img_b_1 - img_c;
img_b_2 = img_b_2 - img_c;
the_angle_rad = the_best_angle / 180 * pi;
p1_1 =  sin(the_angle_rad)*img_a_c + cos(the_angle_rad)*img_b_1;  % 根据旋转后推断的线段位置和长度，反推原始位置 calculate original position
p1_2 =  cos(the_angle_rad)*img_a_c - sin(the_angle_rad)*img_b_1; 
p2_1 =  sin(the_angle_rad)*img_a_c + cos(the_angle_rad)*img_b_2;
p2_2 =  cos(the_angle_rad)*img_a_c - sin(the_angle_rad)*img_b_2;
p1_1 = round(p1_1 + img_c);
p1_2 = round(p1_2 + img_c);
p2_1 = round(p2_1 + img_c);
p2_2 = round(p2_2 + img_c);

p1(1) = p1_1;
p1(2) = p1_2;
p2(1) = p2_1;
p2(2) = p2_2;

end