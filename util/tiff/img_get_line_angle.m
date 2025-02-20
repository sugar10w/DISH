function [ the_best_angle ] = img_get_line_angle(img, angle_list)
% �ҵ�ͼ����������ֱ�ߵĽǶ� Find the angle of the brightest line in the image

if nargin<2
    angle_list = 0:179;
end

% ���¶���ͼ���С redefine image size
if size(img,1)>400
    img = imresize(img, 400/size(img,1));
end

width_list = zeros( 1, numel(angle_list) );

for i = 1:numel(angle_list)
    angle = angle_list(i); 
    img_rot = imrotate(img, -angle, 'crop');
    img_a = sum(img_rot, 2);
    width_list(i) = sum(img_a > max(img_a)*0.1);
end

min_width_list = angle_list( width_list==min(width_list) );
if max(min_width_list)-min(min_width_list)>90
    min_width_list(min_width_list<90) = min_width_list(min_width_list<90)+180;
end
the_best_angle = mean( min_width_list );   %�����С����Щ�Ƕȵ�ƽ��ֵ Average of those angles with the smallest peak width

end

