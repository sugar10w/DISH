function img_proj_list = my_getsideviews(img, angle_list)
% ���ָ�������в����ӽǵ�ͶӰͼ
% Get projections of all specified side views

img = double(img);

img_proj_list = zeros( size(img,2), size(img,3), numel(angle_list)  );

for i = 1:numel(angle_list)
    
    img_temp = imrotate3(img, angle_list(i), [0,0,1], 'crop' );
    
    % �����ӽ�
    img_temp_sideview = sum(img_temp, 1);
    img_temp_sideview = reshape(img_temp_sideview, [size(img_temp,2), size(img_temp,3) ]);
    
    img_proj_list(:,:,i) = img_temp_sideview;
    
end


end