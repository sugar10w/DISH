function temp_mat = my_pb__rot_resize(the_sample, proj_mat, ...
    sample_pixelsize, pixelsize, z_pixelsize, ...
    temp_z_n)
%MY_PB__ROT_RESIZE 在反向传播(三维转二维)时用到的小函数
%  used in back projection (3D -> 2D)

if nargin<6
    temp_z_n = floor( size(the_sample,3) * sample_pixelsize / z_pixelsize ); % z轴分辨率和psf一致 The z-axis resolution is the same as PSF
end

% 计算在新分辨率下的尺寸，然后重新调整尺寸 Calculate size at new resolution, then resize
temp_n_1 = floor( size(the_sample,1) * sample_pixelsize / pixelsize );
temp_n_2 = floor( size(the_sample,2) * sample_pixelsize / pixelsize );

the_parity = mod(size(the_sample, 1),2); 
if mod(temp_n_1, 2) ~= the_parity
    temp_n_1 = temp_n_1 + 1;
end
if mod(temp_n_2, 2) ~= the_parity
    temp_n_2 = temp_n_2 + 1;
end
% if mod(temp_z_n, 2) ~= the_parity
%     temp_z_n = temp_z_n + 1;
% end

% RB2 = imref3d([temp_n_1,temp_n_2,temp_z_n], [-temp_n_1/2,temp_n_1/2], ...
%     [-temp_n_2/2,temp_n_2/2], [-temp_z_n/2,temp_z_n/2]);
% if the_parity==1
%     RB2 = imref3d([temp_n_1,temp_n_2,temp_z_n], [-temp_n_1/2,temp_n_1/2], ...
%         [-temp_n_2/2,temp_n_2/2], [-temp_z_n/2,temp_z_n/2]);
% else
%     RB2 = imref3d([temp_n_1,temp_n_2,temp_z_n], [-temp_n_1/2,temp_n_1/2]-0.5, ...
%         [-temp_n_2/2,temp_n_2/2]-0.5, [-temp_z_n/2,temp_z_n/2]-0.5);
% end

% 对样本进行旋转 rotate the sample
the_sample_rotate = my_rotate_2(the_sample, proj_mat, 'crop');
temp_mat = imresize3(the_sample_rotate, [temp_n_1, temp_n_2, temp_z_n]);

% transformation matrix
% the_matrix = proj_mat * [ ...
%     sample_pixelsize/pixelsize 0 0; ...
%     0 sample_pixelsize/pixelsize 0; ...
%     0 0 sample_pixelsize/z_pixelsize ];
% temp_mat = my_rotate_2(the_sample, the_matrix, 'RB2', RB2);

end

