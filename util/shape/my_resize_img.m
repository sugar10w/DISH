function [img_new] = my_resize_img(img, ratio)
%MY_RESIZE_IMG  图像缩放 Image rescaling
%  这个函数严格地保证中心位置无偏移（且边长奇偶性不变）
%  This function strictly guarantees that the center is not offset (and the parity of the side lengths is unchanged)

ratio = abs(ratio);

% 变换矩阵 transformation matrix
mat_transform = [ ratio, 0; 0, ratio ]; 
proj_mat3  = [mat_transform, [0;0]; 0,0,1 ];
tform = affine2d(proj_mat3);

[ny, nx] = size(img);  
the_parity = mod(nx,2); 

RA = imref2d( [ny,nx], [-nx/2,nx/2], [-ny/2,ny/2] );


nx_new = ceil(nx * ratio);
ny_new = ceil(ny * ratio);
if mod(nx_new, 2)~=the_parity
    nx_new = nx_new + 1;
end
if mod(ny_new, 2)~=the_parity
    ny_new = ny_new + 1;
end
RB = imref2d( [ny_new,nx_new], [-nx_new/2,nx_new/2], [-ny_new/2,ny_new/2] );


[img_new, ~] = imwarp(img, RA, tform, 'OutputView', RB);

end

