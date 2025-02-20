function [img_new] = my_resize_img(img, ratio)
%MY_RESIZE_IMG  ͼ������ Image rescaling
%  ��������ϸ�ر�֤����λ����ƫ�ƣ��ұ߳���ż�Բ��䣩
%  This function strictly guarantees that the center is not offset (and the parity of the side lengths is unchanged)

ratio = abs(ratio);

% �任���� transformation matrix
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

