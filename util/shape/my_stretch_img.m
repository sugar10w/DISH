function [img_new] = my_stretch_img(img, ratio, angle)
% 将二维图像按照指定的方向拉伸
% Stretch a 2D image in the specified direction

if nargin==2
    angle = 0;
end

% 先计算旋转矩阵 rotation matrix
mat_rotate =  [ cos(angle), sin(angle); -sin(angle), cos(angle);  ];
% 计算拉伸矩阵 stretch matrix
mat_stretch = [ ratio, 0; 0, 1;  ];  


% 计算变换矩阵 transformation matrix
mat_transform = inv(mat_rotate) * mat_stretch * mat_rotate;
proj_mat3  = [mat_transform, [0;0]; 0,0,1 ];
tform = affine2d(proj_mat3);

[ny, nx] = size(img);  
the_parity = mod(nx,2);   % 这里假定nx,ny的奇偶性是一样的；It is assumed that the parity of nx, and nz is the same;

RA = imref2d( [ny,nx], [-nx/2,nx/2], [-ny/2,ny/2] );

% 计算拉伸完之后的图像尺寸 calculate the size of the image after stretching
if ratio<1
    RB = RA;
else
    % 将角上的四个点进行变换 Transform four points on the corners
    the_sign_tab = [ 1,1;  1,-1;  -1,-1; -1,1  ];
    the_point_raw = [ nx/2, ny/2 ];
    ny_new_max = 0;
    ny_new_min = 0;
    nx_new_max = 0;
    nx_new_min = 0;
    for i = 1:4
        the_point = the_point_raw .* the_sign_tab(i, :);
        the_point_new = mat_transform * the_point';
        nx_new_max = max(nx_new_max, the_point_new(1));
        nx_new_min = min(nx_new_min, the_point_new(1));
        ny_new_max = max(ny_new_max, the_point_new(2));
        ny_new_min = min(ny_new_min, the_point_new(2));
    end
    ny_new = ceil(ny_new_max - ny_new_min);
    nx_new = ceil(nx_new_max - nx_new_min);
    
    if mod(ny_new,2)~=the_parity
        ny_new = ny_new + 1;
    end
    if mod(nx_new,2)~=the_parity
        nx_new = nx_new + 1;
    end
    RB = imref2d( [ny_new,nx_new], [-nx_new/2,nx_new/2], [-ny_new/2,ny_new/2] );
end

[img_new, ~] = imwarp(img, RA, tform, 'OutputView', RB);


end

