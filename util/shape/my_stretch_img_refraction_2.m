function [img_new] = my_stretch_img_refraction_2(img, the_angle, lambda, pixelsize, n_ri, theta_r, flag_from_air_to_material)
%MY_STRETCH_IMG_REFRACTION 计算折射前后的复振幅
%   在角谱上进行拉伸变换


theta_i = asin(n_ri * sin(theta_r));

%%
the_angle_e = pi*3/2 - the_angle;
img = my_rotate(img, the_angle_e);

[ny, nx] = size(img);

asp = fft2(fftshift(img));  % 先fftshift再fft2，以避免意外的高频成分
asp = fftshift(asp);   % 把低频成分移动到画面中心
%asp = img;

if mod(ny,2)==0
    coordinate_limits_ny = (lambda/ny/pixelsize) * [ -ny/2 ,  ny/2-1 ];
else
    coordinate_limits_ny = (lambda/ny/pixelsize) * [ -(ny-1)/2 ,  (ny-1)/2 ];
end
if mod(nx,2)==0
    coordinate_limits_nx = (lambda/ny/pixelsize) * [ -nx/2 ,  nx/2-1 ];
else
    coordinate_limits_nx = (lambda/ny/pixelsize) * [ -(nx-1)/2 ,  (nx-1)/2 ];
end


RA = imref2d( [ny,nx], coordinate_limits_ny, coordinate_limits_nx);
RB = RA;

%%

inversefn   = @(c)coordinate_trasform_1(c, n_ri, theta_i, theta_r);
inversefn_2 = @(c)coordinate_trasform_2(c, n_ri, theta_i, theta_r);

%%

if flag_from_air_to_material
    % from air to material
    tform = geometricTransform2d(inversefn, inversefn_2);
    asp_new = imwarp(asp, RA, tform, "linear", 'OutputView', RB);
else
    % from material to air
    tform_2 = geometricTransform2d(inversefn_2, inversefn);
    asp_new = imwarp(asp, RA, tform_2, "linear", 'OutputView', RB);
end

% imagesc(abs(asp));
% figure; imagesc(abs(asp_new));

asp_new = ifftshift(asp_new);
img_new = ifftshift(ifft2(asp_new));
%img_new = asp_new;


img_new = my_rotate(img_new, -the_angle_e);


% intensity correction
the_ratio = cos(theta_r)/cos(theta_i);
if flag_from_air_to_material
    img_new = img_new * the_ratio;
else
    img_new = img_new / the_ratio;
end


end


%%

function img_new = my_rotate(img, the_angle)

mat_rotate =  [ cos(the_angle), sin(the_angle); -sin(the_angle), cos(the_angle);  ];
proj_mat3  = [mat_rotate, [0;0]; 0,0,1 ];
tform = affine2d(proj_mat3);

[ny, nx] = size(img);  
RA = imref2d( [ny,nx], [-nx/2,nx/2], [-ny/2,ny/2] );
RB = RA;

[img_new, ~] = imwarp(img, RA, tform, 'OutputView', RB);

end

%% 

function c_new = coordinate_trasform_1(c, n_ri, theta_i, theta_r)

c_new = c / 1;
c = c / n_ri;

temp_2 = + sqrt(1-c(:,1).^2-c(:,2).^2) * n_ri * sin(theta_r)...
         + n_ri * c(:,2) * cos(theta_r);

temp_a = 1 ^ 2;
temp_b = - 2 * temp_2 * 1 * cos(theta_i);
temp_c = temp_2.^2 - (1-c_new(:,1).^2) * 1^2 * sin(theta_i)^2;

c_new(:, 2) = (-temp_b - sqrt(temp_b.^2 - 4*temp_a.*temp_c))./(2*temp_a);

c_new = c_new * 1;

c_new = real(c_new);

end

%%
function c_new = coordinate_trasform_2(c, n_ri, theta_i, theta_r)

c = c / 1;
c_new = c / n_ri;

temp_1 = + sqrt(1-c(:,1).^2-c(:,2).^2) * 1 * sin(theta_i)...
         + 1 * c(:,2) * cos(theta_i);

temp_a = n_ri ^ 2;
temp_b = - 2 * temp_1 * n_ri * cos(theta_r);
temp_c = temp_1.^2 - (1-c_new(:,1).^2) * n_ri^2 * sin(theta_r)^2;

c_new(:, 2) = (-temp_b - sqrt(temp_b.^2 - 4*temp_a.*temp_c))./(2*temp_a);

c_new = c_new * n_ri;

c_new = real(c_new);

end
