function proj_mat = my_get_proj_mat_n(proj_cnt, n_ri, proj_latitude_0, flag_rotate45)
% 根据折射率，获取投影矩阵
%  Get the projection matrix according to the refractive index


if nargin<4
    flag_rotate45 = 1;  % 默认不考虑45度
end
if nargin<3
    proj_latitude_0 = 0.25 * pi;  % 默认入射角度 projection latitude
end
if nargin<2
    n_ri = 1.5;    % 默认折射率 refractive index
end

% 计算折射后的角度 Calculate the angle after refraction
proj_latitude = pi/2 - asin(  sin(pi/2 - proj_latitude_0) * 1 / n_ri ); 

proj_mat = my_get_proj_mat(proj_cnt, proj_latitude, flag_rotate45);

end