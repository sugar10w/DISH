function proj_mat = my_get_proj_mat_n(proj_cnt, n_ri, proj_latitude_0, flag_rotate45)
% ���������ʣ���ȡͶӰ����
%  Get the projection matrix according to the refractive index


if nargin<4
    flag_rotate45 = 1;  % Ĭ�ϲ�����45��
end
if nargin<3
    proj_latitude_0 = 0.25 * pi;  % Ĭ������Ƕ� projection latitude
end
if nargin<2
    n_ri = 1.5;    % Ĭ�������� refractive index
end

% ���������ĽǶ� Calculate the angle after refraction
proj_latitude = pi/2 - asin(  sin(pi/2 - proj_latitude_0) * 1 / n_ri ); 

proj_mat = my_get_proj_mat(proj_cnt, proj_latitude, flag_rotate45);

end