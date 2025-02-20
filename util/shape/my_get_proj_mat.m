function proj_mat = my_get_proj_mat(proj_cnt, proj_latitude, flag_rotate45)
% ֱ�ӻ�ȡͶӰ����
%  Get the projection matrix

if nargin<3
    flag_rotate45 = 1;  % Ĭ�ϲ�����45��
end
if nargin<2
    proj_latitude = 0.25 * pi;      % ͶӰγ��  projection latitude
end

proj_longitude = zeros(1, proj_cnt); % ͶӰ����  projection longitude
for i = 1:proj_cnt
    proj_longitude(i) = - (i-1) * 2 / proj_cnt * pi;  % ���ǵ�ϵͳ�������� Determine the sign according to the actual rotation direction
end

proj_mat = cell(1, proj_cnt);
for proj_i = 1:proj_cnt
    proj_x = [sin(proj_longitude(proj_i)), ...
              cos(proj_longitude(proj_i)), ...
              0];
    proj_y = [-cos(proj_longitude(proj_i)) * sin(proj_latitude), ...
              sin(proj_longitude(proj_i)) * sin(proj_latitude), ...
                                             cos(proj_latitude) ];
    proj_z = [cos(proj_longitude(proj_i)) * cos(proj_latitude),  ...
              -sin(proj_longitude(proj_i)) * cos(proj_latitude), ...
                                            sin(proj_latitude)];
    the_proj_mat = [proj_x; proj_y; proj_z;]';
    
    % ��Ҫ��������һ����ת extra spin
    the_rot_mat = [
        cos(proj_longitude(proj_i)), -sin(proj_longitude(proj_i)), 0;
        sin(proj_longitude(proj_i)),  cos(proj_longitude(proj_i)), 0;
        0,                            0,                           1;
        ];
    the_proj_mat = the_proj_mat * the_rot_mat;
    
    % ����45������
    if flag_rotate45 == 0
        the_dmd_angle = -pi/4;
        the_rot45 = [
            cos(the_dmd_angle), -sin(the_dmd_angle), 0;
            sin(the_dmd_angle),  cos(the_dmd_angle), 0;
            0,                   0,                  1; 
            ];
        the_proj_mat = the_proj_mat * the_rot45; 
    end
    
    
    the_proj_mat (abs(the_proj_mat)<1e-10) = 0;  % ����eps; deal with eps
    
    proj_mat{proj_i} = the_proj_mat;
    
end

end