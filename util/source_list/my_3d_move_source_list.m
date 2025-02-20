function [source_list_new] = my_3d_move_source_list(source_list, moving_center, flag_rotate45)
% Ĭ����б�Ƕ���45�ȣ�����ֻ�ܶ�������source_list_stretchʹ��
% �ƶ������source_list��ʹ��۽����µ���άλ��
%  Transform the input source_list so that it is focused at the new 3D position
% ������ӽǿ� [x y z]�� x�����ң�y�����£�z��ǰ��
%  From the camera perspective, x is left/right, y is up/down, z is front/rear

if nargin<3
    flag_rotate45 = 1;
end


% check moving_center
if ~isequal( size(moving_center), [3,1] )
    moving_center_new = zeros(3,1);
    if numel(moving_center)>=1
        moving_center_new(1) = moving_center(1);
    end
    if numel(moving_center)>=2
        moving_center_new(2) = moving_center(2);
    end
    if numel(moving_center)>=3
        moving_center_new(3) = moving_center(3);
    end
    
    moving_center = moving_center_new;

end

fprintf('[info] my_3d_move_source_list, moving_center y=%.0f, x=%.0f, z=%.0f\n', moving_center(1), moving_center(2), moving_center(3) );


% ��ʼ���� process
proj_cnt = numel(source_list);
proj_mat = my_get_proj_mat(proj_cnt);
source_list_new = cell(1, proj_cnt);

moving_center_ee = moving_center;
moving_center_ee(1) = - moving_center_ee(1);

for proj_i = 1:proj_cnt
    
    img = source_list{proj_i};
    [s1, s2] = size(img);
    
    % ����ƫ������ Calculate offset coordinates
    delta_u = proj_mat{proj_i} * moving_center_ee;
    if flag_rotate45==0
        delta_u1 = round(  ( -delta_u(1)+delta_u(2) ) / sqrt(2)  );
        delta_u2 = round(  ( +delta_u(1)+delta_u(2) ) / sqrt(2)  );
    else
        delta_u1 = -round(delta_u(1));
        delta_u2 = +round(delta_u(2));
    end
    
    fprintf('shift: %d, %d \n', delta_u1, delta_u2);
    
    % �����³ߴ� Calculate new size
    s1_new = s1+2*abs(delta_u1);
    s2_new = s2+2*abs(delta_u2);
    
    img_reshape = my_reshape_img(img, [s1_new, s2_new] );
    img_new     = zeros([s1_new, s2_new], class(img) );
    
    % �����ƶ�ǰ��������ϵ calculate relationship before and after the move
    for x = 1:s1_new
        xx = x + delta_u1;
        if xx<1 || xx>s1_new 
            continue;
        end
        for y = 1:s2_new
            yy = y + delta_u2;
            if yy<1 || yy>s2_new
                continue;
            end
            img_new(x,y) = img_reshape(xx, yy);
        end
    end
    
    source_list_new{proj_i} = img_new;
    
end


end

