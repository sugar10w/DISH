function the_source_air = my_pb_gs_fres_s9(temp_mat, ...
    H_list, psf_n, H_list_air, H_list_material, my_propagation_fun, ...
    proj_i_curr_ratio, stretch_ratio, flag_rotate45, ...
    varargin)
%MY_PB_GS_FRES 
% 梯度下降方法；输入输出(temp_mat, the_source)都要求是能量而不是复数场；这个函数不能用于计算下降方向。
%  input and output (temp_mat, the_source) are both required to be energy/intensity fields rather than complex fields;
%  This function cannot be used to calculate descent direction 
% 考虑折射过程中角谱的扭曲
%  Consider the distortion of angular spectra during refraction 

%%

flag_source_fft_error = 0;
the_source_fft_mask = ones(size(H_list{1}),'single');

flag_binarize = 0;
binarize_ratio = 0;

the_source_intensity_max = max(temp_mat(:));

flag_intensity_map = 0;
the_intensity_map = NaN;

source_margin = 0;

my_stretch_img_refraction_fun = @(the_img, the_angle, flag_from_air_to_material) ...
    my_stretch_img_easy(the_img, the_angle, flag_from_air_to_material, stretch_ratio);

temp_mat_raw_energy = NaN;

alpha = 0.050;
iter_max = 20;

flag_with_init = 0;  % 是否以有初始值的情况进行优化
the_source_material = ones(psf_n, psf_n,'single');
the_source_air = ones(psf_n, psf_n,'single');      % 首先给出初始值，选择全白 initialization

for i = 1:numel(varargin)-1
    if ~isstring(varargin{i}) && ~ischar(varargin{i})
        continue;
    end
    switch varargin{i}
        case 'SourceFftError'
            flag_source_fft_error = 1;
            temp = varargin{i+1};
            the_source_fft_mask = temp;
        case 'BinarizeRatio'
            temp = varargin{i+1};
            if isnumeric(temp) && temp>0
                flag_binarize = 1;
                binarize_ratio = temp;
            end
        case 'MaxIntensity'
            temp = varargin{i+1};
            if isnumeric(temp) && temp>0
                the_source_intensity_max = temp;
            end
        case 'IntensityMap'
            temp = varargin{i+1};
            if ~isnumeric(temp) || ~ismatrix(temp)
                if ~isnan(temp)
                    fprintf(2,'[Error] IntensityMap must be an image.\n');
                end
            else
                flag_intensity_map = 1;
                the_intensity_map = temp;
            end
        case 'Margin'
            temp = varargin{i+1};
            if isnumeric(temp) && temp>0
                source_margin = temp;
            end
        case 'StretchFun'
            temp = varargin{i+1};
            if isa(temp, 'function_handle')
                my_stretch_img_refraction_fun = temp;
            end
        case 'RawEnergy'
            temp = varargin{i+1};
            if isnumeric(temp) && ndims(temp)==3
                temp_mat_raw_energy = temp;
            else
                fprintf(2,'[Error] RawEnergy must be a 3d matrix.\n');
            end
        case 'Alpha'
            temp = varargin{i+1};
            if isnumeric(temp) && temp>0
                alpha = temp;
            end
        case 'MaxIter'
            temp = varargin{i+1};
            if isnumeric(temp) && temp>0
                iter_max = temp;
            end
        case 'InitSource'
            temp = varargin{i+1};
            if isnumeric(temp) && numel(temp)>100
                the_source_air = temp;
                flag_with_init = 1;
                the_source_air = my_reshape_img(the_source_air, [psf_n, psf_n]);
            end
    end
end

%% the_intensity_map

if isnan(the_intensity_map)
    flag_intensity_map = 0;
end
if flag_intensity_map
    the_intensity_map(the_intensity_map<0) = 0;
    the_intensity_map = the_intensity_map .^ 0.5;
end

%%

psf_z_n = numel(H_list);
temp_z_n = size(temp_mat, 3);

% 修改temp_mat的三维尺寸  Modify the 3D size of temp_mat
temp_mat = single(temp_mat);
temp_mat = my_reshape_img(temp_mat, [psf_n, psf_n] );

if ndims(temp_mat_raw_energy)~=3
    temp_mat_raw_energy = temp_mat;  % 备份
end
temp_mat_raw_energy = single(temp_mat_raw_energy);
temp_mat_raw_energy = my_reshape_img(temp_mat_raw_energy, [psf_n, psf_n] );

% 计算直线传播下的强度谱
temp_mat_proj_energy_2d = max(temp_mat_raw_energy, [], 3);
temp_mat_proj_energy = repmat(temp_mat_proj_energy_2d, [1,1,temp_z_n]);

% 强度最大值
the_source_max = the_source_intensity_max.^0.5;

% 提取边缘区域
temp_mat_ignore = (temp_mat < 1e-10);
for z = 2:temp_z_n-1
    for x = 1:psf_n
        for y = 1:psf_n
            if temp_mat(x,y,z)>0 && temp_mat(x,y,z-1)>0 && temp_mat(x,y,z+1)>0
                temp_mat_ignore(x,y,z) = 1;
            end
        end
    end
end
%temp_mat(temp_mat_ignore) = 0;

temp_mat_shell = ~temp_mat_ignore;

temp_mat_proj_energy_2d_shell = imdilate(temp_mat_proj_energy_2d>1e-10, strel('disk',20));
temp_mat_proj_energy_2d_shell(temp_mat_proj_energy_2d>1e-10) = 0;
z_center = round( size(temp_mat_shell,3)/2 + 0.5);
for i = 1:size(temp_mat_shell,1)
    for j = 1:size(temp_mat_shell,2)
        if temp_mat_proj_energy_2d_shell(i,j)
            temp_mat_shell(i,j,z_center) = 1;
        end
    end
end


% 强度调整 Intensity/energy adjustment
temp_mat(temp_mat<0) = 0;
temp_mat = temp_mat .^ 0.5;   % <!> 将能量转化为幅度  Convert energy/intensity fields to complex fields


%% 开始进行迭代  iteration;

% 对应到H_list中的坐标  indexes in H_list
[z_ll, z_rr] = my_get_reshape_range( temp_z_n, psf_z_n);  

% 拉伸参数
the_angle = pi/2 + proj_i_curr_ratio*2*pi ;
if flag_rotate45 == 0    
    the_angle = the_angle + pi/4;  % DMD45度倾斜放置
end

% 仅在限定的范围内进行计算
temp_mat_proj_energy_2d_s = my_stretch_img(temp_mat_proj_energy_2d, stretch_ratio, the_angle);
the_source_mask = imdilate(temp_mat_proj_energy_2d_s > 1e-5, strel('disk', 50));


% 投影函数
from_air_to_material_fun = @(the_source_air)from_air_to_material(the_source_air, ...
    my_propagation_fun, H_list_air, H_list_material, ...
    the_angle, ...
    my_stretch_img_refraction_fun, ...
    the_intensity_map, ...
    flag_source_fft_error, the_source_fft_mask);

from_material_to_air_fun = @(the_source_add)from_material_to_air(the_source_add, ...
    my_propagation_fun, H_list_air, H_list_material, ...
    the_angle, ...
    my_stretch_img_refraction_fun, ...
    the_intensity_map, ...
    flag_source_fft_error, the_source_fft_mask);

%%
% 最优值记录
the_best_MSE_auto = inf;
the_best_MSE_source = the_source_material;

the_last_target_loss = inf;
cnt_ok = 0;
cnt_ok_trigger_num = 4;

%%
time_start = datetime;
out_num = 0;

for iter = 1:iter_max
    
    %imagesc(the_source_air)
    
    % 从the_source_air，计算 the_source_material
    if flag_with_init || iter>1
        the_source_material = from_air_to_material_fun(the_source_air);
    end
    
    if binarize_ratio>0 && max(abs(the_source_air(:)).^2) < the_source_max^2/binarize_ratio
        the_source_air = the_source_air*2;
        fprintf(2, '[Warning] the intensity is abnormal! `proj_intensity` may be too large or `binarize_ratio` may be too small.\n');
    end

    % 传播获得当前情况 current the_image
    if ~flag_binarize || (iter==1 && ~flag_with_init)
        the_image = my_propagation_fun(the_source_material, H_list, z_ll:z_rr);
    else
        the_image_1 = my_propagation_fun(the_source_material, H_list, z_ll:z_rr);
        the_image_b_energy = zeros(size(the_image_1), 'single');
        the_source_list = binarize_multiply_source(the_source_air, the_source_max, binarize_ratio);
        for i = 1:binarize_ratio
            the_source_material_b = from_air_to_material_fun(the_source_list{i});
            the_image_b = my_propagation_fun(the_source_material_b, H_list, z_ll:z_rr);
            the_image_b_energy = the_image_b_energy + abs(the_image_b).^2;
        end
        the_image = the_image_b_energy.^0.5 .*exp(1i*angle(the_image_1));
    end

    % 选择MSE最高
    the_image_temp = abs(the_image).^2;

    temp_mat_proj_energy_temp = temp_mat_proj_energy(temp_mat_shell);
    temp_mat_proj_energy_temp = temp_mat_proj_energy_temp / max(temp_mat_proj_energy_temp(:));
    the_MSE_auto = my_MSE_auto(temp_mat_proj_energy_temp, the_image_temp(temp_mat_shell));
    
    fprintf('gs %d, MSE_auto = %d ; \n', iter,  the_MSE_auto );
    
    if iter>1 || flag_with_init
        
        if the_MSE_auto > the_last_target_loss
            alpha = alpha / 2;
            cnt_ok = 0;
            cnt_ok_trigger_num = cnt_ok_trigger_num + 1;
            %fprintf(2,'failed to update at gs %d \n', iter);
        else
            if the_best_MSE_auto > the_MSE_auto
                the_best_MSE_source = the_source_air;
                the_best_MSE_auto = the_MSE_auto;

                cnt_ok = cnt_ok + 1;
                if cnt_ok>=cnt_ok_trigger_num
                    alpha = alpha*1.2;
                    cnt_ok = 0;
                end
            else
                cnt_ok = 0;
                %fprintf(2,'failed to update at gs %d \n', iter);
            end
            
        end
        the_last_target_loss = the_MSE_auto;

    end

    % 更新相位 update phase
    if (iter==1 && ~flag_with_init)
        temp_mat_ph = temp_mat .*  exp(1i*angle(the_image));
    
        % 回传 update the_source
        the_source_add = zeros(psf_n, 'single');
        for i = 1:temp_z_n
            the_source_sub = my_propagation_fun(temp_mat_ph(:,:,i), H_list, z_ll+i-1, -1);
            the_source_add = the_source_add + the_source_sub;   % 这里选择的策略是直接叠加 direct accumulating
        end
    
        % 回传
        the_source_air_e = from_material_to_air_fun(the_source_add);

        % 灭相位 Convert to real number
        [the_source_air, ~] = my_complex_to_positive_real(the_source_air_e);
    else

        %% 尝试寻找下降方向
        temp_mat_curr = abs(the_image).^2;
        temp_mat_delta = temp_mat_proj_energy - temp_mat_curr;
        
        %temp_mat_delta = temp_mat_delta .* the_image;
        temp_mat_delta = temp_mat_delta .* exp(1j*angle(the_image));
        
        %
        temp_mat_delta(temp_mat_shell==0) = 0;

        %
        temp_mat_delta(isnan(temp_mat_delta)) = 0;

        % 回传
        the_source_add = zeros(psf_n, 'single');
        for i = 1:temp_z_n
            the_source_sub = my_propagation_fun(temp_mat_delta(:,:,i), H_list, z_ll+i-1, -1);
            the_source_add = the_source_add + the_source_sub;   % 这里选择的策略是直接叠加 direct accumulating
        end
        the_source_air_delta = from_material_to_air_fun(the_source_add);

        fprintf('alpha=%d, norm=%d\n', alpha, norm(the_source_air_delta));

        % 尝试更新
        the_source_air_delta = the_source_air_delta ./ norm(the_source_air_delta);
        the_source_air_delta(isnan(the_source_air_delta)) = 0;

        the_source_air_new = the_source_air + norm(the_source_air)*alpha*the_source_air_delta; 
        %alpha = alpha * 0.8;
        
        [the_source_air, ~] = my_complex_to_positive_real(the_source_air_new);

    end

    % 限制幅度
    the_source_air(the_source_air>the_source_max) = the_source_max;
    the_source_air = my_clean_margin(the_source_air, source_margin, 0);
    the_source_air(the_source_mask==0) = 0;

    %imagesc(the_source_air)
    %fprintf(char(ones(1, out_num) * 8));
    %out_num = fprintf('[info][gs] proj %d/%d finished, error=%d, %s\n', iter, iter_max, the_best_MSE_auto, my_predict_time(time_start, iter/iter_max) );

end
%fprintf(char(ones(1, out_num) * 8));

%%
the_source_air = the_best_MSE_source; % 替换成表现最好的那个

the_source_air = the_source_air.^2; % 输出的时候记得平方转化为能量 convert to energy/intensity fields 

end



% 从空气中的共轭面到容器内
function the_source_material = from_air_to_material(the_source_air, ...
    my_propagation_fun, H_list_air, H_list_material, ...
    the_angle, ...
    my_stretch_img_refraction_fun, ...
    the_intensity_map, ...
    flag_source_fft_error, the_source_fft_mask)
if ~isnan(the_intensity_map)
    the_source_air = the_source_air.*the_intensity_map;
end
if flag_source_fft_error
    the_source_air = ifft2( fft2(the_source_air) .* the_source_fft_mask );
end
flag_from_air_to_material = 1;
the_interface_air      = my_propagation_fun(the_source_air, H_list_air); 
the_interface_material = my_stretch_img_refraction_fun(the_interface_air,      the_angle, flag_from_air_to_material ); %the_interface_material 物料侧角谱
the_interface_material = my_reshape_img(the_interface_material, size(the_source_air));
the_source_material    = my_propagation_fun(the_interface_material, H_list_material); %the_source_material 物料内假想的共轭面处的场
end

% 从容器内到空气中
function the_source_air_e = from_material_to_air(the_source_add, ...
    my_propagation_fun, H_list_air, H_list_material, ...
    the_angle, ...
    my_stretch_img_refraction_fun, ...
    the_intensity_map, ...
    flag_source_fft_error, the_source_fft_mask)
flag_from_air_to_material = 0;
the_interface_material = my_propagation_fun(the_source_add, H_list_material, 1, -1);  % 界面处，物料侧角谱
the_interface_air      = my_stretch_img_refraction_fun(the_interface_material, the_angle, flag_from_air_to_material );   % 界面处，物料侧角谱
the_interface_air      = my_reshape_img(the_interface_air, size(the_source_add));
the_source_air_e       = my_propagation_fun(the_interface_air, H_list_air, 1, -1);  % 假想在空气中的DMD面
if flag_source_fft_error
    the_source_air_e = ifft2( fft2(the_source_air_e) .* the_source_fft_mask );
end
if ~isnan(the_intensity_map)
    the_source_air_e = the_source_air_e ./ the_intensity_map;
    the_source_air_e(the_intensity_map<0.2) = 0;
end

end

% 二值化
function the_source_list = binarize_multiply_source(the_source, the_max, the_ratio)
% 换回能量
the_source_energy = abs(the_source).^2;
the_max_energy = the_max ^2;
% 二值化
the_source_list = cell(1, the_ratio);
for i = 1:the_ratio
    img_binary = (the_source_energy > the_max_energy * (i-0.5) / the_ratio);
    img_binary = (the_max_energy/the_ratio)^0.5 * img_binary;

    the_source_list{i} = img_binary;
end

end

