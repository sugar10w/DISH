% 主代码全息优化
%  考虑了折射界面的位置
% main code for Holographic optimization
% The position of the refraction interface is considered

close all; clc;

%% 参数检测和初始化 Parameter detection and initialization
main_initialize;

%% 生成通用后缀 Generate common suffix
str_notes_suffix = sprintf('holo_proj%d_n%.3f_zs%.3f_id%.3f_px%.1f_b%d', proj_cnt, n_ri, z_shift/1e-3, interface_distance/1e-3, pixelsize/1e-6, binarize_ratio);

%% 生成投影矩阵proj_mat和投影函数; Generate projection matrix proj_mat and projection function
proj_mat = my_get_proj_mat_n(proj_cnt, n_ri, proj_latitude_0, flag_rotate45);  

% 使用linear的前传 use my_pf_easy
my_pf_fun_easy = @(the_source)my_pf_easy(the_source, numel(z_H_list));   % 直线传播

my_stretch_img_refraction_fun = @(img, the_angle, flag_from_air_to_material )...
    my_stretch_img_refraction_2(img, the_angle, lambda_0, pixelsize, n_ri, ...
    pi/2-proj_latitude, flag_from_air_to_material);

%% 处理界面问题

z_interface_material = interface_distance / sin(proj_latitude);                     % 在物料中传播的距离
z_interface_air      = (interface_distance / sin(proj_latitude) + z_shift) / n_ri;   % 假定在空气中传播的距离

% 分别生成H矩阵；这两个矩阵应该是很类似的； Generate H matrix respectively; These two matrices should be very similar;
if flag_GPU
    H_list_air      = my_propagation_quick3D_init_GPU(source_n, source_n, -z_interface_air,     pixelsize, lambda_0, 'Fresnel');
    H_list_material = my_propagation_quick3D_init_GPU(source_n, source_n, z_interface_material, pixelsize, lambda,   'Fresnel');
else
    H_list_air      = my_propagation_quick3D_init    (source_n, source_n, -z_interface_air,     pixelsize, lambda_0, 'Fresnel');
    H_list_material = my_propagation_quick3D_init    (source_n, source_n, z_interface_material, pixelsize, lambda,   'Fresnel');
end

%% 根据实际系统参数，修正傅里叶面上的强度
if flag_source_fft_error
    if flag_rotate45 == 1
        fprintf(2, '[Error] flag_source_fft_error only work with flag_rotate45==0\n');
    end

    the_source_fft_error_ratio = sqrt(1000)^(sqrt(2)/2);
    the_source_fft_max = 5; 

    for i = 1:source_n
        for j = 1:source_n
            the_dist = -( (i-ceil(source_n/2)) + (j-ceil(source_n/2))  )/sqrt(2)/source_n * 2;
            the_source_fft_mask(i,j) = the_source_fft_error_ratio ^ the_dist;
        end
    end
    the_source_fft_mask = ifftshift(the_source_fft_mask);

    the_source_fft_mask(the_source_fft_mask>the_source_fft_max) = the_source_fft_max;
    the_source_fft_mask(the_source_fft_mask<1/the_source_fft_max) = 1/the_source_fft_max;
end

%% 生成投影函数  Generate projection function

my_pb_fun_s = @(temp_mat, proj_i_curr_ratio, the_intensity_map, ...
        temp_mat_raw_energy, the_source_new, alpha)...
    my_pb_gs_fres_s9(temp_mat, ...
    H_list, source_n, H_list_air, H_list_material, my_propagation_fun, ...
    proj_i_curr_ratio, stretch_ratio, flag_rotate45, ...
    SourceFftError=the_source_fft_mask, BinarizeRatio=binarize_ratio, ...
    MaxIntensity=proj_intensity, IntensityMap=the_intensity_map, ...
    Margin=source_margin, ...
    StretchFun=my_stretch_img_refraction_fun, ...
    RawEnergy=temp_mat_raw_energy, ...
    InitSource=the_source_new, Alpha=alpha);  % 梯度下降全息优化，考虑场弯曲

my_pf_fun_s = @(the_source, proj_i_curr_ratio)my_pf_gs_fres_s2(the_source, ...
    H_list, z_H_list, H_list_air, H_list_material, ...
    my_propagation_fun, ...
    proj_i_curr_ratio, stretch_ratio, flag_rotate45, ...
    SourceFftError=the_source_fft_mask, ...
    StretchFun=my_stretch_img_refraction_fun);
my_project_forward_fun_s = @(source_list)my_pf__frame_s(source_list, ...
    proj_mat, pixelsize, z_pixelsize, sample_n, sample_pixelsize, ...
    my_pf_fun_s, ...
    AttenuationMask=attenuation_mask, SampleNZ=sample_nz, ...
    NumWorkers=pf_num_workers, IntensityMap=the_intensity_map_list_raw);

%% 处理强度图

if flag_intensity_map
    for i = 1:proj_cnt
        source_list{i} = my_reshape_img(source_list{i}, [source_n, source_n]) .* the_intensity_map_list{i};
    end
end

%% 

% 逐个视角进行优化 optimizing view by view

if ~exist('flag_with_init', 'var')
    flag_with_init = 0;
end
if ~flag_with_init
    source_list_new = cell(1, numel(source_list)); 
end

time_start = datetime;

%%
z_H_list_len = numel(z_H_list);

out_num = 0;

for proj_i = 1:proj_cnt

    proj_i_curr_ratio = (proj_i-1)/proj_cnt;
    
    % preprocess the_sample
    temp_mat_sample = my_pb__rot_resize(the_sample,    proj_mat{proj_i}, sample_pixelsize, pixelsize, z_pixelsize, z_H_list_len);
    %temp_mat_background = temp_mat_sample<0.003;
    temp_mat_background = ~my_select_layered_details(temp_mat_sample);
    
    temp_mat_sample =     my_reshape_img(temp_mat_sample,     [source_n, source_n]);
    temp_mat_background = ~my_reshape_img(~temp_mat_background,[source_n, source_n]);
    
    temp_mat = my_pf_fun_easy( source_list{proj_i} );
    temp_mat = my_reshape_img(temp_mat, [source_n, source_n]);
    temp_mat(temp_mat<1e-9) = 1e-9;
    
    temp_mat_raw_energy = temp_mat;
    temp_mat_raw_energy(temp_mat_sample<0.003) = 0;

    temp_mat(temp_mat_background) = 0;
    
    % iterative optimize
    source_list_new{proj_i} = my_pb_fun_s(temp_mat, proj_i_curr_ratio, the_intensity_map_list_raw{proj_i}, temp_mat_raw_energy, ...
        source_list_new{proj_i}, 0.12 ); %0.010 

    if ispc && mod(proj_i,1)==0 % estimate time
        %fprintf(char(ones(1, out_num) * 8));
        out_num = fprintf('[info][pb_s] proj %d/%d finished, %s\n', proj_i, proj_cnt, my_predict_time(time_start, proj_i/proj_cnt) );
    end
end



%%
save( sprintf('%s/%s_fine_%s.mat', the_dir_savename, sample_id, str_notes_suffix ),  'source_list_new', 'flag_rotate45');  

if binarize_ratio>0
    source_list_new = my_binarize_multiply_source_list(source_list_new, binarize_ratio, proj_intensity);
    save( sprintf('%s/%s_fine_b_%s.mat', the_dir_savename, sample_id, str_notes_suffix ),  'source_list_new', 'flag_rotate45');

    proj_cnt = numel(source_list_new);
    %proj_mat = my_get_proj_mat_n(proj_cnt, n_ri, proj_latitude_0, flag_rotate45);
    my_project_forward_fun_s = @(source_list)my_pf__frame_s(source_list, ...
        proj_mat, pixelsize, z_pixelsize, sample_n, sample_pixelsize, my_pf_fun_s, ...
        NumWorkers=pf_num_workers, IntensityMap=the_intensity_map_list, ...
        SampleNZ=sample_nz, AttenuationMask=attenuation_mask);
end

%% 导出三维图 export 3D volume
if flag_watch_result

    if exist('flag_reference_image', 'var') && binarize_ratio==0 % 对每个视角打分
        if flag_reference_image
            the_reference_image = loadtiff(reference_image_filename);
            the_reference_image = single(the_reference_image);
            the_reference_image = the_reference_image / max(the_reference_image(:));
            the_image = my_pf__frame_s(source_list_new, ...
                proj_mat, pixelsize, z_pixelsize, sample_n, sample_pixelsize, my_pf_fun_s, ...
                NumWorkers=pf_num_workers, IntensityMap=the_intensity_map_list, ...
                SampleNZ=sample_nz, ...
                ReferenceImage=the_reference_image);
        end
    else
        if binarize_ratio>0
            source_list_new = my_3d_multip(source_list_new, proj_intensity);
        end
        the_image = my_project_forward_fun_s(source_list_new);
    end
    
    if numel(the_image)>2e8
        my_save_tiff_data_type = 8;
    else
        my_save_tiff_data_type = 16;
    end
    my_save_complex_tiff(the_image, sprintf('%s/%s_fine_%s.tif', the_dir_savename, sample_id, str_notes_suffix), '10', my_save_tiff_data_type);
    
end


%% 保存原始图案投影用于对比 Save the original projection for comparison
if flag_watch_result
    if exist('source_list_stretch', 'var')
        
        if exist('flag_reference_image', 'var') && binarize_ratio==0  % 对每个视角打分
            if flag_reference_image
                the_reference_image = loadtiff(reference_image_filename);
                the_reference_image = single(the_reference_image);
                the_reference_image = the_reference_image / max(the_reference_image(:));
                fprintf('\n');
                the_image = my_pf__frame_s(source_list_stretch, ...
                    proj_mat, pixelsize, z_pixelsize, sample_n, sample_pixelsize, my_pf_fun_s, ...
                    NumWorkers=pf_num_workers, IntensityMap=the_intensity_map_list, ...
                    SampleNZ=sample_nz,...
                    ReferenceImage=the_reference_image);
            end
        else
            source_list_stretch = my_postprocess_source_list_fun(source_list_stretch);
            if binarize_ratio>0
                source_list_stretch = my_binarize_multiply_source_list(source_list_stretch, binarize_ratio, proj_intensity);
            end
            the_image = my_project_forward_fun_s(source_list_stretch);
        end

        if numel(the_image)>2e8
            my_save_tiff_data_type = 8;
        else
            my_save_tiff_data_type = 16;
        end
        my_save_complex_tiff(the_image, sprintf('%s/%s_before_%s.tif', the_dir_savename, sample_id, str_notes_suffix), '10', my_save_tiff_data_type);
    end
end

