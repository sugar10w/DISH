

%% 参数初始化 Parameter initialization
main_initialize;

%% 生成通用后缀 Generate common suffix
str_notes_suffix = sprintf('prj%d_n%.2f_%s_%s_%s_%s_max%.1f', proj_cnt, n_ri, flag_iter_method, flag_init_method, flag_pb_method, flag_pf_method, proj_intensity);

if proj_latitude_0<1e-5  % 垂直投影 perpencular projection 
    str_notes_suffix = sprintf('%s_ppdcl', str_notes_suffix); 
end

if strcmp(flag_pf_method, 'int')  % 带衰减的仿真 attenuated intensity 
    str_notes_suffix = sprintf('%s_i%.2f_p%.1fum', str_notes_suffix, intensity_attenuation, sample_pixelsize*1e6); 
end

if flag_intensity_map
    str_notes_suffix = sprintf('%s_I', str_notes_suffix);
end

if strcmp(flag_iter_method, 'L2') 
    str_notes_suffix = sprintf('%s_L2', str_notes_suffix);
end

fprintf('[Info] `main2_n` start time: %s\n', datetime);


%% 生成投影矩阵proj_mat和投影函数; Generate projection matrix proj_mat and projection function
proj_mat = my_get_proj_mat_n(proj_cnt, n_ri, proj_latitude_0, flag_rotate45);  
my_project_backward_fun = @(the_sample)my_pb__frame(the_sample, sample_pixelsize, proj_mat, pixelsize, z_pixelsize, my_pb_fun, NumWorkers=pb_num_workers);
my_project_forward_fun = @(source_list)my_pf__frame(source_list, proj_mat, pixelsize, z_pixelsize, sample_n, sample_pixelsize, my_pf_fun, NumWorkers=pf_num_workers, IntensityMap=the_intensity_map_list, AttenuationMask=attenuation_mask);
my_project_backward_fun_allownegative = @(the_sample)my_pb__frame(the_sample, sample_pixelsize, proj_mat, pixelsize, z_pixelsize, my_pb_fun, NumWorkers=pb_num_workers, eps=-1);

%% 初始化source_list; initialize source_list
if exist('source_list_init_filename', 'var')
    load(source_list_init_filename);
    if ~exist('source_list', 'var')
        fprintf(2, '[Error][Init] `source_list` is not found in `%s`\n', source_list_init_filename);
        return;
    end
    if numel(source_list)~=proj_cnt
        fprintf(2, '[Error][Init] `source_list` contains %d elements but `proj_cnt`=%d \n', numel(source_list), proj_cnt);
        return;
    end
    if ~exist('init_ratio', 'var')
        init_ratio = 1;
    end
    if init_ratio ~= 1
        source_list = cellfun(@(x) imresize(x, init_ratio), source_list, "UniformOutput", false);
    end
    source_list = my_postprocess_source_list_fun(source_list);
end

% 对source_list进行预处理
if exist('source_list', 'var')
    fprintf('[Info][Init] `source_list` has been loaded, so `flag_init_method`(%s) would not works.\n', flag_init_method);
else
    source_list = cell(1, proj_cnt); 
    for i = 1:proj_cnt
        source_list{i} = zeros(source_n, 'single');
    end
    the_image = zeros([sample_n, sample_n, sample_nz], 'single');
end
source_list = my_postprocess_source_list_fun(source_list);

for i = 1:proj_cnt
    source_list{i} (source_list{i}>proj_intensity) = proj_intensity;
end

% 生成the_image
if ~exist('the_image', 'var')
    the_image = my_project_forward_fun(source_list);
    
    [Jaccard_index, img2_th_best] = my_Jaccard_index_get_th(the_sample, the_image, thresh_cof, my_Jaccard_index_fun);
    fprintf('[info][easy] jaccard=%d, th=%d\n', Jaccard_index, img2_th_best);
end