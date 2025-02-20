% ������ʼ�� Parameter initialization


%%
addpath('./util/');
addpath('./util/algorithm/');    % �����㷨 iterative algorithms
addpath('./util/project/');      % ͶӰ projection methods
addpath('./util/propagation/');  % ��ѧ���� optical propagation
addpath('./util/shape/');  % ��ת������ߴ�� rotate, resize, etc
addpath('./util/tiff/');   % tiff�ļ���أ�ͼ���ȡ tiff file related
addpath('./util/time/');   % ʱ�� time
addpath('./util/source_list/');  % ͼ������ source list (DMD image sequence)

addpath('./extensions/ParforProgMon/');  % Parfor progress

m = 1;
cm = 1e-2;
mm = 1e-3;
um = 1e-6;
nm = 1e-9;

%% �ļ� file

if ~exist('sample_filename_mat', 'var')
    sample_filename_mat = '';  % ��ȡtif�ļ� read from .tif file
    %sample_filename_mat = 'selected_points';   % ��mat�ļ���ȡ read from .mat file
end

if ~exist('sample_filename_mat', 'var') && ~exist('sample_filename', 'var') && ~exist('pf_only_sample_size', 'var')
    fprintf(2, '[Warning] `sample_filename` and `sample_size` not found! use default filename\n');
    
    sample_filename = 'model/PointsTest/points p-4.tif';
end


%% �������� basic setting
if ~exist('flag_GPU', 'var')
    flag_GPU = 0;          % �Ƿ�ʹ��GPU����FFT Whether to use GPU to calculate FFT
end

if ~exist('flag_rotate45', 'var')
    flag_rotate45 = 1;    % 45����б���õ�DMD
end

if ~exist('pb_num_workers', 'var')
    pb_num_workers = NaN;
end
if ~exist('pf_num_workers', 'var')
    pf_num_workers = NaN;
end

if ~exist('flag_init_method', 'var')
    %flag_init_method = 'zero';   % ȫ���ʼ�� initialize with zeros
    %flag_init_method = 'random'; % �������ʼ�� initialize with random number
    %flag_init_method = 'binary'; % �ö�ֵͶӰ��ʼ�� initialize by binary projection
    flag_init_method = 'auto';   % auto��flag_pb_method�趨��ͶӰ initialize by project_backward which is defined by `flag_pb_method`
end

if ~exist('flag_pb_method', 'var')
    %flag_pb_method = 'gs';     % gs�����㷨���������ڼ����½����� G-S iterative algorithm. This function cannot be used to calculate descent direction
    %flag_pb_method = 'gs0';    % gs�����㷨��Ҫ�󲻽��������������������ڼ����½�����G-S iterative algorithm without energy correction. This function cannot be used to calculate descent direction
    flag_pb_method = 'mean';   % ȡ��ƽ��ֵ get the projected mean
    %flag_pb_method = 'sum';    % ȡ��ͶӰ�� get the projected sum
    %flag_pb_method = 'max';    % ȡ�����ֵ���������ڼ����½����� Get the projected maximum.  This function cannot be used to calculate descent direction
    %flag_pb_method = 'binary'; % ��ֵͶӰ���������ڼ����½����� Get the binarized projection result. This function cannot be used to calculate descent direction
end

if ~exist('flag_pf_method', 'var')
    %flag_pf_method = 'gs';     % fresnel���� fresnel propagation
    flag_pf_method = 'easy';   % ֱ�ߴ��� linear propagation
    %flag_pf_method = 'incoherent';   % ����ɴ��� incoherent propagation
    %flag_pf_method = 'int';          % ��ǿ��˥����ֱ�ߴ��� linear propagation with intensity attenuation
end

%% ����ģ��ϸ�ڵ����� reinforce model details

if ~exist('flag_reinforcement', 'var')  
    flag_reinforcement = 0;
end
if ~exist('reinf_dist', 'var')
    reinf_dist = 150e-6;
end
if ~exist('reinf_max_weight', 'var')
    reinf_max_weight = 2;
end
if ~exist('reinf_down_rate', 'var')
    reinf_down_rate = 1;
end

%% ͶӰ�ӽ�����ͶӰ���ǿ�� number of projection angles, and projection maximum intensity

if ~exist('proj_cnt', 'var')
    proj_cnt = 180;        % ͶӰ�ӽ��� number of projection angles
end
if ~exist('proj_intensity', 'var')
    proj_intensity = 100;  % ͶӰ���ǿ�� projection maximum intensity
end

if ~exist('binarize_ratio','var')
    binarize_ratio = 0;  % ��ֵ������
end

%% �������� iteration parameters

if ~exist('flag_iter_method', 'var')
    flag_iter_method = 'none';  % ������ Do not iterate
    %flag_iter_method = 'L2';    % L2(MSE)��С������  L2(MSE) Minimization Iteration
    %flag_iter_method = 'RL';    % RL����  Richard-Lucy iteration
end

if ~exist('iter_max', 'var')
    iter_max = 0;  % ���������� maximum number of iteration rounds
end
if ~exist('flag_watch_iter', 'var')
    flag_watch_iter = 1;   % �Ƿ񱣴�������� Whether to save the iterative process
end

if ~exist('flag_watch_result', 'var')
    flag_watch_result = 0;   % �Ƿ񵼳���άͼ whether to export final 3D volume
end

%% ����ģ�� load model

if strcmp(sample_filename_mat, '')==1 && ~exist('sample_filename', 'var') && exist('pf_only_sample_size', 'var') 
    % ������ģ�ͣ�����ȡģ�ͳߴ�pf_only_sample_size
    sample_n = pf_only_sample_size; 

    if ~exist('pf_only_sample_size_z', 'var')
        pf_only_sample_size_z = sample_n;
    end
    sample_nz = pf_only_sample_size_z;

    if ~exist('sample_id', 'var')
        sample_id = 'pf_only';
    end
    
elseif strcmp(sample_filename_mat, '')==1 && exist('sample_filename', 'var') 
    % ��tiff�ļ���ȡ load from tiff file
    [~, sample_id, ~] = fileparts(sample_filename);
    sample_id(sample_id=='/') = '_';
    sample_id(sample_id=='\') = '_'; 
    
    the_sample = loadtiff(sample_filename);
    the_sample = single(the_sample);
    
    if max(the_sample(:))>100 && max(the_sample(:))<255
        the_sample = the_sample / 255;
        fprintf(2, 'the_sample = the_sample / 255;\n');
    elseif max(the_sample(:))>10000 && max(the_sample(:))<65535
        the_sample = the_sample / 65535;
        fprintf(2, 'the_sample = the_sample / 65535;\n');
    else
        the_sample = the_sample / max(the_sample(:));
    end
    
    sample_n = max( size(the_sample,1), size(the_sample,2) );
    sample_nz = size(the_sample,3);
    sample_nz = max(sample_nz, sample_n);  % sample_nz>=sample_n
    the_sample = my_reshape_img( the_sample, [sample_n, sample_n, sample_nz] ); % ��ǰ����ά�ȳߴ����
    
else
    % ��mat�ļ���ȡ
    sample_id = sample_filename_mat;
    sample_filename = ['./result_samples/', sample_filename_mat, '.mat'];
    if ~exist(sample_filename, 'file')
        fprintf(2, '\n[Error]\n[Error] the sample `%s` cannot be loaded, or `pf_only_sample_size` is not set. \n\n', sample_id);
        return; 
    end

    load(sample_filename);  
    
    if ~exist('sample_id', 'var')
        [~, sample_id, ~] = fileparts(sample_filename_mat);
        sample_id(sample_id=='/') = '_';
        sample_id(sample_id=='\') = '_';
    end
    
    sample_n = max( size(the_sample,1), size(the_sample,2) );
    sample_nz = size(the_sample,3);
    sample_nz = max(sample_nz, sample_n);  % sample_nz>=sample_n
    the_sample = my_reshape_img( the_sample, [sample_n, sample_n, sample_nz] ); % reshape to cube ��ǰ����ά�ȳߴ����

end

if ~exist('sample_pixelsize', 'var')
    if exist('the_ratio', 'var')
        sample_pixelsize = 5.4e-6 * the_ratio;
    else
        sample_pixelsize = 5.4e-6;  % default resolution is 5.4um
    end
end

if ~exist('the_dir_savename', 'var')
    the_dir_savename = ['./result_samples/', sample_id];
end
if ~exist(the_dir_savename, 'dir')
    mkdir(the_dir_savename);
end

%% ����ģ��ϸ�ڵ����� reinforce model details

if flag_reinforcement
    switch flag_reinforcement_method
        case 'erode'
            % the_dir_savename = [the_dir_savename, '_reinf-erode'];
            fprintf('[info] Sample reinforcement calculating... using erode method.\n');
            sam_bin = imbinarize(the_sample);
            reinf_sample = zeros(size(sam_bin));
            reinf_sample(sam_bin) = 1;
            lay_num = floor(reinf_dist/sample_pixelsize);
            now_model = sam_bin;
            for i=1:lay_num
                temp_model = imerode(now_model, strel("sphere", 1));
                reinf_sample(now_model & ~temp_model) = (1 - (i-1)*sample_pixelsize/reinf_dist) * (reinf_max_weight - 1) + 1;
                now_model = temp_model;
            end
        otherwise
            % the_dir_savename = [the_dir_savename, '_reinf-conv'];
            fprintf('[info] Sample reinforcement calculating... using conv method.\n');
            down_model = imresize(the_sample, 1/reinf_down_rate);

            temp_w = zeros(size(down_model));
            all_w = 0;
            lay_num = floor(reinf_dist/sample_pixelsize/reinf_down_rate);
            for i=1:lay_num
                ker = strel("sphere", i+1).Neighborhood;
                ker2 = strel("sphere", i).Neighborhood;
                ker(2:end-1, 2:end-1, 2:end-1) = logical(ker(2:end-1, 2:end-1, 2:end-1) - ker2);

                w_reinf_use = 1/sqrt(i*reinf_down_rate*sample_pixelsize);
                
                temp_w = temp_w + convn(down_model, ker, "same") / sum(ker, "all") * w_reinf_use;
                all_w = all_w + w_reinf_use;
            end
            
            reinf_sample = the_sample .* (reinf_max_weight - imresize(temp_w / all_w, reinf_down_rate) * (reinf_max_weight - 1));
    end
    my_save_complex_tiff(reinf_sample, fullfile(the_dir_savename, 'reinforcement_sample.tif'), '10');
end

%% �������� parameters in propagation

if ~exist('n_ri', 'var')
    n_ri = 1.478 ;  % ������ refractive index
end

if ~exist('z_shift', 'var')
    z_shift = +0.0e-3;   % ����������ھ۽����ĵ�λ�ã����Ŵ������������The position of the conjugate surface, relative to the focus center (measured along the propagation direction)
end

if ~exist('interface_distance', 'var')
    interface_distance = +4.2e-3;  % ���浽�۽����ĵľ��루��ֱ�ڽ��淽������� Distance from the interface to the focus center (measured perpendicular to the interface direction)
end

if ~exist('pixelsize_ex_ratio', 'var')
    pixelsize_ex_ratio = 1;   % ������ϵ�� Oversampling factor
end

%% �ߴ����� size

if ~exist('z_space', 'var')
    z_space = (sample_pixelsize/(5.4*um)) *  10*mm;  % psf��z����  z-length of PSF
end
if ~exist('z_pixelsize', 'var')
    z_pixelsize = (sample_pixelsize/(5.4*um)) *  200*um; % 12.5 * um;   % psf��z����  z step size of PSF 
end
z_list = (-z_space/2:z_pixelsize:z_space/2);  % psfѡȡ��z����  z coordinate in PSF
if strcmp(flag_pf_method, 'incoherent')==1  
    z_list = z_list - z_shift; 
    fprintf('[info] incoherent pf detected. `z_shift` added. \n');
elseif strcmp(flag_pf_method, 'gs')==1
    z_list = z_list - z_shift; 
    fprintf('[info] gs pf detected. `z_shift` added. \n');
end
z_cnt = numel(z_list);     % psf�Ĳ��� count of layers of PSF
lambda_0 = 405*nm;         % ���������� Laser wavelength
lambda = lambda_0 / n_ri ;   % ����ʱ�ĹⲨ��  wavelength during propagation
pixelsize = sample_pixelsize / pixelsize_ex_ratio; %5.4*um / pixelsize_ex_ratio ; %pixelsize = 2.5*um;  % ����ʱ�����سߴ� pixelsize during propagation

if ~exist('the_ratio', 'var')
    the_ratio = sample_pixelsize / 5.4e-6; 
else
    if abs(sample_pixelsize/5.4e-6/the_ratio -1)>0.01
        fprintf(2, '\n[Error]\n[Error] please check `sample_pixelsize = the_ratio * 5.4e-6`\n\n');
    end
end

if ~exist('source_margin', 'var') || ~exist('source_n', 'var')
    if strcmp(flag_pf_method, 'gs') == 1
        source_margin = 50;  % ����Ԥ�������أ�Ϊ���psfԤ�� pixels on edges, reserved for convolution PSF
    else
        source_margin = 0;
    end
    
    source_n = floor( sqrt(2*sample_n^2+sample_nz^2)*sample_pixelsize/pixelsize ) + 2*source_margin;  % ����ʱ��ʹ������ߴ� img size during propagation
    if mod(source_n, 2) ~= mod(sample_n, 2)
        source_n = source_n + 1;
    end
end

% Ĭ�ϵ�ͶӰγ�� default projection latitude
if ~exist('proj_latitude_0', 'var')
    proj_latitude_0 = 0.25 * pi;
end

proj_latitude = pi/2 - asin(  sin(pi/2 - proj_latitude_0) * 1 / n_ri );
stretch_ratio = cos(pi/2-proj_latitude_0) / cos(pi/2-proj_latitude);

%% ǿ��ͼ

if ~exist('flag_intensity_map', 'var')
    flag_intensity_map = 0; 
end

if ~exist('the_intensity_map', 'var')
    the_intensity_map = NaN; 
end
if ~exist('the_intensity_map_list', 'var')
    the_intensity_map_list = NaN; 
end
if ~exist('the_intensity_map_list_raw', 'var')
    the_intensity_map_list_raw = cell(1, proj_cnt); 
end

if flag_intensity_map && flag_rotate45
    fprintf(2, '\n[Error]\n[Error] flag_intensity_map may not work with flag_rotate45==1. \n\n');
end

if flag_intensity_map && exist('the_intensity_map_filename', 'var')
    if exist(the_intensity_map_filename, 'file')
        the_intensity_map = imread(the_intensity_map_filename);
        the_intensity_map = single(the_intensity_map);
        the_intensity_map = the_intensity_map / max(the_intensity_map(:));
        the_intensity_map = imresize(the_intensity_map, 1/the_ratio);

        % �Զ�����ͼ������
        proj_latitude = pi/2 - asin(  sin(pi/2 - proj_latitude_0) * 1 / n_ri );
        stretch_ratio = cos(pi/2-proj_latitude_0) / cos(pi/2-proj_latitude);
        the_intensity_map_list = cell(1, proj_cnt);
        the_intensity_map_list_raw = cell(1, proj_cnt);
        for i = 1:proj_cnt
            the_intensity_map_temp = the_intensity_map;

            % shift
            dmd_shift_radius = -10 / the_ratio;
            the_angle = (i-1)*2*pi/proj_cnt + pi/4;
            shift_ex = [0,0];
            shift_ex(1) = round( -sin(the_angle) * dmd_shift_radius );
            shift_ex(2) = round( -cos(the_angle) * dmd_shift_radius );
            the_intensity_map_temp = my_imtranslate(the_intensity_map_temp, -shift_ex);
            
            the_intensity_map_list_raw{i} = my_reshape_img(the_intensity_map_temp, [source_n, source_n]);
            the_intensity_map_list_raw{i}(the_intensity_map_list_raw{i}<0.1) = 0.1;

            % stretch
            the_angle = pi/2 + (i-1)*2*pi/proj_cnt ;
            if flag_rotate45 == 0    
                the_angle = the_angle + pi/4;  % DMD45����б����
            end

            the_intensity_map_temp = my_stretch_img(the_intensity_map_temp, 1/stretch_ratio, the_angle);
            the_intensity_map_temp = my_reshape_img(the_intensity_map_temp, [source_n, source_n]);
            the_intensity_map_list{i} = the_intensity_map_temp; 
            the_intensity_map_list{i}(the_intensity_map_list{i}<0.1) = 0.1;

        end
        %the_intensity_map = my_reshape_img(the_intensity_map, [source_n, source_n]);

    else
        fprintf(2,'\n[Error]\n[Error] file `%s` does not exist.\n\n', the_intensity_map_filename);
        return; 
    end
end

if ~exist('intensity_attenuation', 'var')  % ˥��ϵ�� Attenuation coefficient
    intensity_attenuation = 0.5;
end
attenuation_mask = intensity_attenuation .^ (((1:sample_nz)-(sample_nz+1)/2) * sample_pixelsize / sqrt(1 - (sin(pi/2-proj_latitude_0)/n_ri)^2) / cm);


%% Ԥ����PSF;  precompute PSF
if flag_GPU
    H_list = my_propagation_quick3D_init_GPU(source_n, source_n, z_list, pixelsize, lambda, 'Fresnel');
    my_propagation_fun = @my_propagation_quick3D_GPU;
else
    H_list = my_propagation_quick3D_init    (source_n, source_n, z_list, pixelsize, lambda, 'Fresnel');
    my_propagation_fun = @my_propagation_quick3D;
end
% ѡ��H_list���Ӽ�; select a subset of H_list
temp_z_n_pre = round( sqrt(2*sample_n^2+sample_nz^2) * sample_pixelsize / z_pixelsize); 
if temp_z_n_pre<1
    temp_z_n_pre = 1;
end
if temp_z_n_pre<=numel(H_list)
    [ll, rr] = my_get_reshape_range(temp_z_n_pre, numel(H_list) );
    z_H_list =  ll:rr ;
else
    fprintf(2, '[Warning] temp_z_n_pre>numel(H_list), please check `z_space`\n');
    z_H_list = 1:numel(H_list);
end

%% ����Ҷ���ǿ���������

if ~exist('flag_source_fft_error', 'var')
    flag_source_fft_error = 0;
end
if ~exist('the_source_fft_mask', 'var') || size(the_source_fft_mask,1)~=source_n
    the_source_fft_mask = ones(source_n,'single');
end

%% ��������
the_eps_sample = 1e-20; 
the_eps_source = 1e-20;
source_max = 100;

if exist('the_sample', 'var')
    the_sample(the_sample<the_eps_sample) = the_eps_sample; 
end

%% �����趨
switch flag_pb_method
    case 'gs'
        my_pb_fun = @(temp_mat)my_pb_gs_fres(temp_mat, H_list, source_n, my_propagation_fun);  % G-S iterative algorithm
    case 'gs0'
        my_pb_fun = @(temp_mat)my_pb_gs_fres(temp_mat, H_list, source_n, my_propagation_fun, 0);  % G-S iterative algorithm without energy correction.
    case 'mean'
        my_pb_fun = @(temp_mat)my_pb_mean(temp_mat); 
    case 'sum'
        my_pb_fun = @(temp_mat)my_pb_sum(temp_mat, 1);
    case 'max'
        my_pb_fun = @(temp_mat)my_pb_max(temp_mat); 
    case 'binary'
        my_pb_fun = @(temp_mat)my_pb_binary(temp_mat);
    case 'int'
        my_pb_fun = @(temp_mat)my_pb_intensity(temp_mat, intensity_attenuation, z_pixelsize, 1); 
    case 'int2'
        my_pb_fun = @(temp_mat)my_pb_intensity_2(temp_mat, intensity_attenuation, z_pixelsize, 1); 
    case 'intmax'
        my_pb_fun = @(temp_mat)my_pb_intensity_max(temp_mat, intensity_attenuation, z_pixelsize); 
    otherwise
        fprintf(2, '\n[Error]\n[Error] Unknown flag_pb_method `%s`, use `mean`\n', flag_pb_method);
        my_pb_fun = @(temp_mat)my_pb_mean(temp_mat); 
end

switch flag_pf_method
    case 'gs'
        my_pf_fun = @(the_source)my_pf_gs_fres(the_source, H_list, z_H_list, my_propagation_fun);  % fresnel propagation
    case 'easy'
        my_pf_fun = @(the_source)my_pf_easy(the_source, numel(z_H_list));   % linear propagation
    case 'incoherent'
        if ~exist('pixelsize_incoherent', 'var')
            %pixelsize_incoherent = 16.2e-6;
            pixelsize_incoherent = 4.5e-6; 
        end
        if flag_GPU
            H_list_incoherent = my_propagation_quick3D_init_GPU(source_n, source_n, z_list, pixelsize_incoherent, lambda, 'Fresnel');
        else
            H_list_incoherent = my_propagation_quick3D_init    (source_n, source_n, z_list, pixelsize_incoherent, lambda, 'Fresnel');
        end
        F1 = zeros(source_n, 'single');
        F_center = floor((source_n+1)/2);
        F1(F_center, F_center) = 1;
        F1_stack_a = my_propagation_fun(F1, H_list_incoherent, z_H_list); 
        F1_stack_a = abs(F1_stack_a) .^ 2;
        F1_stack_s = F1_stack_a;
        for ii = 1:size(F1_stack_a, 3)
            F1_stack_a1 = F1_stack_a(:, :, ii);
            F1_stack_s1 = imresize(F1_stack_a1, pixelsize_incoherent/pixelsize, 'bilinear'); 
            F1_stack_s1 = my_reshape_img(F1_stack_s1, [source_n, source_n] );
            
            F1_stack_a1_energy = sum(F1_stack_a1, 'all');
            F1_stack_s1_energy = sum(F1_stack_s1, 'all');
            F1_stack_s1 = F1_stack_s1 / F1_stack_s1_energy * F1_stack_a1_energy;
            
            F1_stack_s(:, :, ii) = F1_stack_s1;
        end
        my_pf_fun = @(the_source)my_pf_incoherent_2(the_source, F1_stack_s, z_H_list);
        
    case 'int'
        attenuation_mask = intensity_attenuation .^ (((1:sample_nz)-(sample_nz+1)/2) * sample_pixelsize / sqrt(1 - (sin(pi/2-proj_latitude_0)/n_ri)^2) / cm);
        my_pf_fun = @(the_source)my_pf_easy(the_source, numel(z_H_list));   % linear propagation
    otherwise
        fprintf(2, '\n[Error]\n[Error] Unknown flag_pf_method `%s`, use `easy`\n', flag_pf_method);
        my_pf_fun = @(the_source)my_pf_easy(the_source, numel(z_H_list));   
end

%% Jaccard���

Jaccard_weight = ones([sample_n, sample_n, sample_nz], 'single');  % ��ͨ Jaccard index

if ~exist('flag_Jaccard_weight', 'var')
    flag_Jaccard_weight = 0;
end
if ~exist('Jaccard_weight_r', 'var')
    Jaccard_weight_r = round(0.025 * sample_n);
end
if ~exist('Jaccard_weight_detail_lambda', 'var')
    Jaccard_weight_detail_lambda = 2;
end
if ~exist('thresh_cof', 'var')
    thresh_cof = 1;
end
if flag_Jaccard_weight 
    % ǿ���ض������ Jaccard index    
	if exist('Jaccard_weight_detail_filename', 'var')  
        % ֱ�Ӵ�ָ�����ļ���ȡ������Jaccard_weight
		Jaccard_weight_detail = loadtiff(Jaccard_weight_detail_filename);
		Jaccard_weight_detail = single(Jaccard_weight_detail);
		Jaccard_weight_detail = Jaccard_weight_detail / max(Jaccard_weight_detail(:));
        Jaccard_weight_detail = my_reshape_img(Jaccard_weight_detail, [sample_n, sample_n, sample_nz]);
		Jaccard_weight = Jaccard_weight + Jaccard_weight_detail_lambda * Jaccard_weight_detail;
    else
        % �Զ�����Jaccard_weight
		if sample_n>=500 || sample_nz>=500
			fprintf('[Warning][Jaccard] The sample size is too large. Calculating the desired matrix may take a significant amount of time.\n'); 
		end
		if Jaccard_weight_r>0
			the_sample_temp = imopen(the_sample, strel("sphere", Jaccard_weight_r));
			Jaccard_weight_detail = the_sample - the_sample_temp;
		else
			fprintf('[Warning] Please ensure that Jaccard_weight_r>0 if you intend to use Jaccard_weight.\n');
			Jaccard_weight_detail = Jaccard_weight;
		end
		Jaccard_weight = Jaccard_weight + Jaccard_weight_detail_lambda * Jaccard_weight_detail;
    end
    my_Jaccard_index_fun = @(img1, img2)my_Jaccard_index_weighted_polymer(img1, img2, Jaccard_weight);
elseif thresh_cof > 1
    my_Jaccard_index_fun = @(img1, img2)my_Jaccard_index_weighted_polymer(img1, img2, Jaccard_weight);
else
    my_Jaccard_index_fun = @(img1, img2)jaccard(logical(img1), logical(img2));
end

%% ������ postprocessing function
my_postprocess_source_list_fun = @(source_list)my_postprocess_source_list(source_list, source_n, the_eps_source, source_max, source_margin);
my_postprocess_source_list_fun_allownegative = @(source_list)my_postprocess_source_list(source_list, source_n, -source_max, source_max, source_margin); 


%% ����lossRecorder

lossRecorder = my_LossRecorder();
lossRecorder.addfun(@(img1, img2)my_MSE(img1, img2), 'MSE', '%d');
lossRecorder.addfun2(@(img1, img2)my_Jaccard_index_get_th(img1, img2), 'Jaccard', '%.3f(th=%.3f)');
if thresh_cof>1 || flag_Jaccard_weight
    lossRecorder.addfun2(@(img1, img2)my_Jaccard_index_get_th(img1, img2, thresh_cof, my_Jaccard_index_fun), 'Jaccard+w', '%.3f(th=%.3f)');
end

