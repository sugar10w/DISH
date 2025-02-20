% generate DMD images

addpath('./util/dmd/');
addpath('./util/shape');
addpath('./util/source_list')

%%
the_dir_savename_dmd = [ './result_dmd/' , 'display images dir' ];    % 输出文件夹名 output directory name
source_resize_ratio = 1;    % 图像缩放倍率  image zoom ratio

i_base = 0;

the_radius = -10;   % 黑色件，倾斜投影
dmd_shift_mode =  'auto';         % 补偿电机误差的偏移，定时更新 Compensate the offset of motor error; update it regularly
flag_multi_dmd_mode = 0;  % 是否同时使用多个DMD模式
dmd_shift_mode_list = {'c2A', 'c2B'};

flag_dmd_shift_extra = 0;      % 垂直投影专用，用于移动中心点位置
dmd_shift_extra = [+65, +74];


flag_intensity_re = 0;         % 是否进行强度矫正 whether to perform intensity correction
flag_intensity_re_filename = './data/beam_profile_247.png';   %强度矫正所用的图片名 Image name used for intensity correction


flag_intensity_attenuation = 0;  % 是否进行简易的衰减强度矫正
if flag_intensity_attenuation
    max_intensity_for_attenuation = 1.3;   % 限制衰减强度矫正后的最大强度
    theta_i = 45 / 180*pi;
    theta_r = 30 / 180*pi;
    intensity_attenuation = 0.6;
    intensity_attenuation_pixel = intensity_attenuation ^ ( 5.4e-4/cos(theta_i)/sin(theta_r) );
end


flag_1bit_mode = 1;  % 是否以1bit模式保存 Whether to save in 1bit mode
if flag_1bit_mode
    %flag_binary_method = 'threshold';   %    threshold 阈值化方法
    %img_threshold = 0.5;  % 二值化阈值 binarization threshold

    flag_binary_method = 'random';     % 随机方法 random 
else
    flag_binary_method = 'grayscale';
end

%%

if exist('source_list_new', 'var')
    if numel(source_list_new) == 180
        source_list_new = my_binarize_multiply_source_list_motionblur(source_list_new, 10);
        fprintf('my_binarize_multiply_source_list\n');
    end
elseif exist('source_list_stretch','var')
    if numel(source_list_stretch) == 180
        source_list_stretch = my_binarize_multiply_source_list_motionblur(source_list_stretch, 10);
        fprintf('my_binarize_multiply_source_list\n');
    end
end

%%


if ~exist('dmd_shift_mode', 'var')
    fprintf('[Info] get_dmd_shift_ex: use [none] shift \n');
    dmd_shift_mode = 'none';
end

if ~exist('flag_dmd_shift_extra', 'var')
    flag_dmd_shift_extra = 0;
    dmd_shift_extra = [0, 0];
end

if ~exist('flag_multi_dmd_mode', 'var')
    flag_multi_dmd_mode = 0;
    dmd_shift_mode_list = {dmd_shift_mode};
end

if ~exist('move_3d_vector', 'var')
    move_3d_vector = [0,0,0];
end

% 将source_list中的内容转化为dmd可读取的形式

if exist('source_list_output', 'var')
    source_list = source_list_output;
    fprintf('[Info] use source_list_output\n');
elseif exist('source_list_new', 'var') % 优先使用 source_list_new
    source_list = source_list_new; 
    fprintf('[Info] use source_list_new\n');
elseif exist('source_list_stretch', 'var') % 优先使用 source_list_stretch
    source_list = source_list_stretch; 
    fprintf('[Info] use source_list_stretch\n');
end



if ~exist('flag_rotate45', 'var')
    flag_rotate45 = 1;   % 是否额外旋转45度 Whether to rotate 45 degrees
end


if ~exist('flag_intensity_re', 'var')
    flag_intensity_re = 0;   % 是否通过DMD进行强度修正  whether to perform intensity correction
    if ~exist('flag_intensity_re_filename', 'var')
        flag_intensity_re_filename = 'img_intensity_re_0.2.bmp';
    end
end

if ~exist('flag_1bit_mode', 'var')
    flag_1bit_mode = 1;       % 是否保存为1bit图像  Whether to save in 1bit mode
end


% 如果没有source_list，则生成默认图像
flag_auto_genrate_point_test = 0;
if ~exist('source_list', 'var')  % 生成默认图像
    


    %rot_motor_img_shift_deg = 0;
    
    flag_auto_genrate_point_test = 1;
    
    img = zeros(1080, 1080);
    
    img(539:541, 539:541) = 1; % 3pix
    
    source_list = cell(1,1800);
    for i = 1:numel(source_list)
        source_list{i} = img;
    end 
    
    source_resize_ratio = 1;
    
    the_dir_savename_dmd = [ './result_dmd/' , 'point test' ];    % 输出文件夹名 output directory name
    if norm(move_3d_vector)>0
        the_dir_savename_dmd = sprintf('%s move%.0f', the_dir_savename_dmd, norm(move_3d_vector)); 
    end
    
    
    flag_intensity_attenuation = 0;

    fprintf(2,'[Warning] 没有发现source_list, 自动生成中心点序列\n');
end




%%
proj_cnt = numel(source_list);

for i = 1:proj_cnt
    if isscalar(source_list{i}) && isnan(source_list{i})
        source_list{i} = 0;
    elseif isempty(source_list{i})
        source_list{i} = 0;
    end
end

flag_test = 0;   % 是否只输出第一张 whether to output only the first image

if ~exist('source_resize_ratio', 'var')
    source_resize_ratio = 1;  % 图像缩放倍率  image zoom ratio
end

if ~exist('rot_motor_img_shift_deg', 'var')
    rot_motor_img_shift_deg =  6.55; %  0; %   1.3; %        % 1.3deg 对应是T=4s转速， 6.5deg 对应是T=0.8s转速
end
rot_motor_img_shift = round(rot_motor_img_shift_deg / (360/proj_cnt));  % 为了配合电机的位置，DMD需要显示前n个时刻的图像；In order to match the position of the motor, the DMD needs to display the images at n moments in advance;



if ~exist('flag_binary_method', 'var')
    %flag_binary_method = 'threshold';   %    threshold 阈值化方法
    flag_binary_method = 'random';     % random 随机方法
end




%%
if flag_intensity_re
    img_intensity_re = imread(flag_intensity_re_filename); 
    img_intensity_re = double(img_intensity_re);
    img_intensity_re = img_intensity_re / max( img_intensity_re(:) );
end

%%
if ~exist('i_base', 'var')
    i_base = 0;           % 输出起始序号 output starting sequence number
end
sample_id = '';           % 输出图像前缀名 output image prefix

if ~exist('the_dir_savename_dmd', 'var')   % 自动生成文件夹名 generate directory name if necessary
    if exist('z', 'var')
        the_dir_savename_dmd = sprintf('./result_dmd/temp z%d/', z);
    else
        the_dir_savename_dmd = './result_dmd/temp/';
    end
end


the_dir_savename_dmd = sprintf('%s x%.1f proj%d rr%d r%d %s', ...
    the_dir_savename_dmd, source_resize_ratio, proj_cnt, rot_motor_img_shift, the_radius, flag_binary_method);
if flag_multi_dmd_mode
    the_dir_savename_dmd = sprintf('%s %s+%s', ...
        the_dir_savename_dmd, dmd_shift_mode_list{1}, dmd_shift_mode_list{2});
else
    the_dir_savename_dmd = sprintf('%s %s', ...
        the_dir_savename_dmd, dmd_shift_mode);
end
if flag_intensity_attenuation
    the_dir_savename_dmd = sprintf('%s att%.2fmax%.2f', ...
        the_dir_savename_dmd, intensity_attenuation, max_intensity_for_attenuation);
end
if flag_intensity_re
    [~, flag_intensity_re_filename_sub, ~] = fileparts(flag_intensity_re_filename);
    the_dir_savename_dmd = sprintf('%s re[%s]', ...
        the_dir_savename_dmd, flag_intensity_re_filename_sub); 
end
if flag_dmd_shift_extra
    the_dir_savename_dmd = sprintf('%s [Z]', ...
        the_dir_savename_dmd); 
end

%the_dir_savename_dmd = [the_dir_savename_dmd, ' ', dmd_shift_mode, ' x', num2str(source_resize_ratio), ' proj', num2str(proj_cnt), ' rr', num2str(rot_motor_img_shift)];

if the_dir_savename_dmd(numel(the_dir_savename_dmd)) ~= '/'  
    the_dir_savename_dmd = [the_dir_savename_dmd, '/'];
end

if ~exist(the_dir_savename_dmd, 'dir')  % 生成文件夹 make directory
    mkdir(the_dir_savename_dmd);
end
filename_output_format = [the_dir_savename_dmd, sample_id, '%04d.png']; % 输出路径 output filename


%%

size_dmd   = [1080, 1920];
size_dmd_c = [1080, 1080];

% 重新规划大小
if source_resize_ratio~=1 && source_resize_ratio~=0
    fprintf('source_resize_ratio=%d \n', source_resize_ratio);
    source_list_resize = cell(1, proj_cnt);
    
    img_temp_suggested_size = round(( max(size_dmd_c)+100 )/source_resize_ratio); % 这个大小的图像就足够用了
    for i = 1:proj_cnt
        
        img_temp = source_list{i};
        if max(size(img_temp)) > img_temp_suggested_size
            img_temp = my_reshape_img(img_temp, [img_temp_suggested_size, img_temp_suggested_size]);
        end
        
        source_list_resize{i} = imresize(img_temp, source_resize_ratio, 'bilinear');
    end
else
    source_list_resize = source_list;
end

% 调整亮度等级 brightness
if ~exist('max_intensity', 'var')
    max_intensity = 0;
end
if max_intensity==0
    for i = 1:proj_cnt
        if max( source_list_resize{i}(:))>0
            max_intensity = max(max_intensity, max( source_list_resize{i}(:)) );
        end
    end
end
for i = 1:proj_cnt
    source_list_resize{i} = source_list_resize{i} / max_intensity;
end
fprintf('max intensity = %d \n', max_intensity);

if norm(move_3d_vector)>1
    fprintf(2, 'moving source list by [%d %d %d]\n', move_3d_vector);
    source_list_resize = my_3d_move_source_list(source_list_resize, move_3d_vector, flag_rotate45);
end



%%

i_max = proj_cnt;
if flag_test   % 是否只输出第一张 whether to output only the first image
    i_max = 1;
end

parfor i = 1:i_max
    
    if mod(i,20)==0
        fprintf('%d / %d \n', i, proj_cnt);
    end
    
    img = source_list_resize{i};
    img = single(img);
    img(isnan(img)) = 0;

    img = my_reshape_img(img, size_dmd_c);   % 调整尺寸 reshape
    if flag_rotate45
        img = imrotate(img, -45  );            % DMD是倾斜放置的，大部分时候都要转 The DMD is placed obliquely and needs to be turned 45deg most of the time
    end
    img = my_reshape_img(img, size_dmd);
    
    
    % 进行简易衰减修正
    if flag_intensity_attenuation
        dmd_center_point = size_dmd / 2;
        
        the_theta = (i-1)/proj_cnt * 2*pi;
        the_angle_vector = [-cos(the_theta+pi/4), sin(the_theta+pi/4) ];
        
        the_attenuation_mask = ones(size_dmd);
        for x = 1:size_dmd(1)
            for y = 1:size_dmd(2)
                the_attenuation_mask(x, y) = intensity_attenuation_pixel ^ (([x,y]-dmd_center_point)*the_angle_vector');
            end
        end
        
        img = img .* the_attenuation_mask / max_intensity_for_attenuation;
        img(img>1) = 1;
        
    end
    
    
    % 适配DMD位置误差，手动额外修正！
    if flag_multi_dmd_mode
        img_accum = zeros(size(img));
        for dmd_shift_mode_i = 1:numel(dmd_shift_mode_list)
            dmd_shift_mode_temp = dmd_shift_mode_list{dmd_shift_mode_i};
            img_temp = my_imtranslate(img, get_dmd_shift_ex_auto( (i-1)/proj_cnt , dmd_shift_mode_temp, the_radius)); 
            img_accum = img_accum + img_temp; 
        end
        img = img_accum / numel(dmd_shift_mode_list); 
        if flag_1bit_mode && strcmp(flag_binary_method, 'threshold')==0
            fprintf('[Warning] `flag_multi_dmd_mode` works with thresholding or grayscales'); 
        end
    else
        img = my_imtranslate(img, get_dmd_shift_ex_auto( (i-1)/proj_cnt , dmd_shift_mode, the_radius)); 
    end
    
    % 垂直投影可能需要额外位移
    if flag_dmd_shift_extra
        img = my_imtranslate(img, dmd_shift_extra);
    end
    
    
    % 投影到DMD上的光强可能不平整
    if flag_intensity_re
        img = img .* img_intensity_re;
    end
    
    if flag_1bit_mode
        
        if strcmp(flag_binary_method, 'threshold')
            img = img >= img_threshold;
        elseif strcmp(flag_binary_method, 'random')
            img = ( sqrt(0.58^2+1.68*img)-0.58 ) / 0.84; % 随机化强度可能有非线性?? The randomization strength may be nonlinear ??
            img = img > rand(size(img));
        else
            img = img >= img_threshold;
        end
        
        img = 1 - img;        % 颠倒颜色 invert colors
        img = logical(img);
    else
        img = 1 - img;    % 颠倒颜色 invert colors
    end
    
    
    ii_rot_motor_shift = mod(i+rot_motor_img_shift-1, proj_cnt) +1;  % 电机滞后矫正; Motor lag correction
    
    imwrite(img, sprintf(filename_output_format, ii_rot_motor_shift+i_base  ));
    
end


fprintf('rr %d\n', rot_motor_img_shift );
fprintf('done \n');


if flag_auto_genrate_point_test
    clear; 
end
