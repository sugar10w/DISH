% Fine optimization with holography

clc; 

addpath('util\source_list\');

time_start = datetime; 

%% input filename

if ~exist('sample_filename', 'var')
    sample_filename = './model/PointsTest/points p-3.tif';
end
if ~exist('source_list', 'var')
    load('./result_samples/points p-3/points p-3_coarse_prj180_n1.45_none_auto_sum_int_max4.0_i0.60_p5.4um_I_J1.00th1.51.mat')
end
proj_intensity = my_get_source_list_max(source_list);
source_margin=40; source_n=130; the_ratio = 1;


% if ~exist('sample_filename', 'var')
%     sample_filename = './model/Misc/Bird300.tif';
% end
% if ~exist('source_list', 'var')
%     load('./result_samples/Bird300/Bird300_coarse_prj180_n1.45_none_auto_sum_int_max8.0_i0.60_p5.4um_I_J0.97th1.30.mat')
% end
% proj_intensity = my_get_source_list_max(source_list);
% source_margin=40; source_n=400; the_ratio = 1;


%%
binarize_ratio = 10;


%%

n_ri = 1.45; 

intensity_attenuation = 0.6; 


proj_cnt = numel(source_list);

sample_pixelsize = 5.4e-6 * the_ratio;

%%

interface_distance = 5.0e-3;   % 打印中心到前表面的距离 
z_shift = +0.0e-3;     % 共轭面到打印中心的传播距离 


%%
flag_rotate45 = 0; % 1 -> Ignoring the 45 degree placement of the DMD in the projection matrix

flag_watch_iter = 0; 
flag_watch_result = 1;

%
pixelsize_ex_ratio = 1;  
if pixelsize_ex_ratio~=1
    for proj_i = 1:proj_cnt
        source_list{proj_i} = my_resize_img(source_list{proj_i}, pixelsize_ex_ratio);
    end
end

pb_num_workers = 1;
pf_num_workers = 1;

flag_source_fft_error = 0;  

% the_intensity_map
flag_intensity_map = 1;

the_intensity_map_filename = './data/beam_profile_247.png'; 

%% main
main2_n_gs_s;


%%

time_lapse = datetime - time_start; 
fprintf('[Info] Total run time is %s\n', char(time_lapse));
