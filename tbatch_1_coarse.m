% Coarse optimization without considering refraction and holography

%% clean workspace
if exist('the_dir_savename', 'var') && exist('str_notes_suffix', 'var')
    source_list_init_filename = sprintf('%s/PMAr_source_%s.mat', the_dir_savename, str_notes_suffix );
    if exist(source_list_init_filename, "file")
        keep_vars = {'source_list_init_filename'};
        expr = sprintf('clearvars -except %s', strjoin(keep_vars, ' '));
        evalin('base', expr);
    else
        clear;
    end
end

%%
addpath('util/shape/');
addpath('util/tiff/');

time_start = datetime; 

%% input filename

sample_filename = './model/PointsTest/points p-3.tif';  % input filename 
source_margin=0; source_n=60;    % pattern size
the_ratio = 1;                   % sample_pixelsize = 5.4e-6 * the_ratio;
jaccard_max_iter = 20;           % max iteration of the coarse optimization
proj_intensity = 4.0;            % relative intensity of the light source * exposure time

% sample_filename = './model/Misc/Bird300.tif';
% source_margin=0; source_n=300; 
% the_ratio = 1; 
% jaccard_max_iter = 20;
% proj_intensity = 8.0;


%%
the_sample = loadtiff(sample_filename);
if mod(size(the_sample,1)+source_n, 2) ~= 0
    source_n = source_n + 1;
end

%%

n_ri = 1.45;    % refractive index

intensity_attenuation =  0.6;  % remaining energy after passing through 1cm of material

%%
flag_iter_method = 'none';
flag_init_method = 'auto';
flag_pb_method = 'sum'; 
flag_pf_method = 'int'; 
flag_process = 1;

if ~exist('the_ratio', 'var')
    the_ratio = 1;
end
sample_pixelsize = 5.4e-6 * the_ratio;  %  default pixelsize is 5.4um


flag_rotate45 = 0; %  0->Considering the 45 degree placement of the DMD in the projection matrix

flag_watch_iter = 0;    
flag_watch_result = 1;  

pb_num_workers = 1;
pf_num_workers = 1;

flag_intensity_map = 1;
the_intensity_map_filename = './data/beam_profile_247.png';  % beam profile


%% initialization
main2_easy;

if exist('the_image', 'var')
    [Jaccard_index, img2_th_best] = my_Jaccard_index_get_th(the_sample, the_image, thresh_cof, my_Jaccard_index_fun);
    if img2_th_best>2
        for i = 1:numel(source_list)
            source_list{i} = source_list{i} / img2_th_best;
        end
        the_image = the_image / img2_th_best;
    end
end

%% main
if flag_process && jaccard_max_iter>0
    main_Adam_reinf;
    main_save;
end


%%
time_lapse = datetime - time_start; 
fprintf('[Info] Total run time is %s\n', char(time_lapse));

