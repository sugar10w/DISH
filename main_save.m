%% save source_list

% 追加阈值参数
str_notes_suffix = sprintf('%s_J%.2fth%.2f', str_notes_suffix, sel_Jaccard_index, sel_th_best);

%source_list_temp = my_3d_plus(source_list, my_3d_multip(source_list_residual, best_residual_ratio));
save( sprintf('%s/%s_coarse_%s.mat', the_dir_savename, sample_id, str_notes_suffix ),  'source_list', 'flag_rotate45');  

%% 拉伸后的结果 Save the stretched results 
if abs(proj_latitude_0)>pi/360
    source_list_stretch = my_stretch_source_list(source_list, proj_latitude_0, n_ri, flag_rotate45);
    % save( sprintf('%s/%s_coarse_strh_%s.mat', the_dir_savename, sample_id, str_notes_suffix ),  'source_list_stretch', 'flag_rotate45');  
end

%% 

weight_str = '';
if flag_Jaccard_weight
    weight_str = sprintf('_weight-r=%d', Jaccard_weight_r);
end

thresh_str = '';
if thresh_cof > 1
    thresh_str = sprintf('_thresh-cof=%.2f', thresh_cof);
end

%%
if exist('flag_watch_result', 'var') && flag_watch_result
    the_image = my_project_forward_fun(source_list);
    if numel(the_image)>2e8
        my_save_tiff_data_type = 8;
    else
        my_save_tiff_data_type = 16;
    end
    my_save_complex_tiff(the_image, sprintf('%s/%s_coarse_%s%s%s.tif', the_dir_savename, sample_id, str_notes_suffix, weight_str, thresh_str), '10', my_save_tiff_data_type);
end

