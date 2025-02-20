% Traditional interation without considering refraction and holography

addpath('util/source_list/');

%% 使用 Adam 优化器对 Jaccard index 进行非线性优化

if ~exist('flag_PM_Adam_run', 'var')
    flag_PM_Adam_run = true;

    alpha_ori = 0.3;  % alpha=初始搜索步长best_residual_ratio*alpha_ori
    beta1 = 0.3;
    beta2 = 0.7;
    epsilon = 1e-8;
    m = cellfun(@(x) single(zeros(size(x))), source_list, "UniformOutput", false);
    v = 0;

    lr_rate = 0.5;  % 学习率更新比率
    min_lr = alpha_ori * lr_rate^2;
    % thresh_ctrl = 1;
    thresh_ctrl = 1.5;
    
    tgt_th = 1.2;
    max_rate = 1.3;
    lr_cnt = 0;
    best_jac = 0;
    best_th = 0;
    best_sv_jac = 0;

    sam_bin = imbinarize(the_sample);

    if ~flag_reinforcement
        reinf_sample = the_sample;
    end

    str_notes_suffix_raw = str_notes_suffix;
else
    source_list = temp_source_list;
    the_image = temp_the_image;
    str_notes_suffix = str_notes_suffix_raw;
    source_list_sv = source_list;
end

if ~exist('jaccard_max_iter', 'var')
    jaccard_max_iter = 50;
end

fprintf('[info] Using Adam optimizer for PM reinforcement, α=%f, β1=%f, β2=%f, ε=%d\n', alpha_ori, beta1, beta2, epsilon);



%%

the_target_loss_fun = @(the_image)get_target_loss(the_image, reinf_sample, tgt_th, thresh_ctrl, sam_bin);
the_residual_fun = @(the_image)get_residual(the_image, reinf_sample, tgt_th, thresh_ctrl, sam_bin);


%%

Jaccard_list = nan(1, jaccard_max_iter+1);
th_best_list = nan(1, jaccard_max_iter+1);
loss_MSE_list = nan(1, jaccard_max_iter+1);
target_loss_list = nan(1, jaccard_max_iter+1);

% 原本的jaccard index

best_target_loss = the_target_loss_fun(the_image);

[Jaccard_index, img2_th_best] = my_Jaccard_index_get_th(the_sample, the_image, thresh_cof, my_Jaccard_index_fun);
if ~exist('sel_Jaccard_index', 'var')
    sel_Jaccard_index = Jaccard_index;
    sel_th_best = img2_th_best;
end

raw_Jaccard_index = Jaccard_index;
raw_img2_th_best = img2_th_best;

Jaccard_list(1) = Jaccard_index;
th_best_list(1) = img2_th_best;
loss_MSE_list(1) = my_MSE(the_sample * max(1, img2_th_best), the_image);
target_loss_list(1) = best_target_loss;

last_target_loss = best_target_loss;
cnt_ok = 0;
cnt_ok_trigger_num = 4;

%%
if exist('i__', 'var')
    iter_st = i__+1;
else
    iter_st = 1;
end

time_start = datetime;

for i__ = iter_st:(jaccard_max_iter+iter_st-1)

    fprintf('\n=-=-=-=-=-=-=-=-=-=-PMAr Jaccard: %02d-=-=-=-=-=-=-=-=-=-=\n', i__);

    the_residual = the_residual_fun(the_image);

    %%
    source_list_residual = my_project_backward_fun_allownegative(the_residual);
    source_list_residual = my_postprocess_source_list_fun_allownegative(source_list_residual);  % 下降方向
    
    m = cellfun(@(mm, g) beta1 * mm + (1 - beta1) * g, m, source_list_residual, "UniformOutput", false);
    v = beta2 * v + (1 - beta2) * sum(cell2mat(cellfun(@(g) sum(g.^2, "all"), source_list_residual, "UniformOutput", false)), "all");
    m_hat = cellfun(@(mm) mm / (1 - beta1^i__), m, "UniformOutput", false);
    v_hat = v / (1 - beta2^i__);
    real_dir = cellfun(@(mh) mh ./ ((v_hat.^.5)+epsilon), m_hat, "UniformOutput", false);

    %% 假定随着ratio的增加，target_loss只会出现一个峰值；在三个等间距点的区间上更新；
    if i__ == 1
        the_image_residual = my_project_forward_fun(real_dir);
        best_target_loss_0 = -inf;
        best_residual_ratio = 0;
        best_img2_th_best = raw_img2_th_best;
    
        the_ratio_region = 0.05;
        
        l_ratio = 0;
        
        l_Ja = -the_target_loss_fun( the_image_residual * l_ratio + the_image );
        r_ratio = 1;
        
        r_Ja = -the_target_loss_fun( the_image_residual * r_ratio + the_image );
        l_th=tgt_th; r_th=tgt_th;
        for ext_i = 1:20
            m_Ja = r_Ja;
            m_th = r_th;
            m_ratio = r_ratio;
    
            r_ratio = m_ratio * 2;
            
            r_Ja = -the_target_loss_fun(the_image_residual * r_ratio + the_image);

            if r_Ja < m_Ja
                break
            end
        end
        fprintf('[Info] Final l_ratio: %f, l_Ja: %f, m_ratio: %f, m_Ja: %f, r_ratio: %f, r_Ja: %f.\n', l_ratio, l_Ja, m_ratio, m_Ja, r_ratio, r_Ja);
        
        ratio_list_temp = [
            l_ratio, (l_ratio + m_ratio) / 2, m_ratio, (m_ratio + r_ratio) / 2, r_ratio;
            l_Ja,                          0,    m_Ja,                       0,    r_Ja;
            l_th,                          0,    m_th,                       0,    r_th;
            ];
        cnt_temp = 0;
        lm_th=tgt_th; mr_th=tgt_th;
        while ratio_list_temp(1, 5) - ratio_list_temp(1, 1) > the_ratio_region
            cnt_temp = cnt_temp + 1;
    
            ratio_list_temp(1, 2) = mean(ratio_list_temp(1, [1, 3]));
            ratio_list_temp(1, 4) = mean(ratio_list_temp(1, [3, 5]));
            
            lm_Ja = -the_target_loss_fun(the_image_residual * ratio_list_temp(1, 2) + the_image);
            mr_Ja = -the_target_loss_fun(the_image_residual * ratio_list_temp(1, 4) + the_image);
            ratio_list_temp(2, 2) = lm_Ja;
            ratio_list_temp(3, 2) = lm_th;
            ratio_list_temp(2, 4) = mr_Ja;
            ratio_list_temp(3, 4) = mr_th;
    
            [max_Ja, max_ind] = max(ratio_list_temp(2, :));
    
            if best_target_loss_0 < max_Ja
                best_target_loss_0 = max_Ja;
                best_residual_ratio = ratio_list_temp(1, max_ind);
                best_img2_th_best = ratio_list_temp(3, max_ind);
            end
    
            m_ind = min(4, max(2, max_ind));
            ratio_list_temp(:, 1) = ratio_list_temp(:, m_ind-1);
            ratio_list_temp(:, 5) = ratio_list_temp(:, m_ind+1);
            ratio_list_temp(:, 3) = ratio_list_temp(:, m_ind);
        end

        alpha = alpha_ori * max(best_residual_ratio, 300);
    end
    
    source_list = cellfun(@(x, rd) min(max(0, x + alpha*rd), proj_intensity), source_list, real_dir, "UniformOutput", false);
    
    the_image = my_project_forward_fun(source_list);
    [Jaccard_index, img2_th_best] = my_Jaccard_index_get_th(the_sample, the_image, thresh_cof, my_Jaccard_index_fun);
    
    fprintf('[Info][PMAr][%d] before:   Jac=%d, th=%d\n', i__, raw_Jaccard_index, raw_img2_th_best);
    % fprintf('[Info][PMAr][%d] expected: Jac=%d, th=%d\n', i__, best_Jaccard_index, best_img2_th_best);
    fprintf('[Info][PMAr][%d] after:    Jac=%d, th=%d\n', i__, Jaccard_index, img2_th_best);
    fprintf('[Info][PMAr][%d] ratio=%d\n', i__, alpha);
    
    % flag_watch_iter=1;
    if flag_watch_iter
        save(sprintf('%s/mid_[%d]_PMAr_all.mat', the_dir_savename, i__));  
        
        temp_str_notes_suffix = sprintf('%s_J%.2fth%.2f', str_notes_suffix, Jaccard_index, img2_th_best);
        save( sprintf('%s/mid_[%d]_PMAr_source_%s.mat', the_dir_savename, i__, temp_str_notes_suffix ),  'source_list', 'flag_rotate45');  
        my_save_complex_tiff(the_image, sprintf('%s/mid_[%d]_PMAr_output_%s.tif', the_dir_savename, i__, temp_str_notes_suffix), '10');
        if abs(proj_latitude_0)>pi/360
            source_list_stretch = my_stretch_source_list(source_list, proj_latitude_0, n_ri, flag_rotate45);
            save( sprintf('%s/mid_[%d]_PMAr_source_strh_%s.mat', the_dir_savename, i__, temp_str_notes_suffix ),  'source_list_stretch', 'flag_rotate45');  
        end
    end
    


    %%

    curr_target_loss = the_target_loss_fun(the_image);
    
    fprintf('[Info][PMAr] current target loss = %d\n', curr_target_loss);

    if curr_target_loss>last_target_loss
        alpha = alpha * lr_rate;
        cnt_ok = 0;
        cnt_ok_trigger_num = cnt_ok_trigger_num + 1;
    else
        
        if curr_target_loss<best_target_loss
            best_target_loss = curr_target_loss;
            source_list_sv = source_list;  % update source_list_sv
    
            cnt_ok = cnt_ok + 1;
            if cnt_ok>=cnt_ok_trigger_num
                alpha = alpha*1.5;
                cnt_ok = 0;
            end
        else
            cnt_ok = 0;
        end

    end

    last_target_loss = curr_target_loss;

    %%

    if (Jaccard_index + .01 * img2_th_best)>best_sv_jac
        %source_list_sv = source_list;
        sel_Jaccard_index = Jaccard_index;
        sel_th_best = img2_th_best;
        best_sv_jac = Jaccard_index + .01 * img2_th_best;
    end
    
    if Jaccard_index > best_jac
        best_jac = Jaccard_index;
    end
    if img2_th_best > best_th
        best_th = img2_th_best;
    end

    raw_Jaccard_index = Jaccard_index;
    raw_img2_th_best = img2_th_best;
    Jaccard_list(i__-iter_st+2) = Jaccard_index;
    th_best_list(i__-iter_st+2) = img2_th_best;
    loss_MSE_list(i__-iter_st+2) = my_MSE(the_sample * max(1, img2_th_best), the_image);
    target_loss_list(i__-iter_st+2) = curr_target_loss;

    %%
    fprintf('[info][PMAr] round %d(%d) finished, %s\n', i__, jaccard_max_iter+iter_st-1, my_predict_time(time_start, (i__-iter_st+1) / jaccard_max_iter));
end

temp_the_image = the_image;
temp_source_list = source_list;
source_list = source_list_sv;

%% 展示迭代效果

% figure(NumberTitle="off", Name="MSE");
% plot((iter_st-1):i__, loss_MSE_list(1:(i__-iter_st+2)), DisplayName='PMAr', LineWidth=1); hold on; legend;
% title('MSE');

figure(1);
plot((iter_st-1):i__, Jaccard_list(1:(i__-iter_st+2)), DisplayName='PMAr', LineWidth=1); hold on; %legend;
title('jaccard');

figure(2);
plot((iter_st-1):i__, th_best_list(1:(i__-iter_st+2)), DisplayName='PMAr', LineWidth=1); hold on; %legend;
title('image threshold');

figure(3);
plot((iter_st-1):i__, target_loss_list(1:(i__-iter_st+2)), DisplayName='PMAr', LineWidth=1); hold on; %legend;
title('target loss');

%%

function the_target_loss = get_target_loss(the_image, reinf_sample, tgt_th, thresh_ctrl, sam_bin)
img_low_bin = imbinarize(the_image, tgt_th / thresh_ctrl);
the_target_loss = mean( sam_bin .* (reinf_sample * thresh_ctrl * tgt_th - the_image).^2 ...
            + ~sam_bin .* img_low_bin .* (the_image - tgt_th / thresh_ctrl).^2, 'all');
end

function the_residual = get_residual(the_image, reinf_sample, tgt_th, thresh_ctrl, sam_bin)
img_low_bin = imbinarize(the_image, tgt_th / thresh_ctrl);
the_residual = sam_bin .*                 (reinf_sample * thresh_ctrl * tgt_th - the_image) ...
             - ~sam_bin .* img_low_bin .* (the_image - tgt_th / thresh_ctrl);
end

