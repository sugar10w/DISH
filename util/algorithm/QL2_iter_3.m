function X = QL2_iter_3(X, Y, my_project_forward_fun, my_project_backward_fun, my_postprocess_source_list_fun, ...
    iter_max, sample_id, proj_intensity, flag_watch_iter, lossRecorder)
% QL2_ITER    min|Y-AX|, 0<=X<=max

if nargin<10
    lossRecorder = my_LossRecorder();
    lossRecorder.addfun(my_MSE, 'MSE', '%d');
end

Y = Y / max(Y(:));  % 最大值设置为1，方便分析 The maximum value is set to 1
loss_MSE_list = zeros(1, iter_max);
% loss_MSE_auto_list = zeros(1, iter_max);
% loss_Jac_auto_list = zeros(1, iter_max);
% loss_Jac_auto_th_list = zeros(1, iter_max);

time_start = datetime;

X_std = my_project_backward_fun(Y);  % 先获得理想的投影图 ideal projection
X_std = my_postprocess_source_list_fun(X_std);

l_search_min = 5e-4;
l_search_max = Inf;

for iter = 1:iter_max
    fprintf('=-=-=-=-=-=-%d=-=-=-=-=-=-=-=\n', iter);
    Y_curr = my_project_forward_fun(X);
    
    % 记录目前情况 current situation
    fprintf('target max = %d, curr max = %d \n', max(Y(:)), max(Y_curr(:)) );
    % [~, img2_th_best] = my_Jaccard_index_get_th(Y, Y_curr);
    loss_MSE_list(iter) = my_MSE(Y, Y_curr);
    % loss_MSE_auto_list(iter) = my_MSE_auto(Y, Y_curr);
    % [loss_Jac_auto_list(iter), loss_Jac_auto_th_list(iter)] = my_Jaccard_index_get_th(Y, Y_curr);
    % fprintf('MSE loss = %d  ( auto: %d ) \n', loss_MSE_list(iter), loss_MSE_auto_list(iter)  );
    % fprintf('Jaccard_index_auto = %.3f  \n',  loss_Jac_auto_list(iter)  );
    lossRecorder.update(Y, Y_curr);

    % 用类似梯度下降的方式迭代  Iterate in a gradient descent-like fashion
    X_curr = my_project_backward_fun(Y_curr);
    X_curr = my_postprocess_source_list_fun(X_curr);
    
    % 梯度方向就是(X_curr{i} - X_std{i})   The gradient direction is (X_curr{i} - X_std{i})
    A = loss_MSE_list(iter);
    
    grad_dir = cellfun(@(x,y) x - y, X_std, X_curr, "UniformOutput", false);
    
    for l_search = [1, 2]  % 搜索步长 search step
        X_temp = cellfun(@(x, gd) min(max(0, x + l_search*gd), proj_intensity), X, grad_dir, "UniformOutput", false);

        % 给出loss
        Y_temp = my_project_forward_fun(X_temp);
        loss_l2_temp = my_MSE(Y, Y_temp);
        fprintf('l_search=%d, MSE_temp=%d\n', l_search, loss_l2_temp);

        % 
        if l_search==1
            B = loss_l2_temp;
        elseif l_search==2
            C = loss_l2_temp;
        end

    end
    
    % 忽略约束条件会出问题，注意调整搜索步长或梯度方向的系数 Ignoring constraints may cause problems, try adjusting the search step size or the coefficients in the gradient direction
    l_search_best = - ( - 1.5*A + 2*B  - 0.5*C ) / ( A - 2*B + C );    
    if isnan(l_search_best)
        l_search_best = 0;
    end
    if l_search_best<l_search_min
        fprintf('[Info][QL2_iter] cannot update, stop here.\n');
        break;
    end

    l_search_best = min(l_search_max, l_search_best);
    X = cellfun(@(x, gd) min(max(0, x + l_search_best*gd), proj_intensity), X, grad_dir, "UniformOutput", false);
    %     Y_temp = my_project_forward_fun(X);
    %     loss_l2_temp = sqrt ( mean( abs( Y - Y_temp ).^2 , [1,2,3]) );
    %     fprintf('l_search=%d, loss_l2_temp=%d\n', l_search_best, loss_l2_temp);
    fprintf('l_search=%d\n', l_search_best);
    
    % 保存中间过程 save intermediate process
    if flag_watch_iter
        % my_save_complex_tiff(Y_curr, sprintf('./result_samples/%s/%03d.tif', sample_id, iter), '10');
        save(['./result_samples/', sample_id, '/', num2str(iter), '_source.mat'], 'X');  
    end
    
    fprintf('[QL2 iter] %s\n', my_predict_time(time_start, iter/iter_max) );
    
end

% fprintf('MSE\tMSE auto\tJac auto\n');
% for iter = 1:iter_max
%     fprintf('%d\t%d\t%d(%d)\n', loss_MSE_list(iter), loss_MSE_auto_list(iter), loss_Jac_auto_list(iter), loss_Jac_auto_th_list(iter));
% end
lossRecorder.review();

% figure(NumberTitle="off", Name="MSE");
% plot(loss_MSE_list(1:iter), DisplayName='L2', LineWidth=1); hold on; legend;
% title('MSE');

end

