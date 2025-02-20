function X = Richardson_Lucy_iter(X, Y, my_project_forward_fun, my_project_backward_fun, my_postprocess_source_list_fun, ...
    iter_max, sample_id, proj_intensity, flag_watch_iter)
% Richardson_Lucy 迭代算法
%  Richardson_Lucy iterative algorithm

Y = Y / max(Y(:));  % 最大值设置为1，方便分析 The maximum value is set to 1
the_eps_sample = 1e-3; 
Y(Y<the_eps_sample) = the_eps_sample;  % 为了算法稳定性考虑，暗值不能太低 For stability, the value of dark part cannot be too low

loss_MSE_list = zeros(1, iter_max);
loss_MSE_auto_list = zeros(1, iter_max);
loss_Jac_auto_list = zeros(1, iter_max);
loss_Jac_auto_th_list = zeros(1, iter_max);

time_start = datetime;

for iter = 1:iter_max
    fprintf('=-=-=-=-=-=-%d=-=-=-=-=-=-=-=\n', iter);
    Y_curr = my_project_forward_fun(X);
    
    % 输出目前情况 current situation
    fprintf('target max = %d, curr max = %d \n', max(Y(:)), max(Y_curr(:)) );
    loss_MSE_list(iter) = my_MSE(Y, Y_curr);
    loss_MSE_auto_list(iter) = my_MSE_auto(Y, Y_curr);
    [loss_Jac_auto_list(iter), loss_Jac_auto_th_list(iter)] = my_Jaccard_index_get_th(Y, Y_curr);
    fprintf('MSE loss = %d  ( auto: %d ) \n', loss_MSE_list(iter), loss_MSE_auto_list(iter)  );
    fprintf('Jaccard_index_auto = %.3f  \n',  loss_Jac_auto_list(iter)  );
    
    % Richardson–Lucy 
    the_error = Y ./ Y_curr;
    the_error(the_error>10000) = 10000;            % 最大值抑制 Max
    the_error = my_clean_margin(the_error, 5, 1);  % 清理边缘 clean up the edges
    the_error(the_error>1 & Y<=2*the_eps_sample) = 1;  % 暗部处理 adjust the value of dark part
    
    the_error_X = my_project_backward_fun(the_error);
    the_error_X = my_postprocess_source_list_fun(the_error_X);

    % 更新source_list; update source_list
    for i = 1:numel(X)
        X{i} = X{i} .* the_error_X{i};
        X{i}(X{i} > proj_intensity) = proj_intensity;
    end
    X = my_postprocess_source_list_fun(X);
    
    % 特殊操作:非线性变换 nonlinear transformations
    % X_max = my_get_source_list_max(X); 
    % if proj_intensity>0 && X_max>proj_intensity   % 如果强度太低的话，不要进行非线性变换；防止整张图都被压制为0；If the intensity is too low, do not perform nonlinear transformations; or the entire image may be suppressed to 0;
    %     fprintf('[RL iter] Nonlinear process...\n');
    %     for i = 1:numel(X)
    %         X{i} = my_nonlinear_01( X{i}, proj_intensity);
    %     end
    % end
    
    % 保存中间过程 save intermediate process
    if flag_watch_iter
        my_save_complex_tiff(Y_curr, sprintf('./result_samples/%s/%03d.tif', sample_id, iter), '10');
        save(['./result_samples/', sample_id, '/', num2str(iter), '_source.mat'], 'X');  
        
        the_error(the_error>1000)=1000; 
        my_save_complex_tiff(the_error, sprintf('./result_samples/%s/error_%03d.tif', sample_id, iter), '10');
        save(['./result_samples/', sample_id, '/', num2str(iter) ,'_error_source.mat'], 'the_error_X');  
    end
    
    fprintf('[RL iter] %s\n', my_predict_time(time_start, iter/iter_max) );
    
end


fprintf('MSE\tMSE auto\tJac auto\n');
for iter = 1:iter_max
    fprintf('%d\t%d\t%d(%d)\n', loss_MSE_list(iter), loss_MSE_auto_list(iter), loss_Jac_auto_list(iter), loss_Jac_auto_th_list(iter));
end


end