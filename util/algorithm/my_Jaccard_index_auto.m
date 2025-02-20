function [Jaccard_index] = my_Jaccard_index_auto(img1, img2, my_Jaccard_index_fun)
%my_Jaccard_index  穷举阈值，获得最优的 Jaccard index ； 要求img1已经是二值的；
% Enumerate thresholds to obtain the optimal Jaccard index; img1 should be binary;

if nargin<3
    my_Jaccard_index_fun = @(img1, img2)my_Jaccard_index(img1, img2);
end

% if max(img1(:)) ~= 1
%     fprintf('[Info] unacceptable img1 (max(img1(:)) ~= 1) . `my_Jaccard_index_auto` does not work   \n');
%     Jaccard_index = my_Jaccard_index_fun(img1, img2);
%     return;
% end
% 
% img1 = img1 > 0.5;
% 
% img2_max = max(img2(:));
% 
% Jaccard_index_best = 0;
% img2_th_best = 0;
% 
% for th = 0:0.05:1
%     img2_th = th * img2_max;
%     temp = my_Jaccard_index_fun(img1, img2>img2_th);
%     
%     if Jaccard_index_best<temp
%         Jaccard_index_best = temp;
%         img2_th_best = img2_th;
%     end
%     
% end
% 
% Jaccard_index = Jaccard_index_best;
% fprintf('[Info] my_Jaccard_index_auto img2 > %.2f \n', img2_th_best);

[Jaccard_index, ~] = my_Jaccard_index_get_th(img1, img2, my_Jaccard_index_fun);

end

