function [threshold_list, Jaccard_index_list] = my_Jaccard_index_histogram(img1, img2, threshold_list, my_Jaccard_index_fun)
%my_Jaccard_index_histogram  穷举阈值，返回 Jaccard index 序列 ； 要求img1已经是二值的；
% Enumerate thresholds to obtain the optimal Jaccard index; img1 should be binary;

if nargin<4
    my_Jaccard_index_fun = @(img1, img2)my_Jaccard_index(img1, img2);
end

if nargin<3
    threshold_list = 0:0.05:1;
end
img2_th_best = zeros(size(threshold_list));

if max(img1(:)) ~= 1
    fprintf('[Info] unacceptable img1 (max(img1(:)) ~= 1) . `my_Jaccard_index_auto` does not work   \n');
    return;
end

img1 = img1 > 0.5;

img2_max = max(img2(:));


for i = 1:numel(threshold_list)
    th = threshold_list(i);
    img2_th = th * img2_max;
    temp = my_Jaccard_index_fun(img1, img2>img2_th);
    Jaccard_index_list(i) = temp;

end


end