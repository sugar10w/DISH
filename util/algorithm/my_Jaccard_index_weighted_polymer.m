function [Jaccard_index] = my_Jaccard_index_weighted_polymer(img1, img2, w)
%my_Jaccard_index  处理一个二值图像和一个灰度图像（0~1），计算Jaccard index
% Calculate Jaccard index of a binary image and a grayscale image

img1_bin = logical(img1);
% img2_bin = logical(img2);
img2 = min(img2, 1);
img2 = max(img2, 0);

Jaccard_index = sum(img1_bin .* img2 .* w, 'all') / sum(max(img1_bin, img2) .* w, 'all');

% Jaccard_index = jaccard(img1, img2);

end