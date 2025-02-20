function [Jaccard_index] = my_Jaccard_index_weighted(img1, img2, w)
%my_Jaccard_index  处理两个二值图像，计算Jaccard index
% Calculate Jaccard index of two binary images

img1_bin = logical(img1);
img2_bin = logical(img2);

Jaccard_index = sum((img1_bin & img2_bin) .* w, 'all') / sum((img1_bin | img2_bin) .* w, 'all');

% Jaccard_index = jaccard(img1, img2);

end

