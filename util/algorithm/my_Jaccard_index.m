function [Jaccard_index] = my_Jaccard_index(img1, img2)
%my_Jaccard_index  处理两个二值图像，计算Jaccard index
% Calculate Jaccard index of two binary images

img1 = logical(img1);
img2 = logical(img2);

Jaccard_index = jaccard(img1, img2);

end

