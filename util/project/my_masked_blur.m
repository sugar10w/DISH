function img_blur = my_masked_blur(img, mask, kernel)
%MY_MASKED_BLUR 带mask的blur

mask = single(mask);

if max(mask(:)) ~= 1
    fprintf(2, '[Error][my_masked_blur] max(mask(:)) ~= 1\n');
end

img_1 = img .* mask;
img_1 = imfilter(img_1, kernel);

img_2 = imfilter(mask, kernel);

img_blur = img_1 ./ img_2;
img_blur(mask==0) = 0;

end

