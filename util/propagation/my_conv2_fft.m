function [img_conv2] = my_conv2_fft(img_a, img_b)
%my_conv2_fft 使用FFT，计算卷积 
% Use FFT to compute convolution

[ra, ca]=size(img_a(:,:,1));
[rb, cb]=size(img_b(:,:,1));

r = ra+rb-1;
c = ca+cb-1;

img_aa = zeros(r,c, 'single');
img_bb = zeros(r,c, 'single');

img_aa(1:ra, 1:ca) = img_a;
img_bb(1:rb, 1:cb) = img_b;
img_conv2_fft = fft2(img_aa) .* fft2(img_bb);

img_conv2 = ifft2(img_conv2_fft);

end

