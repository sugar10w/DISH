function [MSE_auto] = my_MSE_auto(img1, img2)
%my_MSE_a 计算 MSE ，允许img2乘上系数以获得更小的MSE
% Compute MSE; allowing img2 to be multiplied by a factor to get a smaller MSE

if max(img1(:)) ~= 1
    fprintf('[Warning] max(img1)~=1. my_MSE_auto may not work, \n');
end

img2_t = sum(img1.*img2, 'all') / sum( img2.^2, 'all' );
%fprintf('[Info] my_MSE_auto: img2 / %d \n', 1/img2_t);

if ~isnan(img2_t)
    MSE_auto = my_MSE( img1, img2*img2_t );
else
    fprintf(2, '[Warning][MSE_auto] img2==0\n');
    MSE_auto = my_MSE( img1, img2 );
end

end

