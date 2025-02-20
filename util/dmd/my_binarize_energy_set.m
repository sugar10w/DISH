function [img_bw_best, img_th_best] = my_binarize_energy_set(img, energy)
% 进行二值化 binarize
%  要求二值化后的图像能量尽可能符合所设定的能量energy
%  The image energy after binarization is required to match `energy`


% 如果`energy<1`，那改成像素个数 If `energy<1`, change it to the number of pixels
if energy<1
    energy = energy * numel(img);
end
energy = ceil(energy);
energy = max(0, energy);
energy = min(energy, numel(img) );

% 初始化 initialization
img_th_l = 0;
img_th_r = max( img(:) );

% 固定次数的二分法 dichotomy
for iter = 1:20  
    
    img_th_new = (img_th_l + img_th_r)/2;
    energy_new = sum( img>img_th_new, 'all' );
    
    if energy_new==energy
        img_th_l = img_th_new;
        img_th_r = img_th_new;
        break;
    elseif energy_new<energy
        % 说明img_th_new太大了 img_th_new is too large
        img_th_r = img_th_new;
    elseif energy_new>energy
        % 说明img_th_new太小了 img_th_new is too small
        img_th_l = img_th_new;
    end
    
end


% 决定最终结果 decide the final result
energy_l = sum( img>img_th_l, 'all' );
energy_r = sum( img>img_th_r, 'all' );

if abs(energy-energy_l)<abs(energy-energy_r)
    img_th_best = img_th_l;
else
    img_th_best = img_th_r;
end
img_bw_best = img > img_th_best;


end

