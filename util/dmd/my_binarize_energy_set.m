function [img_bw_best, img_th_best] = my_binarize_energy_set(img, energy)
% ���ж�ֵ�� binarize
%  Ҫ���ֵ�����ͼ�����������ܷ������趨������energy
%  The image energy after binarization is required to match `energy`


% ���`energy<1`���Ǹĳ����ظ��� If `energy<1`, change it to the number of pixels
if energy<1
    energy = energy * numel(img);
end
energy = ceil(energy);
energy = max(0, energy);
energy = min(energy, numel(img) );

% ��ʼ�� initialization
img_th_l = 0;
img_th_r = max( img(:) );

% �̶������Ķ��ַ� dichotomy
for iter = 1:20  
    
    img_th_new = (img_th_l + img_th_r)/2;
    energy_new = sum( img>img_th_new, 'all' );
    
    if energy_new==energy
        img_th_l = img_th_new;
        img_th_r = img_th_new;
        break;
    elseif energy_new<energy
        % ˵��img_th_new̫���� img_th_new is too large
        img_th_r = img_th_new;
    elseif energy_new>energy
        % ˵��img_th_new̫С�� img_th_new is too small
        img_th_l = img_th_new;
    end
    
end


% �������ս�� decide the final result
energy_l = sum( img>img_th_l, 'all' );
energy_r = sum( img>img_th_r, 'all' );

if abs(energy-energy_l)<abs(energy-energy_r)
    img_th_best = img_th_l;
else
    img_th_best = img_th_r;
end
img_bw_best = img > img_th_best;


end

