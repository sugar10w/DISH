function source_list_new = my_binarize_source_list_multistep(source_list, proj_ratio, max_intensity)
% 这个函数并不真的产生二值化图像，而是产生与之等效的阶梯化的图案

if nargin<3
    max_intensity = my_get_source_list_max(source_list);
end
if nargin<2
    proj_ratio = 3;
end

proj_cnt = numel(source_list);

source_list_new = cell(1, proj_cnt);
for i = 1:proj_cnt
    img = source_list{i} / max_intensity;
    
    img_0 = zeros(size(img), 'single');
    for ii = 1:proj_ratio
        img_0 = img_0 + ( img > ((ii-0.5)/proj_ratio) );
    end
    img_0 = img_0 / proj_ratio;
    
    source_list_new{i} = img_0;
end


end
