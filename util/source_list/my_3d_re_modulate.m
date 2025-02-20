function [source_list_new] = my_3d_re_modulate(source_list, mode1, mode2)
%MY_3D_RE_MODULATE  调整source_list的总强度
%  Adjust the energy of each source_list

if nargin<3
    mode2 = 'none';
end


proj_cnt = numel(source_list);


% get the_intensity
the_intensity = my_get_re_modulate_intensity(proj_cnt, mode1, mode2);

% get source_list_new
source_list_new = cell(1, proj_cnt);
for i = 1:proj_cnt
    img = source_list{i} * the_intensity(i);
    source_list_new{i} = single( rand(size(img)) < img); 
end

end

