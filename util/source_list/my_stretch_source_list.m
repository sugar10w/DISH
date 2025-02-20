function source_list_stretch = my_stretch_source_list(source_list, proj_latitude_0, n_ri, flag_rotate45, ...
    varargin)

flag_reverse = 0;

for i = 1:numel(varargin)
    switch varargin{i}
        case 'reverse'
            flag_reverse = 1;
        otherwise
            continue;
    end
end

%%

proj_cnt = numel(source_list);


%%
proj_latitude = pi/2 - asin(  sin(pi/2 - proj_latitude_0) * 1 / n_ri );
stretch_ratio = cos(pi/2-proj_latitude_0) / cos(pi/2-proj_latitude);

if flag_reverse  % source_list_stretch -> source_list
    stretch_ratio = 1/stretch_ratio;
end


source_list_stretch = cell(1, proj_cnt);
for i = 1:proj_cnt
    the_angle = pi/2 + (i-1)*2*pi/proj_cnt ;
    if flag_rotate45 == 0    
        the_angle = the_angle + pi/4;  % DMD45∂»«„–±∑≈÷√
    end
    source_list_stretch{i} = my_stretch_img(source_list{i}, stretch_ratio, the_angle );
end



end