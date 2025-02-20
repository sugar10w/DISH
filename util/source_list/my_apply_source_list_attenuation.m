function source_list_new = my_apply_source_list_attenuation(source_list, sample_pixelsize, intensity_attenuation, n_ri, flag_mode, flag_rotate45)

if nargin<6
    flag_rotate45 = 1; 
end

if nargin<5
    flag_mode = 'surface';
end


proj_cnt = numel(source_list);
source_list_new = cell(1, proj_cnt);


theta_i = pi/4;
theta_r = asin( sin(theta_i)/n_ri );

intensity_attenuation_pixel = intensity_attenuation ^ (sample_pixelsize/1e-2);

if strcmp(flag_mode, 'surface')  % 只矫正表面附近的衰减（从空气中的等相位面到溶液中的等相位面）
    intensity_attenuation_pixel = intensity_attenuation_pixel ^ ( 1/cos(theta_i)*sin(theta_r) );
elseif strcmp(flag_mode, 'surface_zaxial') % 从空气中的等相位面，到z轴柱状物，仅用于近似处理
    intensity_attenuation_pixel = intensity_attenuation_pixel ^ ( 1/cos(theta_i)/sin(theta_r) );
elseif strcmp(flag_mode, 'zaxial')  % 从溶液中的等相位面，到z轴柱状物，仅用于近似处理
    intensity_attenuation_pixel = intensity_attenuation_pixel ^ ( 1/tan(theta_r) );
else
    intensity_attenuation_pixel = intensity_attenuation_pixel ^ ( 1/cos(theta_i)*sin(theta_r) );
end


for i = 1:proj_cnt

    img = source_list{i};
    
    the_theta = (i-1)/proj_cnt * 2*pi;
    if flag_rotate45==0
        the_theta = the_theta + pi/4;
    end
    the_angle_vector = [-cos(the_theta), sin(the_theta)];

    the_attenuation_mask = ones(size(img));
    the_center_point = size(img)/2;
    for x = 1:size(the_attenuation_mask, 1)
        for y = 1:size(the_attenuation_mask, 2)
            the_attenuation_mask(x,y) = intensity_attenuation_pixel ^ (([x,y]-the_center_point)*the_angle_vector');
        end
    end
    
    img = img .* the_attenuation_mask;
    
    source_list_new{i} = img;

end


end