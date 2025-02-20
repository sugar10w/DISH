function the_image_0 = my_pf_gs_fres_s2(the_source_air, ...
    H_list, z_H_list, H_list_air, H_list_material, ...
    my_propagation_fun, proj_i_curr_ratio, stretch_ratio, flag_rotate45, ...
    varargin)
%MY_PF_GS_FRES 向前传播，fresnel传播; Fresnel propagation
%  the_source 和 the_image 都要求是能量
%  考虑了折射界面的位置
%  `the_source` and `the_image` are required to be energy/intensity fields
%  The position of the refraction interface is considered

flag_source_fft_error = 0;
the_source_fft_mask = ones(size(H_list{1}),'single');

% pixelsize = 5.4e-6;
% theta_r = 28.5823 / 180 * pi;
% lambda = 405e-9;

my_stretch_img_refraction_fun = @(the_img, the_angle, flag_from_air_to_material) ...
    my_stretch_img_easy(the_img, the_angle, flag_from_air_to_material, stretch_ratio);

for i = 1:numel(varargin)-1
    if ~isstring(varargin{i}) && ~ischar(varargin{i})
        continue;
    end
    switch varargin{i}
        case 'SourceFftError'
            flag_source_fft_error = 1;
            temp = varargin{i+1};
            the_source_fft_mask = temp;
        case 'StretchFun'
            temp = varargin{i+1};
            if isa(temp, 'function_handle')
                my_stretch_img_refraction_fun = temp;
            end
    end
end

%%

if nargin<9
    flag_rotate45 = 1;  % 默认不进行修正
end

%%
the_angle = pi/2 + proj_i_curr_ratio*2*pi ;   % 拉伸角度
if flag_rotate45 == 0
    the_angle = the_angle + pi/4;  % DMD45度倾斜放置
end

the_source_air = the_source_air.^0.5;       % 空气中的DMD（假想的共轭面上）
if flag_source_fft_error
    the_source_air_e = ifft2( fft2(the_source_air) .* the_source_fft_mask );
else
    the_source_air_e = the_source_air; 
end

the_interface_air = my_propagation_fun(the_source_air_e, H_list_air);   %the_interface_air 空气侧角谱

flag_from_air_to_material = 1;
the_interface_material = my_stretch_img_refraction_fun(the_interface_air, the_angle, flag_from_air_to_material ); %the_interface_material 物料侧角谱
the_interface_material = my_reshape_img(the_interface_material, size(the_source_air));

the_source_material = my_propagation_fun(the_interface_material, H_list_material); %the_source_material 物料内假想的共轭面处的场

the_image_0 = my_propagation_fun(the_source_material, H_list, z_H_list);
the_image_0 = abs(the_image_0) .^ 2;

end

