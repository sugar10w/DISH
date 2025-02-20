function the_source = my_pb_intensity_2(temp_mat, intensity_attenuation, z_pixelsize, n)
%MY_PB_INTENSITY 先进行强度衰减矫正，然后得到投影和，最后除以固定系数
%   Perform intensity attenuation correction first,
% then obtain the projection sum, 
% and finally divide by a fixed coefficient `n`

if nargin<4
    n = 50;
end

temp_z_n = size(temp_mat, 3);
temp_z_center = (temp_z_n + 1)/2;
for i = 1:temp_z_n
    the_attenuation = intensity_attenuation ^ ( (-i+temp_z_center) * z_pixelsize / 1e-2 ) ;
    temp_mat(:, :, i) = temp_mat(:, :, i) * the_attenuation;
end

the_source = sum(temp_mat, 3) / n;

end

