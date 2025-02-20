function the_source = my_pb_intensity_max(temp_mat, intensity_attenuation, z_pixelsize)
%MY_PB_INTENSITY �Ƚ���ǿ��˥��������Ȼ��õ�ͶӰ�ͣ������Թ̶�ϵ��
%   Perform intensity attenuation correction first,
% then obtain the projection sum, 
% and finally divide by a fixed coefficient `n`

temp_z_n = size(temp_mat, 3);
temp_z_center = floor( (temp_z_n + 1)/2 );
for i = 1:temp_z_n
    the_attenuation = intensity_attenuation ^ ( (-i+temp_z_center) * z_pixelsize / 1e-2 ) ;
    temp_mat(:, :, i) = temp_mat(:, :, i) * the_attenuation;
end


the_source = max(temp_mat, [], 3);

end

