function the_image_0 = my_pf_intensity(the_source, size_z, intensity_list)
%my_pf_intensity�� ֱ�ߴ����������ض�ǿ�ȵ��� 
%  linear propagation with specific intensity modulation

the_image_0 = repmat(the_source, [1,1,size_z]); 

for z = 1:size_z
    the_image_0(:,:,z) = the_image_0(:,:,z) * intensity_list(z);
end

end

