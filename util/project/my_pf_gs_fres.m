function the_image_0 = my_pf_gs_fres(the_source, H_list, z_H_list, my_propagation_fun)
%MY_PF_GS_FRES 向前传播，fresnel传播; Fresnel propagation
%  the_source 和 the_image 都要求是能量
%  `the_source` and `the_image` are required to be energy/intensity fields

the_image_0 = my_propagation_fun(the_source.^0.5, H_list, z_H_list);
the_image_0 = abs(the_image_0) .^ 2;

end

