function the_image_0 = my_pf_gs_fres(the_source, H_list, z_H_list, my_propagation_fun)
%MY_PF_GS_FRES ��ǰ������fresnel����; Fresnel propagation
%  the_source �� the_image ��Ҫ��������
%  `the_source` and `the_image` are required to be energy/intensity fields

the_image_0 = my_propagation_fun(the_source.^0.5, H_list, z_H_list);
the_image_0 = abs(the_image_0) .^ 2;

end

