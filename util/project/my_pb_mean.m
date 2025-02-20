function the_source = my_pb_mean(temp_mat)
%MY_PB_MEAN ���ͶӰƽ��ֵ
%  get the projected mean

temp_mat_no_zero = sum(temp_mat~=0, 3);
the_source = sum(temp_mat, 3) ./ temp_mat_no_zero;  % ƽ������;  mean

the_source(isinf(the_source)) = 0;
the_source(isnan(the_source)) = 0;

end

