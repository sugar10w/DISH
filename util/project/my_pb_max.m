function the_source = my_pb_max(temp_mat)
%MY_PB_SUM ѡȡͶӰ���ֵ����������������ڵ�����
%  Get the projected maximum. This function cannot be used to calculate descent direction.

the_source = max(temp_mat, [], 3);

end
