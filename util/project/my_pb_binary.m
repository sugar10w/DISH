function the_source = my_pb_binary(temp_mat)
%my_pb_binary ��ö�ֵ����ͶӰ�������������������ڼ����½�����
%  Get the binarized projection result. This function cannot be used to calculate descent direction.

the_source = sum(temp_mat, 3) > 0.02;  %�����ֵ���ֶ����õģ�����ֵ��������� Threshold can be adjusted as needed

%the_source = sum(temp_mat, 3) > 0.274;
%the_source = sum(temp_mat, 3) > 0.4;

end

