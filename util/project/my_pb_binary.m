function the_source = my_pb_binary(temp_mat)
%my_pb_binary 获得二值化的投影结果。这个函数不能用于计算下降方向。
%  Get the binarized projection result. This function cannot be used to calculate descent direction.

the_source = sum(temp_mat, 3) > 0.02;  %这个阈值是手动设置的，其他值或许更合适 Threshold can be adjusted as needed

%the_source = sum(temp_mat, 3) > 0.274;
%the_source = sum(temp_mat, 3) > 0.4;

end

