function the_source = my_pb_max(temp_mat)
%MY_PB_SUM 选取投影最大值。这个函数不能用于迭代。
%  Get the projected maximum. This function cannot be used to calculate descent direction.

the_source = max(temp_mat, [], 3);

end
