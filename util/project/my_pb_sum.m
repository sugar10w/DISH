function the_source = my_pb_sum(temp_mat, n)
%MY_PB_SUM ѡȡͶӰ�ͣ�Ȼ����Թ̶�ϵ��
%  Get the projected sum, then divide it by a fixed factor

if nargin<2
    n = 50;
end

the_source = sum(temp_mat, 3) / n;

end
