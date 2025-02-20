function info_str = my_predict_time(time_start, ratio)
% ͨ����ʼʱ��͵�ǰ���ȣ��������ʱ���
%  Calculate the end time point based on the start time and the current progress

time_elapsed = datetime - time_start;
time_predict = time_elapsed / ratio;
time_end = time_start + time_predict;

info_str = sprintf('duration= %s; estimated end time= %s;', char(time_elapsed), char(time_end));

end

