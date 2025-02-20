function source_list = my_postprocess_source_list(source_list, source_n, the_eps_source, source_max, source_margin)
%MY_POSTPROCESS_SOURCE_LIST ����
% post-process

for i = 1:numel(source_list)
    % �������óߴ� reshape
    source_list{i} = my_reshape_img( source_list{i}, [source_n, source_n]);
    % �����Ч�������� avoid invalid data types
    if islogical(source_list{i})
        source_list{i} = single(source_list{i});
    end
    % ǿ��ת��Ϊ��ʵ�� coerce to positive real number
    % source_list{i} = abs(source_list{i});
    source_list{i} ( source_list{i}<the_eps_source ) = the_eps_source; % ��Сֵ���� min
    source_list{i} ( source_list{i}>source_max )     = source_max;     % ���ֵ���� max
    % ������յ�����psf��λ�� clear the edge and leave room for PSF
    source_list{i} = my_clean_margin(source_list{i}, source_margin, the_eps_source);
end

end

