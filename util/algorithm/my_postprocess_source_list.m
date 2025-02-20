function source_list = my_postprocess_source_list(source_list, source_n, the_eps_source, source_max, source_margin)
%MY_POSTPROCESS_SOURCE_LIST 后处理
% post-process

for i = 1:numel(source_list)
    % 重新设置尺寸 reshape
    source_list{i} = my_reshape_img( source_list{i}, [source_n, source_n]);
    % 规避无效数据类型 avoid invalid data types
    if islogical(source_list{i})
        source_list{i} = single(source_list{i});
    end
    % 强制转化为正实数 coerce to positive real number
    % source_list{i} = abs(source_list{i});
    source_list{i} ( source_list{i}<the_eps_source ) = the_eps_source; % 最小值限制 min
    source_list{i} ( source_list{i}>source_max )     = source_max;     % 最大值限制 max
    % 边上清空掉，给psf留位置 clear the edge and leave room for PSF
    source_list{i} = my_clean_margin(source_list{i}, source_margin, the_eps_source);
end

end

