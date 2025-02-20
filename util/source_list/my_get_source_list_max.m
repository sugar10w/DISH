function [source_list_max] = my_get_source_list_max(source_list)
% 获得source_list的最大值  get maximum of source_list

source_list_max = 0;
for i = 1:numel(source_list)
    curr_max = max( source_list{i}, [], 'all' );
    if source_list_max < curr_max
        source_list_max = curr_max;
    end
end

end

