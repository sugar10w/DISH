function [source_list_new] = my_binarize_multiply_source_list_2( source_list, proj_ratio, max_intensity )
% ��source_list�е�ÿһ�ŻҶ�ͼ��ɶ����ֵͼ(�ͻҶ�ͼ)
%   Split each grayscale image into multiple binary images (and grayscale images)

if nargin<3
    max_intensity = my_get_source_list_max(source_list);
end

if nargin<2
    proj_ratio = 3;   % ��proj_cnt����proj_ratio����magnification
end

proj_cnt = numel(source_list);


%% ���¸ı��С
for i = 1:proj_cnt
    source_list{i} = source_list{i} / max_intensity;
end

%% ����С���У���������������չ�������� Generates small sequences, i.e. sequences spread out from the center to two wings

relative_idx = zeros(1, proj_ratio);
for i = 1:proj_ratio
    relative_idx(i) = (  -mod(i,2)*2+1  ) * ( floor(i/2) );
end

%%
source_list_new = cell(1, proj_cnt*proj_ratio);
source_list_sub = cell(1, proj_ratio);

for proj_i = 1:proj_cnt
    for gray_i = 1:proj_ratio
        source_list_sub{gray_i} = source_list{proj_i} > ((gray_i-0.5) / proj_ratio ) ;
    end

    % ��һ�ŻҶ�ͼ
    source_list_sub{proj_ratio} = max(source_list{proj_i}*proj_ratio - (proj_ratio-1), 0);
    
    proj_i_new_c = (proj_i-1)*proj_ratio + 1;
    for gray_i = 1:proj_ratio
        proj_i_new = mod( proj_i_new_c + relative_idx(gray_i) - 1, proj_cnt*proj_ratio ) + 1;
        source_list_new{ proj_i_new } = source_list_sub{gray_i};
    end
end


end