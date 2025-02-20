function binaryMatrix = my_select_layered_details(temp_mat_sample)
% 三维矩阵，选出细节区域

the_temp_mat = temp_mat_sample;

if size(the_temp_mat, 3)==1
    binaryMatrix = true(size(temp_mat_sample));
    return;
end


%% 用高通滤波去提取高频成分

for z = 1:size(the_temp_mat, 3)
    the_layer = the_temp_mat(:, :, z);
    the_layer_blurred = imgaussfilt(the_layer, 1);
    the_layer_detail = abs(the_layer - the_layer_blurred);
    the_temp_mat(:, :, z) = max(the_layer_detail, imgaussfilt(the_layer_detail, 10));
end

the_temp_mat(temp_mat_sample<0.003) = 0;

z_center = ceil(1 + size(the_temp_mat, 3)/2);
the_temp_mat(:, :, z_center) = the_temp_mat(:, :, z_center) + 0.01;

%% 找出高频成分所在的层

[~, maxIndices] = max(the_temp_mat, [], 3);
binaryMatrix = false(size(the_temp_mat));

[I, J] = ndgrid(1:size(the_temp_mat,1), 1:size(the_temp_mat,2));
linearIndices = sub2ind(size(the_temp_mat), I, J, maxIndices);
binaryMatrix(linearIndices) = true;


%% 删除完全没有光强的区域

the_temp_mat_proj = max(temp_mat_sample, [], 3);
the_temp_mat_proj = imdilate(the_temp_mat_proj>1e-10, strel('disk',20));

binaryMatrix_center = binaryMatrix(:, :, z_center);
binaryMatrix_center(the_temp_mat_proj==0) = 0;
binaryMatrix(:, :, z_center) = binaryMatrix_center;

end