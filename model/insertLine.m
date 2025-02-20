function out = insertLine(model_or_size, start_pnt, end_pnt, width, type)
%INSERTLINE 此处显示有关此函数的摘要
%   此处显示详细说明
if ndims(model_or_size) == 3
    out = model_or_size;
    model_size = size(model_or_size);
else
    out = false(model_or_size);
    model_size = model_or_size;
end
real_start_pnt = start_pnt;
real_end_pnt = end_pnt;
ext = ceil(width);
lef_but = floor(min(real_start_pnt, real_end_pnt));
start_pnt = real_start_pnt - lef_but + ext + 1;
end_pnt = real_end_pnt - lef_but + ext + 1;

temp_model = false(2*ext + ceil(max(real_start_pnt - lef_but, real_end_pnt - lef_but)) + 1);
model_size = size(temp_model);
[y_mat, x_mat, z_mat] = meshgrid(1:model_size(2), 1:model_size(1), 1:model_size(3));

vect1 = [x_mat(:),y_mat(:),z_mat(:)] - start_pnt;
vect2 = [x_mat(:),y_mat(:),z_mat(:)] - end_pnt;
cro = sqrt(sum(cross(vect1, repmat((end_pnt - start_pnt) / norm(end_pnt - start_pnt),...
    size(x_mat(:), 1), 1), 2) .^ 2, 2)) <= width;
seg1 = dot(vect1, repmat((end_pnt - start_pnt), size(x_mat(:), 1), 1), 2) >= 0;
seg2 = dot(vect2, repmat((start_pnt - end_pnt), size(x_mat(:), 1), 1), 2) >= 0;
temp_model(cro & seg1 & seg2) = true;

if nargin > 4 && strcmp(type, 'ball')
    start_ball = sum(vect1.^2, 2) <= width^2;
    end_ball = sum(vect2.^2, 2) <= width^2;
    temp_model(start_ball | end_ball) = true;
end

m_st = max([1, 1, 1], lef_but - ext - 1);
m_ed = min(size(out), lef_but - ext + model_size - 2);
t_st = ext + 2 - lef_but + m_st;
t_ed = ext + 2 - lef_but + m_ed;

out(m_st(1):m_ed(1), m_st(2):m_ed(2), m_st(3):m_ed(3)) = ...
    out(m_st(1):m_ed(1), m_st(2):m_ed(2), m_st(3):m_ed(3)) | ...
    temp_model(t_st(1):t_ed(1), t_st(2):t_ed(2), t_st(3):t_ed(3));

end

