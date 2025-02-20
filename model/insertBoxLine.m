function out = insertBoxLine(model_or_size, st_pnt, ed_pnt, lw, v_pl, depth)
%INSERTLINE 此处显示有关此函数的摘要
%   此处显示详细说明
if ndims(model_or_size) == 3
    out = model_or_size;
    model_size = size(model_or_size);
else
    out = false(model_or_size);
    model_size = model_or_size;
end
real_start_pnt = st_pnt;
real_end_pnt = ed_pnt;
ext = ceil(max(lw, depth));
lef_but = floor(min(real_start_pnt, real_end_pnt));
st_pnt = real_start_pnt - lef_but + ext + 1;
ed_pnt = real_end_pnt - lef_but + ext + 1;

temp_model = false(2*ext + ceil(max(real_start_pnt - lef_but, real_end_pnt - lef_but)) + 1);
model_size = size(temp_model);
[y_mat, x_mat, z_mat] = meshgrid(1:model_size(2), 1:model_size(1), 1:model_size(3));

vect = [x_mat(:),y_mat(:),z_mat(:)] - (st_pnt + ed_pnt) / 2;

v_pnt = (ed_pnt - st_pnt);

v_n = cross(v_pnt, v_pl);
v_pl = cross(v_pnt, v_n);
dist1 = abs(dot(vect, repmat(v_n, size(x_mat(:), 1), 1), 2) / norm(v_n));
dist2 = abs(dot(vect, repmat(v_pl, size(x_mat(:), 1), 1), 2) / norm(v_pl));
dist3 = abs(dot(vect, repmat(v_pnt, size(x_mat(:), 1), 1), 2) / norm(v_pnt));
seg_dist1 = dist1 < (depth/2);
seg_dist2 = dist2 < (lw/2);
seg_dist3 = dist3 <= (norm(v_pnt)/2);

temp_model(seg_dist1 & seg_dist2 & seg_dist3) = true;

m_st = max([1, 1, 1], lef_but - ext - 1);
m_ed = min(size(out), lef_but - ext + model_size - 2);
t_st = ext + 2 - lef_but + m_st;
t_ed = ext + 2 - lef_but + m_ed;

out(m_st(1):m_ed(1), m_st(2):m_ed(2), m_st(3):m_ed(3)) = ...
    out(m_st(1):m_ed(1), m_st(2):m_ed(2), m_st(3):m_ed(3)) | ...
    temp_model(t_st(1):t_ed(1), t_st(2):t_ed(2), t_st(3):t_ed(3));

end

