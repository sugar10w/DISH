function model_out = modelRotate(model_in, rot_mat, pnt_ori, pnt_tgt)
%MODELROTATE ��ת��άģ�͡�Rotate a 3D model.
%   ����ԭʼģ�͡�3*3�������4*4��ת����ԭͼê�㡢Ŀ��ê��λ�ã�
%   �����ת���ģ�͡�
%   Input an origin model, 3*3 or augmented 4*4 rotation matrix,
%   anchor point in the origin model, anchor point in the output model,
%   output rotated model.

% �����������
% Manage input args.
if (nargin < 2)
    fprintf('Error:input args not enough\n')
    return;
end
if (nargin < 3)
    pnt_ori = round(size(model_in) / 2);
end
if (nargin < 4)
    pnt_tgt = pnt_ori;
end

% ���ê��
% mark anchor point
model_in = single(model_in ~= 0);
try
    pnt_val = model_in(pnt_ori(1), pnt_ori(2), pnt_ori(3));
    model_in(pnt_ori(1), pnt_ori(2), pnt_ori(3)) = 10000;
catch
    fprintf('Error:input origin anchor point illegal\n');
    model_out = false;
    return;
end

% ��תģ��
% Rotate model
if size(rot_mat) == [3, 3]
    rot_mat = [[rot_mat, [0;0;0]]; 0,0,0,1];
elseif size(rot_mat) == [4, 4]
else
    fprintf('Error:input rotation matrix shape illegal\n');
    return;
end
aff_form = affine3d(rot_mat);
model_rot = imwarp(model_in, aff_form);
[~, ind] = max(model_rot, [], 'all', 'linear');
[px, py, pz] = ind2sub(size(model_rot), ind);
pnt_coord = [py, px, pz];
model_in(pnt_ori(1), pnt_ori(2), pnt_ori(3)) = pnt_val;
model_rot = gather(imwarp(model_in, aff_form) ~= 0);

% ƽ��ģ��
% Translate model
model_res = imtranslate(model_rot, pnt_tgt(:, [2, 1, 3]) - pnt_coord);
model_out = model_res(1:size(model_in, 1), 1:size(model_in, 2),...
                      1:size(model_in, 3));

end

