function [the_sample_rotate] = my_rotate_2(the_sample, proj_mat, bbox, RB2)
%MY_ROTATE 旋转三维模型 Rotate 3D volume
%  这个函数严格地保证中心位置无偏移（且边长奇偶性不变）
%  This function strictly guarantees that the center is not offset (and the parity of the side lengths is unchanged)

if nargin<3
    bbox = 'crop';
end
if nargin<4
    RB2 = NaN;
end

proj_mat4 = [ proj_mat, [0;0;0]; 0,0,0,1];
tform = affine3d(proj_mat4);

[ny, nx, nz] = size(the_sample);  
the_parity = mod(nx,2);   % 这里假定nx,ny,nz的奇偶性是一样的；It is assumed that the parity of nx, ny, and nz is the same;

RA = imref3d([ny,nx,nz], [-nx/2,nx/2], [-ny/2,ny/2], [-nz/2,nz/2]);
% if the_parity==1
%     RA = imref3d([ny,nx,nz], [-nx/2,nx/2], [-ny/2,ny/2], [-nz/2,nz/2]);
% else
%     RA = imref3d([ny,nx,nz], [-nx/2,nx/2]-0.5, [-ny/2,ny/2]-0.5, [-nz/2,nz/2]-0.5);
% end

switch bbox
    case 'loose'
        RB1 = images.spatialref.internal.applyGeometricTransformToSpatialRef(RA, tform);
        nx_r = RB1.ImageExtentInWorldX;
        ny_r = RB1.ImageExtentInWorldY;
        nz_r = RB1.ImageExtentInWorldZ;
        if mod(nx_r,2)~=the_parity
            nx_r = nx_r + 1;
        end
        if mod(ny_r,2)~=the_parity
            ny_r = ny_r + 1;
        end
        if mod(nz_r,2)~=the_parity
            nz_r = nz_r + 1;
        end

        RB2 = imref3d([ny_r,nx_r,nz_r], [-nx_r/2,nx_r/2], [-ny_r/2,ny_r/2], [-nz_r/2,nz_r/2]);
        % if the_parity==1
        %     RB2 = imref3d([ny_r,nx_r,nz_r], [-nx_r/2,nx_r/2], [-ny_r/2,ny_r/2], [-nz_r/2,nz_r/2]);
        % else
        %     RB2 = imref3d([ny_r,nx_r,nz_r], [-nx_r/2,nx_r/2]-0.5, [-ny_r/2,ny_r/2]-0.5, [-nz_r/2,nz_r/2]-0.5);
        % end
    case 'crop'
        RB2 = RA;
    case 'RB2'
        %fprintf('[Info][my_rotate_2] Use the specified size of RB2 \n');
    otherwise
        fprintf(2, '[Error][my_rotate_2] Unknown bbox `%s`\n', bbox);
end

[the_sample_rotate, ~] = imwarp(the_sample, RA, tform, "linear", 'OutputView',RB2);

end