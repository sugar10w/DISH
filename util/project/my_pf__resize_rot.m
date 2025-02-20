function the_image_0r = my_pf__resize_rot(the_image_0, proj_mat, ...
    sample_pixelsize, pixelsize, z_pixelsize, ...
    varargin)
%MY_PF_RESIZE_ROT 在正向传播(二维转三维)时用到的小函数
%  used in forward projection (2D -> 3D)

flag_new_shape = 0;
the_new_shape = NaN;

for i=1:(numel(varargin)-1)
    if ~isstring(varargin{i}) && ~ischar(varargin{i})
        continue;
    end
    switch varargin{i}
        case 'NewShape'
            temp = varargin{i+1};
            if isnumeric(temp)
                flag_new_shape = 1;
                the_new_shape=temp;
            end
    end
end

% the new shape
if flag_new_shape
    temp_n_1 = the_new_shape(1);
    temp_n_2 = the_new_shape(2);
    temp_z_n = the_new_shape(3);
    the_parity = mod(temp_n_1,2); 
else
    temp_n_1 = floor( size(the_image_0, 1)   * pixelsize / sample_pixelsize );
    temp_n_2 = floor( size(the_image_0, 2)   * pixelsize / sample_pixelsize );
    temp_z_n = floor( size(the_image_0, 3) * z_pixelsize / sample_pixelsize );
    the_parity = mod(size(the_image_0, 1),2); 
    if mod(temp_n_1, 2) ~= the_parity
        temp_n_1 = temp_n_1 + 1;
    end
    if mod(temp_n_2, 2) ~= the_parity
        temp_n_2 = temp_n_2 + 1;
    end
    if mod(temp_z_n, 2) ~= the_parity
        temp_z_n = temp_z_n + 1;
    end
end

RB2 = imref3d([temp_n_1,temp_n_2,temp_z_n], [-temp_n_1/2,temp_n_1/2], ...
    [-temp_n_2/2,temp_n_2/2], [-temp_z_n/2,temp_z_n/2]);
% if the_parity==1
%     RB2 = imref3d([temp_n_1,temp_n_2,temp_z_n], [-temp_n_1/2,temp_n_1/2], ...
%         [-temp_n_2/2,temp_n_2/2], [-temp_z_n/2,temp_z_n/2]);
% else
%     RB2 = imref3d([temp_n_1,temp_n_2,temp_z_n], [-temp_n_1/2,temp_n_1/2]-0.5, ...
%         [-temp_n_2/2,temp_n_2/2]-0.5, [-temp_z_n/2,temp_z_n/2]-0.5);
% end


% % resize
% the_image_0s = imresize3(the_image_0, [temp_n_1, temp_n_2, temp_z_n], 'linear');

% transformation matrix
the_matrix = [ ...
    pixelsize/sample_pixelsize 0 0; ...
    0 pixelsize/sample_pixelsize 0; ...
    0 0 z_pixelsize/sample_pixelsize ] * inv(proj_mat);
% rotate
the_image_0r = my_rotate_2( the_image_0, the_matrix, 'RB2', RB2);

end

