function img_reshape = my_reshape_img(img, the_shape_new)
% 通过裁剪的方式，重新设置图像尺寸
% 奇数的对称中心是(n+1)/2，偶数的对称中心是n/2+1
%  Reshape the image by cropping
%  The center of odd numbers is (n+1)/2, and the center of even numbers is n/2+1

data_class = class(img); % 数据类型，single/double等

if numel(the_shape_new)==2 && ndims(img)==2 % 二维情况  2D images
    the_shape = size(img);
    if the_shape(1)==the_shape_new(1) && the_shape(2)==the_shape_new(2)
        img_reshape = img;
        return;
    end

    % 先改第一个维度  1st dimension
    [ll,rr] = my_get_reshape_range(the_shape(1), the_shape_new(1));
    if the_shape_new(1)==the_shape(1)
        img_temp = img;
    elseif the_shape_new(1)<the_shape(1)
        img_temp = img(ll:rr, : );
    else
        img_temp = zeros(the_shape_new(1), the_shape(2), data_class);
        img_temp(ll:rr ,:) = img;
    end
    
    % 再改第二个维度  2nd dimension
    img = img_temp;
    [ll,rr] = my_get_reshape_range(the_shape(2), the_shape_new(2));
    if the_shape_new(2)==the_shape(2)
        img_temp = img;
    elseif the_shape_new(2)<the_shape(2)
        img_temp = img(:, ll:rr );
    else
        img_temp = zeros(the_shape_new(1), the_shape_new(2), data_class);
        img_temp(:, ll:rr) = img;
    end
    
    
    % 结果
    img_reshape = img_temp;
    
else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 三维情况  3D images
    if ndims(img)==2
        img = reshape(img, size(img,1), size(img,2), 1);
    end
    if numel(the_shape_new)==2
        the_shape_new(3) = size(img,3);
    end
    
    the_shape = size(img);
    if the_shape(1)==the_shape_new(1) ...
            && the_shape(2)==the_shape_new(2) ...
            && the_shape(3)==the_shape_new(3)
        img_reshape = img;
        return;
    end
    
    % 先改第一个维度   1st dimension
    [ll1,rr1] = my_get_reshape_range(the_shape(1), the_shape_new(1));
    [ll2,rr2] = my_get_reshape_range(the_shape(2), the_shape_new(2));
    [ll3,rr3] = my_get_reshape_range(the_shape(3), the_shape_new(3));

    if all(the_shape_new<the_shape)
        img_reshape = img(ll1:rr1, ll2:rr2, ll3:rr3);
        return;
    end

    if the_shape_new(1)==the_shape(1)
        img_temp = img;
    elseif the_shape_new(1)<the_shape(1)
        img_temp = img(ll1:rr1, :, : );
    else
        img_temp = zeros(the_shape_new(1), the_shape(2), the_shape(3), data_class);
        img_temp(ll1:rr1 ,:, :) = img;
    end

    % 再改第二个维度  2nd dimension
    img = img_temp;
    if the_shape_new(2)==the_shape(2)
        img_temp = img;
    elseif the_shape_new(2)<the_shape(2)
        img_temp = img(:, ll2:rr2, : );
    else
        img_temp = zeros(the_shape_new(1), the_shape_new(2), the_shape(3), data_class);
        img_temp(:, ll2:rr2, :) = img;
    end
    
    % 再改第三个维度  3rd dimention
    img = img_temp;
    if the_shape_new(3)==the_shape(3)
        img_temp = img;
    elseif the_shape_new(3)<the_shape(3)
        img_temp = img(:, :, ll3:rr3 );
    else
        img_temp = zeros(the_shape_new(1), the_shape_new(2), the_shape_new(3), data_class);
        img_temp(:, :, ll3:rr3) = img;
    end
    
    % 结果 result
    img_reshape = img_temp;

end
    
end
