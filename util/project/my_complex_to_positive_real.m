function [img_positive_real, best_angle] = my_complex_to_positive_real(img_complex)
%MY_COMPLEX_TO_POSITIVE_REAL 将复数图像转化为正实数图像
% convert complex image to positive real image

the_vector = sum(img_complex, 'all');
if abs(the_vector)==0
    the_vector = 1;
end
the_vector = the_vector / abs(the_vector);

best_angle = angle(the_vector);
img_positive_real = my_complex_to_positive_real_project(img_complex, best_angle);
best_loss = my_get_loss(img_complex, img_positive_real.*exp(1j*best_angle));

the_angle_list = (0:1:360)/180*pi; 
for the_angle = the_angle_list
    img_positive_real_curr = my_complex_to_positive_real_project(img_complex, the_angle);
    the_loss_curr = my_get_loss(img_complex, img_positive_real_curr.*exp(1j*the_angle));
    if the_loss_curr < best_loss
        best_angle = the_angle; 
        best_loss = the_loss_curr;
        img_positive_real = img_positive_real_curr;
    end
end

end

function [img_positive_real] = my_complex_to_positive_real_project(img_complex, the_angle)
the_vector = exp( 1j*the_angle );
img_positive_real = real( img_complex .* conj(the_vector) );
img_positive_real(img_positive_real<0) = 0;
end

function the_loss = my_get_loss(A, B)
C = abs(A-B);
the_loss = C(:)' * C(:);
end