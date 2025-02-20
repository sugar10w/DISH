function FF_stack = my_pf_incoherent_2(the_source, F1_stack, z_H_list)
%MY_PF_INCOHERENT_2 非相干传播的仿真

N = size(the_source, 1);

% 卷积 convolution
test_img = my_conv2_fft(the_source, F1_stack(:,:,1));
[img_size_1, img_size_2] = size(test_img);
FF_stack = zeros(img_size_1, img_size_2, numel(z_H_list));
for i = 1:numel(z_H_list)
    FF_stack(:,:,i) = my_conv2_fft(the_source, F1_stack(:,:,i));
end
FF_stack = my_reshape_img(FF_stack, [N,N]);

end

