function the_stack = my_propagation_quick3D_GPU(U0, H_list, z_list, the_direction)
%MY_PROPAGATION_QUICK2D_GPU 传播z_lsit中标注的z轴位置；使用前先检查轴向分辨率
%  Propagating the z-axis position annotated in z_lsit; check the axial resolution before use

% 缺省参数 Default parameter
if nargin<4
    the_direction = 1;
end
if nargin<3
    z_list = 1:numel(H_list);
end

% 转成GPU格式 Convert to gpuArray
if ~isa(U0, 'gpuArray')
    U0 = gpuArray(U0);
end
FU0=fft2(U0);
[N, M] = size(U0);

clear U0;

if the_direction == 1
    % 卷积(FFT乘) convolution (multiply on FFT)
    if numel(z_list)==1
        the_stack = gather( ifft2( FU0 .* H_list{z_list} ));
    else
        the_stack = zeros(N,M, numel(z_list), 'single');
        for i = 1:numel(z_list)
            the_stack(:,:,i) = gather( ifft2( FU0 .* H_list{z_list(i)} ) );
        end
    end
else
    % 反方向，使用FFT除 deconvolution (divide on FFT)
    if numel(z_list)==1
        the_stack = gather( ifft2( FU0 ./ H_list{z_list} ));
    else
        the_stack = zeros(N,M, numel(z_list), 'single');
        for i = 1:numel(z_list)
            the_stack(:,:,i) = gather( ifft2( FU0 ./ H_list{z_list(i)} ) );
        end
    end
end


end

