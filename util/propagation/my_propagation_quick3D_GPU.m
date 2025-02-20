function the_stack = my_propagation_quick3D_GPU(U0, H_list, z_list, the_direction)
%MY_PROPAGATION_QUICK2D_GPU ����z_lsit�б�ע��z��λ�ã�ʹ��ǰ�ȼ������ֱ���
%  Propagating the z-axis position annotated in z_lsit; check the axial resolution before use

% ȱʡ���� Default parameter
if nargin<4
    the_direction = 1;
end
if nargin<3
    z_list = 1:numel(H_list);
end

% ת��GPU��ʽ Convert to gpuArray
if ~isa(U0, 'gpuArray')
    U0 = gpuArray(U0);
end
FU0=fft2(U0);
[N, M] = size(U0);

clear U0;

if the_direction == 1
    % ���(FFT��) convolution (multiply on FFT)
    if numel(z_list)==1
        the_stack = gather( ifft2( FU0 .* H_list{z_list} ));
    else
        the_stack = zeros(N,M, numel(z_list), 'single');
        for i = 1:numel(z_list)
            the_stack(:,:,i) = gather( ifft2( FU0 .* H_list{z_list(i)} ) );
        end
    end
else
    % ������ʹ��FFT�� deconvolution (divide on FFT)
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

