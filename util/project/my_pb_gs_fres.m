function the_source = my_pb_gs_fres(temp_mat, H_list, psf_n, my_propagation_fun, flag_energy)
%MY_PB_GS_FRES 
% GS�������������(temp_mat, the_source)��Ҫ�������������Ǹ���������������������ڼ����½�����
%  GS iteration; input and output (temp_mat, the_source) are both required to be energy/intensity fields rather than complex fields;
%  This function cannot be used to calculate descent direction 

if nargin<5
    flag_energy = 1;  % Ҫ������������� energy correction required
end

psf_z_n = numel(H_list);
temp_z_n = size(temp_mat, 3);

% �޸�temp_mat����ά�ߴ�  Modify the 3D size of temp_mat
temp_mat = my_reshape_img(temp_mat, [psf_n, psf_n] );
temp_mat_raw_energy = temp_mat;  % ����

% ǿ�ȵ��� Intensity/energy adjustment
if flag_energy  % ����ǿ�Ƚ��� energy correction
    temp_mat_no_zero = sum(temp_mat>1e-5, 3);
    temp_mat_no_zero(temp_mat_no_zero<1) = 1;
    for i_ = 1:size(temp_mat, 3)
        temp_mat(:,:,i_) = temp_mat(:,:,i_) ./ (temp_mat_no_zero.^2); 
    end
end
temp_mat(temp_mat<0) = 0;
temp_mat = temp_mat .^ 0.5;   % <!> ������ת��Ϊ����  Convert energy/intensity fields to complex fields

if ispc
    %my_save_complex_tiff(temp_mat, 'temp_mat.tif');
end


% ��ʼ����GS����  GS iteration;
[z_ll, z_rr] = my_get_reshape_range( temp_z_n, psf_z_n);  % ��Ӧ��H_list�е�����  indexes in H_list


% ���ȸ�����ʼֵ��ѡ��ͶӰ��������ƽ������ initialization
if flag_energy
    the_source = sum(temp_mat, 3);
else
    temp_mat_no_zero = sum(temp_mat>1e-5, 3);
    temp_mat_no_zero(temp_mat_no_zero<1) = 1;
    the_source = sum(temp_mat, 3) ./ temp_mat_no_zero; 
end

for iter = 1:10
    % ������õ�ǰ��� current the_image
    the_image = my_propagation_fun(the_source, H_list, z_ll:z_rr);
    if ispc
        the_MSE      = my_MSE     (temp_mat_raw_energy, abs(the_image).^2);
        the_MSE_auto = my_MSE_auto(temp_mat_raw_energy, abs(the_image).^2);
        fprintf('gs %d, MSE = %d, MSE_auto = %d\n', iter, the_MSE, the_MSE_auto );
        %my_save_complex_tiff( the_image, sprintf('image_%02d.tif', iter) );
        %my_save_complex_tiff( the_source, sprintf('image_%02d_r.tif', iter) );
    end

    % ������λ update phase
    temp_mat_ph = temp_mat .*  exp(1i*angle(the_image));
    % �ش� update the_source
    the_source = zeros(psf_n, 'single');
    for i = 1:temp_z_n
        the_source_sub = my_propagation_fun(temp_mat_ph(:,:,i), H_list, z_ll+i-1, -1);
        the_source = the_source + the_source_sub;   % ����ѡ��Ĳ�����ֱ�ӵ��� direct accumulating
    end
    
    % ǿ������λ Convert to real number
    %the_source = abs(the_source);  % ��Ҫֱ��ȡģ��  Do not take the absolute value!
    
    the_source_vector = sum(the_source, 'all');  % ͨ����ȡ�ɷ� Extract the optimal phase
    the_source_vector = the_source_vector / abs(the_source_vector);
    the_source = real(  the_source .* conj(the_source_vector) );
    the_source(the_source<0) = 0;
end

%%
the_source = the_source.^2; % �����ʱ��ǵ�ƽ��ת��Ϊ���� convert to energy/intensity fields 


end

