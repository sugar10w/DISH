function my_save_complex_tiff(mat, filename_output_stack, flag, data_type)
% 把一个复数三维矩阵保存到tif；一个存强度，一个存相位
% save a complex 3D matrix to tif; One for strength, one for phase

if nargin <=2
    flag = 'max';
end

if nargin<=3
    data_type = 16;
end

if data_type == 8
    data_max = 255;
    data_max_1 = 2.55;
    data_mid = 128;
    data_10  = 100;
else
    data_max = 65535;
    data_max_1 = 6.5535;
    data_mid = 32768;
    data_10  = 10000;
end

if ~isreal(mat)  % 保存三维复数矩阵 3D complex matrix
    mat_proj = abs(mat).^2;
    mat_angle =  uint16( data_mid + angle(mat)/pi*(data_mid-1) );
    
    mmax = max(mat_proj(:));
    switch flag
        case 'max'
            % 不操作 no operation
        case '10'
            mmax = 10^ceil( log(mmax/data_max_1)/log(10) ) / data_10 * data_max ;
        otherwise
            mmax = 1;
    end
    mat_proj_norm = uint16(data_max*mat_proj /mmax);

    options.overwrite = true;

    if data_type==8
        mat_proj_norm = uint8(mat_proj_norm);
        mat_angle = uint8(mat_angle);
    end
    saveastiff(mat_proj_norm, filename_output_stack, options);
    saveastiff(mat_angle,     [filename_output_stack, '.angle.tif'], options);
    
else % 保存三维实数矩阵 3D real matrix
    mat = single(mat);
    mmax = max(mat(:));
    
    switch flag
        case 'max'
            % 不操作 no operation
        case '10'
            mmax = 10^ceil( log(mmax/data_max_1)/log(10) ) / data_10 * data_max ;
        otherwise
            mmax = 1;
    end
    
    if min(mat(:))>mmax*-0.01
        psf_norm = uint16(data_max*mat/mmax);
    else
        psf_norm = data_max*mat/mmax;
        psf_norm = uint16(psf_norm/2 + data_mid);
    end
    
    options.overwrite = true;
    
    if data_type==8
        psf_norm = uint8(psf_norm);
    end
    saveastiff(psf_norm, filename_output_stack, options);
end
    
end
