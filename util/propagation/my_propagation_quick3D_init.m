function H_list = my_propagation_quick3D_init(N,M, z_list, pixelsize, lambda, method)
% 在使用 my_propagation_quick3D 前使用，针对z生成H_list矩阵列
% Use before my_propagation_quick3D to generate H_list for z_list

x=single(1:M);
y=single(1:N);
L0X=pixelsize*M;
L0Y=pixelsize*N;
k=2*pi/lambda;

u=lambda*(-M/L0X/2+1/L0X*(x-1));
v=lambda*(-N/L0Y/2+1/L0Y*(y-1));
clear x y;
[uu,vv] = meshgrid(u,v); 
clear u v;
temp = uu.^2+vv.^2;
clear uu vv;

H_list = cell(1, numel(z_list));

for i = 1:numel(z_list)
    z = z_list(i);
    if(strcmp(method,'Angular Spectrum'))    
        H_list{i} = exp(1i*k*z*sqrt(1-temp)); % Angular Specturm method 
        H_list{i} = ifftshift(H_list{i});
    elseif(strcmp(method,'Fresnel'))
        H_list{i} = exp(1i*k*z*(1-temp/2)); % Fresnel method
        H_list{i} = ifftshift(H_list{i});
    else
        errordlg('Type of transfer function must be <Angular Spectrum> or <Fresnel>','Error');
    end
end
clear temp;

end