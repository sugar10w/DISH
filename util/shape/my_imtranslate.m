function img_out = my_imtranslate(img, v, filling)
%my_imtranslate 平移图像（不改变图像大小）
%  Pan the image (without resizing the image)

if nargin==2
    filling=0;
end

[s1, s2] = size(img);

img_out = ones(s1, s2, class(img)) * filling;


v = round(v);


if v(1)>=0
    ro1 = 1+v(1)  :  s1      ;
    r1  = 1       :  s1-v(1) ;
else
    ro1 = 1       :  s1+v(1) ;
    r1  = 1-v(1)  :  s1      ; 
end

if v(2)>=0
    ro2 = 1+v(2)  :  s2      ;
    r2  = 1       :  s2-v(2) ;
else
    ro2 = 1       :  s2+v(2) ;
    r2  = 1-v(2)  :  s2      ; 
end

img_out(ro1, ro2) = img(r1, r2);

end

