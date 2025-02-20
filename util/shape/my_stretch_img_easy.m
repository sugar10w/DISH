function img_new = my_stretch_img_easy(img, the_angle, flag_from_air_to_material, stretch_ratio)

% 简易拉伸
if flag_from_air_to_material
    img_new = my_stretch_img(img, 1/stretch_ratio, the_angle );
else
    img_new = my_stretch_img(img, stretch_ratio, the_angle );
end

end