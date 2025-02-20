function [img_shell] = model_get_shell(img, thickness)
%MODEL_GET_SHELL 给定三维模型，输出三维模型的壳

if nargin<2
    thickness = 1;
end

kernel = strel('sphere', thickness);

img_erode = imerode(img, kernel);

img_shell = img;
img_shell(img_erode>0) = 0;

end

