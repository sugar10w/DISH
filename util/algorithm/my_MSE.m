function [MSE] = my_MSE(img1, img2)
%my_MSE �򵥵ؼ��� MSE 
% calculate MSE

 MSE = mean( ( img1 - img2 ).^2 , 'all' );


end

