function [MSE] = my_MSE(img1, img2)
%my_MSE ºÚµ•µÿº∆À„ MSE 
% calculate MSE

 MSE = mean( ( img1 - img2 ).^2 , 'all' );


end

