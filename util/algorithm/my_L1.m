function the_l1 = my_L1(img1, img2, lower_th, upper_th)

if nargin<4
    upper_th = 1;
end
if nargin<3
    lower_th = 0.5;
end

img1(img1>=upper_th) = upper_th;
img1(img1<=lower_th) = lower_th;

img2(img2>=upper_th & img1==upper_th) = upper_th;
img2(img2<=lower_th & img1==lower_th) = lower_th;

the_l1 = mean(abs(img1-img2), "all");

end