function [ll, rr] = my_get_reshape_range(shape1, shape2)
%MY_GET_MARGIN 将尺寸从shape1到shape2时，给出合适的对应范围
% 不论放大和缩小多少次，只要回到原始尺寸，坐标就不会变
% 奇数的中心点是(n+1)/2，偶数的中心点是n/2+1
%  When changing the size from shape1 to shape2, give the appropriate corresponding range
%  No matter how you zoom in and out, as long as you return to the original size, the index will not change
%  The center point of odd numbers is (n+1)/2, and the center point of even numbers is n/2+1

if shape2 < shape1 % 缩小 reduce
    mm = floor( (shape1 - shape2 )/2 );
    if mod(shape1,2)==0 && mod(shape2,2)==1  % 仅当偶数变小成奇数时，左边缘+1; 
        mm=mm+1;
    end
    ll = mm+1;
    rr = mm+shape2;
else % 扩大 expand
    mm = floor( (shape2 - shape1)/2 );
    if mod(shape1,2)==1 && mod(shape2,2)==0 % 仅当奇数变大成偶数时，左边缘+1; only if even is expanded to odd
        mm=mm+1;
    end
    ll = mm+1;
    rr = mm+shape1;
end


end

