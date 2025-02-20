function V = my_clean_margin(V, mm, a)
%MY_CLEAN_MARGIN 将边界处设置成指定的值
%  Set the boundary to the specified value

if ndims(V)==2
    [n1,n2] = size(V);
    V (1:mm,       :) = a;
    V (n1-mm+1:n1, :) = a;
    V (:,       1:mm) = a;
    V (:, n2-mm+1:n2) = a;
elseif ndims(V)==3
    [n1,n2,n3] = size(V);
    V (1:mm,       :, :) = a;
    V (n1-mm+1:n1, :, :) = a;
    V (:,       1:mm, :) = a;
    V (:, n2-mm+1:n2, :) = a;
    V (:,          :,       1:mm) = a;
    V (:,          :, n3-mm+1:n3) = a;
else
    fprintf('[Error] [clean margin] ndims(V)=%d \n', ndims(V));
end

end

