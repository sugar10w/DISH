function y = my_nonlinear_01(x, the_scale)
%MY_NONLINEAR_01 
% custumed sigmoid

if nargin==1
    the_scale = 2;
end

x = x / the_scale;
y = ( 1 ./( 1+exp( (-x+0.5)*5 )  ) - 0.5) / 0.848283639957513 + 0.5;
y = y * the_scale;

end

