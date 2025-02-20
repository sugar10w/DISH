function [the_intensity] = my_get_re_modulate_intensity(proj_cnt, mode1, mode2)
% 配合my_3d_re_modulate使用
% use with my_3d_re_modulate

if nargin<3
    mode2 = 'none';
end

the_intensity1 = ones(1, proj_cnt);
switch mode1
    case 'line1'
        ct1 = proj_cnt * 0.25 + 1;
        ct2 = proj_cnt * 0.75 + 1;
        for i = 1:proj_cnt
            the_intensity1(i) = 1 - min( abs(i-ct1), abs(i-ct2) ) / (proj_cnt*0.25); 
        end
        the_intensity1(the_intensity1<0.2) = 0.2;
   case 'line1 half1'
        ct1 = proj_cnt * 0.25 + 1;
        ct2 = proj_cnt * 0.75 + 1;
        for i = 1:proj_cnt
            the_intensity1(i) = 1 - min( abs(i-ct1), abs(i-ct2) ) / (proj_cnt*0.25); 
        end
        the_intensity1(the_intensity1<0.2) = 0.2;
        the_intensity1(round(proj_cnt/2)+1:proj_cnt) = 0;
   case 'line1 half2'
        ct1 = proj_cnt * 0.25 + 1;
        ct2 = proj_cnt * 0.75 + 1;
        for i = 1:proj_cnt
            the_intensity1(i) = 1 - min( abs(i-ct1), abs(i-ct2) ) / (proj_cnt*0.25); 
        end
        the_intensity1(the_intensity1<0.2) = 0.2;
        the_intensity1(1:round(proj_cnt/2)) = 0;    
        
    case 'line2'
        ct1 = proj_cnt * 0.25 + 1;
        ct2 = proj_cnt * 0.75 + 1;
        for i = 1:proj_cnt
            the_intensity1(i) = min( abs(i-ct1), abs(i-ct2) ) / (proj_cnt*0.25); 
        end
        the_intensity1(the_intensity1<0.2) = 0.2;
    case 'line2 half1'
        ct1 = proj_cnt * 0.25 + 1;
        ct2 = proj_cnt * 0.75 + 1;
        for i = 1:proj_cnt
            the_intensity1(i) = min( abs(i-ct1), abs(i-ct2) ) / (proj_cnt*0.25); 
        end
        the_intensity1(the_intensity1<0.2) = 0.2;
        the_intensity1(ct1:ct2) = 0;
    case 'line2 half2'
        ct1 = proj_cnt * 0.25 + 1;
        ct2 = proj_cnt * 0.75 + 1;
        for i = 1:proj_cnt
            the_intensity1(i) = min( abs(i-ct1), abs(i-ct2) ) / (proj_cnt*0.25); 
        end
        the_intensity1(the_intensity1<0.2) = 0.2;
        the_intensity1(1:ct1) = 0;   
        the_intensity1(ct2:proj_cnt) = 0;

    case 'cos^2'
        for i = 1:proj_cnt
            the_intensity1(i) = (  cos (i/proj_cnt*2*pi)  ) ^2 ;
        end
    case 'sin^2'
        for i = 1:proj_cnt
            the_intensity1(i) = (  sin (i/proj_cnt*2*pi)  ) ^2 ;
        end
    case 'sin'
        for i = 1:proj_cnt
            the_intensity1(i) = sin (i/proj_cnt*2*pi) ;
        end
    otherwise
        for i = 1:proj_cnt
            the_intensity1(i) = 1; % 不调制 do not modulate
        end
        fprintf('[info][my_get_re_modulate_intensity] no modulation\n');
end


% 若有必要，使用二次调制  secondary modulation if necessary
the_intensity2 = ones(1, proj_cnt);
switch mode2
    case 'ci1'
        u = 0.8;
        for i = 1:proj_cnt
            xi = (i - 1) / proj_cnt; 
            yi = ( min(0.25, abs(xi-0.5)) - 0.25 ) / 0.25 ;
            the_intensity2(i) = 1 + (1-u)*yi;
        end
end

% 合并 merge
the_intensity = the_intensity1 .* the_intensity2; 

end

