function [ll, rr] = my_get_reshape_range(shape1, shape2)
%MY_GET_MARGIN ���ߴ��shape1��shape2ʱ���������ʵĶ�Ӧ��Χ
% ���۷Ŵ����С���ٴΣ�ֻҪ�ص�ԭʼ�ߴ磬����Ͳ����
% ���������ĵ���(n+1)/2��ż�������ĵ���n/2+1
%  When changing the size from shape1 to shape2, give the appropriate corresponding range
%  No matter how you zoom in and out, as long as you return to the original size, the index will not change
%  The center point of odd numbers is (n+1)/2, and the center point of even numbers is n/2+1

if shape2 < shape1 % ��С reduce
    mm = floor( (shape1 - shape2 )/2 );
    if mod(shape1,2)==0 && mod(shape2,2)==1  % ����ż����С������ʱ�����Ե+1; 
        mm=mm+1;
    end
    ll = mm+1;
    rr = mm+shape2;
else % ���� expand
    mm = floor( (shape2 - shape1)/2 );
    if mod(shape1,2)==1 && mod(shape2,2)==0 % ������������ż��ʱ�����Ե+1; only if even is expanded to odd
        mm=mm+1;
    end
    ll = mm+1;
    rr = mm+shape1;
end


end

