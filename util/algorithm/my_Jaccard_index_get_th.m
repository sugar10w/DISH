function [Jaccard_index, img2_th_best] = my_Jaccard_index_get_th(img1, img2, thresh_cof, my_Jaccard_index_fun)
%my_Jaccard_index  穷举阈值，获得最优的 Jaccard index ； 要求img1已经是二值的；
% Enumerate thresholds to obtain the optimal Jaccard index; img1 should be binary;

if nargin<4
    my_Jaccard_index_fun = @(img1, img2)my_Jaccard_index(img1, img2);
end
if nargin<3
    thresh_cof = 1;
end

if thresh_cof < 1
    thresh_cof = 1;
end

if max(img1(:)) ~= 1
    fprintf('[Info] unacceptable img1 (max(img1(:)) ~= 1) . `my_Jaccard_index_auto` may not work   \n');
    %Jaccard_index = my_Jaccard_index_fun(img1, img2);
    %return;
end

img1 = single(img1);
img2 = single(img2);

% img1 = img1 > 0.5;
img1 = imbinarize(img1);

img2_max = max(img2(:));
img2_min = min(img2(:));

Jaccard_index = 0;
img2_th_best = 0;
the_jaccard_region = 1e-4;
the_th_region = 1e-3;

l_Ja = my_Jaccard_index_fun(img1, imbinarize(img2, img2_min));
r_Ja = my_Jaccard_index_fun(img1, imbinarize(img2, img2_max));
m_Ja = my_Jaccard_index_fun(img1, imbinarize(img2, .5*img2_max + .5*img2_min));

th_list_temp = [
    0,  .25,   .5, .75,  1;
    l_Ja, 0, m_Ja,  0, r_Ja;
];
cnt_temp = 0;
while max(th_list_temp(2, :)) - min(th_list_temp(2, :)) > the_jaccard_region ...
    && th_list_temp(1, 5) - th_list_temp(1, 1) > the_th_region
    cnt_temp = cnt_temp + 1;

    th_list_temp(1, 2) = mean(th_list_temp(1, [1, 3]));
    th_list_temp(1, 4) = mean(th_list_temp(1, [3, 5]));
    lm_Ja = my_Jaccard_index_fun(img1, imbinarize(img2, th_list_temp(1, 2)*img2_max + (1-th_list_temp(1, 2))*img2_min));
    mr_Ja = my_Jaccard_index_fun(img1, imbinarize(img2, th_list_temp(1, 4)*img2_max + (1-th_list_temp(1, 4))*img2_min));
    th_list_temp(2, 2) = lm_Ja;
    th_list_temp(2, 4) = mr_Ja;

    [max_Ja, max_ind] = max(th_list_temp(2, :));

    if Jaccard_index < max_Ja
        Jaccard_index = max_Ja;
        img2_th_best = th_list_temp(1, max_ind)*img2_max + (1-th_list_temp(1, max_ind))*img2_min;
    end

    m_ind = min(4, max(2, max_ind));
    th_list_temp(:, 1) = th_list_temp(:, m_ind-1);
    th_list_temp(:, 5) = th_list_temp(:, m_ind+1);
    th_list_temp(:, 3) = th_list_temp(:, m_ind);
end

% for th = 0:0.05:1
%     img2_th = th * img2_max + (1-th) * img2_min;
%     if thresh_cof <= 1
%         img2_temp = img2>img2_th;
%     else
%         img2_temp = (img2 - img2_th) / ((thresh_cof - 1) * img2_th);
%         img2_temp = min(img2_temp, 1);
%         img2_temp = max(img2_temp, 0);
%     end
%     temp = my_Jaccard_index_fun(img1, img2_temp);
% 
%     if Jaccard_index_best<temp
%         Jaccard_index_best = temp;
%         img2_th_best = img2_th;
%     end
% 
% end

% Jaccard_index = Jaccard_index_best;
%fprintf('[Info] my_Jaccard_index_auto img2 > %.2f \n', img2_th_best);

end

