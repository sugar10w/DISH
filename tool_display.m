% Display source_list

close all; 
 
if exist('X', 'var')
    fprintf('Display X \n');
    source_list_display = X;
elseif exist('source_list_new', 'var')
    fprintf('Display  source_list_new\n');
    source_list_display = source_list_new;         
elseif exist('source_list_stretch', 'var')
    fprintf('Display  source_list_stretch\n');
    source_list_display = source_list_stretch ;
else
    fprintf('Display  source _list\n');
    source_list_display = source_list;
end 

figure;
temp_sum = sum(cat(3, source_list_display{:}), 3);
fprintf('sum=%d\n', sum(temp_sum(:)));
imagesc(temp_sum); colorbar;
% clim([0,1e-6]);


figure;
for i = 1:numel(source_list_display)

    imagesc(source_list_display{i});colorbar;
    
    title(int2str(i))
    
    fprintf('%d\t%f\t%f\n', i, max( source_list_display{i}(:) ), mean( source_list_display{i}(:) ) );
    
    waitforbuttonpress; 
    
end
