function the_image = my_pf__frame_s(...
    source_list, proj_mat, pixelsize, z_pixelsize, ...
    sample_n, sample_pixelsize, ...
    my_pf_fun, varargin)
% my_project_forward  project forward �޸��˵�ͳһ��ܣ����ڴ������������; ��άת��ά 
%  project forward general framework (3D -> 2D)
%   ����Ҫ��ȫ1��ǰ����Ҫ���ȫ1�������ж�����
%   Additional requirement: images of all 1 -> volume of all 1, which is convenient to compare the energy

proj_cnt = numel(source_list);

%%
par = parcluster('local');
max_num_workers = par.NumWorkers;
num_workers = NaN;

flag_intensity_map = 0;
the_intensity_map_list = NaN; 

sample_nz = sample_n;
attenuation_mask = NaN;

flag_reference_image = 0;
reference_image = NaN;

for i=1:(numel(varargin)-1)
    if ~isstring(varargin{i}) && ~ischar(varargin{i})
        continue;
    end
    switch varargin{i}
        case 'NumWorkers'
            temp = varargin{i+1};
            if isnumeric(temp)
                if temp > 0 && temp < max_num_workers
                    num_workers = temp;
                else
                    num_workers = max_num_workers;
                end
            end
        case 'IntensityMap'
            temp = varargin{i+1};
            if ~iscell(temp) || mod(proj_cnt, numel(temp))~=0
                if iscell(temp) || isnumeric(temp)
                    fprintf(2,'[Error] IntensityMap must be a cell of `proj_cnt` images.\n');
                end
            else
                flag_intensity_map = 1;
                the_intensity_map_list = temp;
            end
        case 'SampleNZ'
            temp = varargin{i+1};
            if isnumeric(temp)
                sample_nz = temp;
                if isnan(attenuation_mask)
                    attenuation_mask = ones(1,1,sample_nz);
                elseif numel(attenuation_mask) ~= sample_nz
                    fprintf(2, '[Error] numel(temp) ~= sample_nz\n');
                end
            end
        case 'AttenuationMask'
            temp = varargin{i+1};
            if isnumeric(temp) 
                sample_nz = numel(temp);
                attenuation_mask = ones(1,1,sample_nz);
                attenuation_mask(1,1,:) = temp;
            end
        case 'ReferenceImage'
            temp = varargin{i+1};
            flag_reference_image = 1;
            reference_image = temp;
        otherwise
            continue
    end
end

if numel(attenuation_mask)==1 && isnan(attenuation_mask)
    attenuation_mask = ones(1,1,sample_nz);
end

%% ��ò��в������ڴ�ռ����ޣ�û���򵥵���parfor�� get parallel parameters; memory space is limited, thus `parfor` cannot be used simply;
if num_workers==1
    fprintf('[Info][my_pf__frame] Since `num_workers==1`, `parfor` would not be utilized in `my_pf__frame` function. \n');
elseif isempty(gcp("nocreate"))
    if ~isnan(num_workers)
        parpool('local', num_workers);
    else
        parpool('local');
        num_workers = max_num_workers;
    end
else
    temp = gcp("nocreate");
    now_num_workers = temp.NumWorkers;
    
    if ~isnan(num_workers) && now_num_workers ~= num_workers
        delete(gcp("nocreate"));
        parpool('local', num_workers);
    else
        num_workers = now_num_workers;
    end
end
par_num = num_workers;

fprintf('[info][pf_frame_s] Using %d workers for forward propagation.\n', num_workers);

%% ����source_list�е�Ԫ��������proj_mat����������
% ����my_binarize_multiply_source_list���������ϵ
proj_cnt_j = numel(proj_mat);

proj_i_to_j = 1:proj_cnt;
if proj_cnt ~= proj_cnt_j
    fprintf('[info][pf_frame_s] numel(source_list) ~= proj_cnt_j \n');
    if mod(proj_cnt, proj_cnt_j)~=0
        fprintf(2, '[Error][pf_frame_s] mod(proj_cnt, proj_cnt_j)~=0 \n');
    end

    for proj_i = 1:proj_cnt
        proj_j = floor((proj_i-1)/proj_cnt*proj_cnt_j+0.4999)+1;
        if proj_j>proj_cnt_j
            proj_j = 1;
        end
        proj_i_to_j(proj_i) = proj_j;
    end
end

%% ���� the_intensity_map

if flag_intensity_map
    for proj_i = 1:proj_cnt
        proj_j = proj_i_to_j(proj_i);
        source_list{proj_i} = source_list{proj_i} .* the_intensity_map_list{proj_j};
    end
end


%%
tic

pix_cnt = proj_cnt * sample_n^2*sample_nz; 
ppm_flag = pix_cnt > 1e9; 

the_image = zeros(sample_n, sample_n, sample_nz, 'single');  % �洢���ս�� final result

time_start = datetime;


if num_workers == 1
    out_num = 0;
    for proj_i = 1:proj_cnt
        
        proj_j = proj_i_to_j(proj_i);
        proj_i_curr_ratio = (proj_j-1)/proj_cnt_j;

        % ���� projection
        the_image_0 = my_pf_fun( source_list{proj_i}, proj_i_curr_ratio  ); 
    
        % ��ת rotate 
        the_image_0r = my_pf__resize_rot(the_image_0, proj_mat{proj_j}, sample_pixelsize, pixelsize, z_pixelsize, ...
            NewShape=[sample_n, sample_n, sample_nz]);
    
        % �޼� reshape
        the_image_sub = my_reshape_img(the_image_0r, [sample_n, sample_n, sample_nz]);

        % ���۵�ǰ�ӽ�
        if flag_reference_image
            the_score = sum( the_image_sub.*reference_image, 'all') / sum(the_image_sub, 'all');
            fprintf('%d\t%d\n', proj_i, the_score);
        end

        % �ӵ����ս���� add to the final result
        the_image = the_image + the_image_sub;

        if ispc && mod(proj_i, 1)==0 && ppm_flag && ~flag_reference_image % ���Ԥ��ʱ�� estimate time
            fprintf(char(ones(1, out_num) * 8));
            out_num = fprintf('[info][pf_s] proj %d/%d finished, %s\n', proj_i, proj_cnt, my_predict_time(time_start, proj_i / proj_cnt));
        end
    end
else
    
    the_image_sub_list = cell(1, par_num);  % �洢��ʱ��� temporary results
    
    par_my_pf_fun = cell(1, par_num);  % ��my_pf_fun ��䵽 cell �У�����㲥���� Avoid broadcasting variables
    for i = 1:par_num
        par_my_pf_fun{i} = my_pf_fun;
    end

    out_num = 0;
    for proj_i = 1:par_num:proj_cnt
        
        % ȷ����ǰ�������м����index; current indexes to be computed
        par_proj_i_list = proj_i : min(proj_i+par_num-1, proj_cnt);
        
        % ��ԭ������ȡ Extract from original data
        par_source_list = source_list(par_proj_i_list);

        par_proj_j_list = proj_i_to_j(par_proj_i_list);
        par_proj_mat    = proj_mat   (par_proj_j_list);

        parfor par_proj_i = 1:numel(par_proj_i_list)
            
            proj_i_curr = par_proj_i_list(par_proj_i);
            proj_i_curr_ratio = (par_proj_j_list(par_proj_i)-1)/proj_cnt_j;
            
            % ���� projection
            % ͬʱָ����ǰ��λ�� specify the current location
            the_image_0 = par_my_pf_fun{par_proj_i}( par_source_list{par_proj_i}, proj_i_curr_ratio ); 
    
            % ��ת rotate 
            the_image_0r = my_pf__resize_rot(the_image_0, par_proj_mat{par_proj_i}, sample_pixelsize, pixelsize, z_pixelsize, ...
                NewShape=[sample_n, sample_n, sample_nz]);
    
            % �޼� reshape
            the_image_sub = my_reshape_img(the_image_0r, [sample_n, sample_n, sample_nz]);
            the_image_sub_list{par_proj_i} = the_image_sub;
        end
        
        % �ӵ����ս���� add to the final result
        for i = 1:numel(par_proj_i_list)  
            the_image = the_image + the_image_sub_list{i};
        end
        
        if ispc && ppm_flag % ���Ԥ��ʱ�� estimate time
            fprintf(char(ones(1, out_num) * 8));
            out_num = fprintf('[info][pf_s,par] proj %d/%d finished, %s\n', par_proj_i_list(end), proj_cnt, my_predict_time(time_start, (proj_i-1+par_num)/proj_cnt) );
        end
        
        
    end
end

toc

% ˥�� attenuation
the_image = the_image .* attenuation_mask;

% �������� adjust energy
the_image = the_image / proj_cnt;  

end

