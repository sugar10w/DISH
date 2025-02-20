function the_image = my_pf__frame(...
    source_list, proj_mat, pixelsize, z_pixelsize, ...
    sample_n, sample_pixelsize, ...
    my_pf_fun, varargin)
% my_project_forward  project forward ͳһ���; ��άת��ά 
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

for i=1:(numel(varargin)-1)
    if ~ischar(varargin{i}) && ~isstring(varargin{i})
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
            if ~iscell(temp) || numel(temp)~=proj_cnt 
                if ~isnan(temp)
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

fprintf('[info] Using %d workers for forward propagation.\n', num_workers);

%% ���� the_intensity_map

if flag_intensity_map
    for i = 1:proj_cnt
        source_list{i} = source_list{i} .* the_intensity_map_list{i};
    end
end

%%

tic

pix_cnt = proj_cnt * sample_n^2 *sample_nz; 
ppm_flag = pix_cnt > 1e9; 

the_image = zeros(sample_n, sample_n, sample_nz, 'single');  % �洢���ս�� final result

time_start = datetime;

if num_workers == 1
    out_num = 0;
    for proj_i = 1:proj_cnt
        
        % ���� projection
        the_image_0 = my_pf_fun( source_list{proj_i} ); 
    
        % ��ת rotate 
        the_image_0r = my_pf__resize_rot(the_image_0, proj_mat{proj_i}, sample_pixelsize, pixelsize, z_pixelsize, ...
            NewShape=[sample_n, sample_n, sample_nz]);
        
        % �޼� reshape
        the_image_sub = my_reshape_img(the_image_0r, [sample_n, sample_n, sample_nz]);
        
        % �ӵ����ս���� add to the final result
        the_image = the_image + the_image_sub;

        if ispc && mod(proj_i, 1)==0 && ppm_flag % ���Ԥ��ʱ�� estimate time
            fprintf(char(ones(1, out_num) * 8));
            out_num = fprintf('[info][pf] proj %d/%d finished, %s\n', proj_i, proj_cnt, my_predict_time(time_start, proj_i / proj_cnt));
        end
    end
else
   
    the_image_sub_list = cell(1, num_workers);  % �洢��ʱ��� temporary results
    
    par_my_pf_fun = cell(1, num_workers);  % ��my_pf_fun ��䵽 cell �У�����㲥���� Avoid broadcasting variables
    for i = 1:num_workers
        par_my_pf_fun{i} = my_pf_fun;
    end
    
    out_num = 0;
    for proj_i = 1:num_workers:proj_cnt
        
        % ȷ����ǰ�������м����index; current indexes to be computed
        par_proj_i_list = proj_i : min(proj_i+num_workers-1, proj_cnt);
        
        % ��ԭ������ȡ Extract from original data
        par_source_list = source_list(par_proj_i_list);
        par_proj_mat    = proj_mat   (par_proj_i_list);
        
        parfor par_proj_i = 1:numel(par_proj_i_list)
            % ���� projection
            the_image_0 = par_my_pf_fun{par_proj_i}( par_source_list{par_proj_i} ); 
            
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
            out_num = fprintf('[info][pf,par] proj %d/%d finished, %s\n', par_proj_i_list(end), proj_cnt, my_predict_time(time_start, (proj_i-1+num_workers)/proj_cnt) );
        end
        
    end
end

toc

% ˥�� attenuation
the_image = the_image .* attenuation_mask;

% �������� adjust energy
the_image = the_image / proj_cnt;  

end

