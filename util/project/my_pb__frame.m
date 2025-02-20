function source_list = my_pb__frame(...
    the_sample, sample_pixelsize,  ...
    proj_mat, pixelsize, z_pixelsize, ...
    my_pb_fun, varargin)
%%
par = parcluster('local');
max_num_workers = par.NumWorkers;
num_workers = 1;
the_eps = 0;

for i=1:(numel(varargin)-1)
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
        case 'eps'
            temp = varargin{i+1};
            if isnumeric(temp)
                the_eps = temp;
            end
        otherwise
            continue
    end
end

%% 获得并行参数；内存空间有限，没法简单地用parfor； get parallel parameters; memory space is limited, thus `parfor` cannot be used simply;
if num_workers==1
    fprintf('[Info][my_pb__frame] Since `num_workers==1`, `parfor` would not be utilized in `my_pb__frame` function. \n');
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

fprintf('[Info] Using %d workers for backward propagation.\n', num_workers);

%% MY_PB_FRAME  project backward 统一框架（三维转二维）
%  project backward general framework (3D -> 2D)

the_sample(the_sample<the_eps) = the_eps;

proj_cnt = numel(proj_mat);
source_list = cell(1, proj_cnt);

pix_cnt = proj_cnt*numel(the_sample);
ppm_flag = pix_cnt > 1e10;

%%
tic

if num_workers==1
    time_start = datetime;
    out_num = 0;
    for proj_i = 1:proj_cnt
        temp_mat = my_pb__rot_resize(the_sample, proj_mat{proj_i}, sample_pixelsize, pixelsize, z_pixelsize); % 旋转
        source_list{proj_i} = my_pb_fun(temp_mat); % 传播
        if ispc && mod(proj_i,1)==0
            fprintf(char(ones(1, out_num) * 8));
            out_num = fprintf('[info][pb] projection %d/%d, %s \n', proj_i, numel(proj_mat), my_predict_time(time_start, proj_i/proj_cnt));
        end
    end 
else

    if ppm_flag
        warning('off','instrument:udp:ClassToBeRemoved');
        ppm = ParforProgressbar(proj_cnt);
    else
        ppm = 0;
    end

    parfor proj_i = 1:proj_cnt
        warning('off','instrument:udp:ClassToBeRemoved');
        temp_mat = my_pb__rot_resize(the_sample, proj_mat{proj_i}, sample_pixelsize, pixelsize, z_pixelsize); % 旋转
        source_list{proj_i} = my_pb_fun(temp_mat); % 传播
        if ppm_flag
            ppm.increment();
        end
    end

    if ppm_flag
        delete(ppm);
    end
end
toc

warning('on');
end