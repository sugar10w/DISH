classdef my_LossRecorder < handle
    %MY_LOSSRECORDER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        fun_list
        fun_name_list
        output_style_list
        data_list

        fun2_list
        fun2_name_list
        output2_style_list
        data2_list

        round_cnt

        stop_add
    end
    
    methods
        function obj = my_LossRecorder()
            %MY_LOSSRECORDER 构造此类的实例
            obj.fun_list = {};
            obj.fun_name_list = {};
            obj.output_style_list = {};
            obj.data_list = {};

            obj.fun2_list = {};
            obj.fun2_name_list = {};
            obj.output2_style_list = {};
            obj.data2_list = {};

            obj.round_cnt = 0;
            obj.stop_add = 0;
        end
        
        function addfun(obj, fun, fun_name, output_style)
            %METHOD1 此处显示有关此方法的摘要
            if obj.stop_add
                fprintf('[Error] my_LossRecorder has been updated and no function can be added.\n');
                return;
            end

            obj.fun_list{end+1}          = fun;
            obj.fun_name_list{end+1}     = fun_name; 
            obj.output_style_list{end+1} = output_style;
            obj.data_list{end+1} = {};
        end
        function addfun2(obj, fun, fun_name, output_style)
            %METHOD1 此处显示有关此方法的摘要
            if obj.stop_add
                fprintf('[Error] my_LossRecorder has been updated and no function can be added.\n');
                return;
            end

            obj.fun2_list{end+1}          = fun;
            obj.fun2_name_list{end+1}     = fun_name; 
            obj.output2_style_list{end+1} = output_style;
            obj.data2_list{end+1} = {};
        end


        function update(obj, Y, Y_curr)
            obj.stop_add = 1;
            obj.round_cnt = obj.round_cnt + 1;

            for i = 1:numel(obj.fun_list)
                a = obj.fun_list{i}(Y, Y_curr); 
                obj.data_list{i}{obj.round_cnt} = a;
                fprintf(['[%d]%s:\t', obj.output_style_list{i} , '\n'], obj.round_cnt, ...
                    obj.fun_name_list{i}, a);
            end
            
            for i = 1:numel(obj.fun2_list)
                [a,b] =  obj.fun2_list{i}(Y, Y_curr); 
                obj.data2_list{i}{obj.round_cnt} = [a, b];
                fprintf(['[%d]%s:\t', obj.output2_style_list{i} , '\n'], obj.round_cnt, ...
                    obj.fun2_name_list{i}, a, b);
            end
        end

        function review(obj)
            
            for i = 1:numel(obj.fun_list)
                fprintf('%s\t', obj.fun_name_list{i});
            end
            for i = 1:numel(obj.fun2_list)
                fprintf('%s\t', obj.fun2_name_list{i});
            end
            fprintf('\n');

            for j = 1:obj.round_cnt
                fprintf('[%d]\t', j);
                for i = 1:numel(obj.fun_list)
                    fprintf([obj.output_style_list{i}, '\t'], obj.data_list{i}{j});
                end
                for i = 1:numel(obj.fun2_list)
                    fprintf([obj.output2_style_list{i}, '\t'], ...
                        obj.data2_list{i}{j}(1), obj.data2_list{i}{j}(2));
                end
                fprintf('\n');
            end

        end


    end
end

