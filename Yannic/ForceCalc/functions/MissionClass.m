classdef MissionClass
     properties
        time;
        idx = 1;
        gamma;
        V_K;
        distance;
        altitude;
%         configuration;
        i_configuration
        params;
     end
     methods
         
        function obj = initMissionProperties(obj, n)
            initVar = zeros(n,1);
            obj.time = initVar;
            obj.gamma = initVar;
            obj.V_K = initVar;
            obj.distance = initVar;
            obj.altitude = obj.params.H_0 * ones(n,1);
%             obj.configuration = cell(n,1);
%             obj.configuration{1} = obj.params.configuration{1};
            obj.i_configuration = ones(n,1);
        end
        
%         function obj = computeNumberIterations(obj)
%             obj.params.n_iter = 1;
%             for i = 1:length(obj.params.opt_param_range(:,1))
%                 n_iter_tmp = length(obj.params.opt_param_range(i,1):...
%                     obj.params.opt_param_delta(i):...
%                     obj.params.opt_param_range(i,2));
%                 obj.params.n_iter = obj.params.n_iter * n_iter_tmp;
%             end
%         end
        
        function obj = pickMissionState(obj, idx)
            obj.time = obj.time(idx);
            obj.gamma = obj.gamma(idx);
            obj.V_K = obj.V_K(idx);
            obj.distance = obj.distance(idx);
            obj.altitude = obj.altitude(idx);
            obj.i_configuration = obj.i_configuration(idx);
        end
         
        function objStore = storeMissionState(objStore, obj, idx)
            objStore.time(idx) = obj.time;
            objStore.gamma(idx) = obj.gamma;
            objStore.V_K(idx) = obj.V_K;
            objStore.distance(idx) = obj.distance;
            objStore.altitude(idx) = obj.altitude;
            objStore.i_configuration(idx) = obj.i_configuration;
        end
        
        function dt = getFlightTime( obj )
            if obj.params.isDynamic
                dt = obj.params.dt;
            else
                % determine flight time
                switch obj.params.configuration{obj.i_configuration(obj.idx)}
                    case 'loiter'
                        dt = obj.params.paramForConfig(obj.i_configuration(obj.idx));
                    case 'climb'
                        dt = obj.params.paramForConfig(obj.i_configuration(obj.idx)) / ...
                            ( obj.V_K(obj.idx) * sin(obj.gamma(obj.idx)) );
                    case 'cruise'
                        dt = obj.params.paramForConfig(obj.i_configuration(obj.idx)) / ...
                            ( obj.V_K(obj.idx) * cos(obj.gamma(obj.idx)) );
                end
            end
        end
        
        function obj = incrementMission( obj, dt )
            obj.idx = obj.idx + 1;
            obj.time(obj.idx) = obj.time(max(obj.idx-1,1)) + dt;

        end
        
        function obj = integratePosition(obj)
            dt = obj.time( obj.idx ) - obj.time( max( obj.idx - 1, 1 ) );
            obj.distance(obj.idx) = obj.V_K(obj.idx-1) * cos(obj.gamma(obj.idx-1)) * ...
                dt + obj.distance( max( obj.idx - 1, 1 ) );
            obj.altitude(obj.idx) = obj.V_K(obj.idx-1) * sin(obj.gamma(obj.idx-1)) * ...
                dt + obj.altitude( max( obj.idx - 1, 1 ) );
        end
        
        function obj = incrementConfiguration(obj)
            % while previous configuration is equal to the current
            % configurationlength(obj.time)
            for i = 1:obj.idx-1
                if obj.i_configuration( max( obj.idx-i, 1 ) ) ~= ...
                        obj.i_configuration(obj.idx)
                    break;
                end
            end
            % compute the indicator (condition) if the configuration can be incremented
            switch obj.params.configuration{obj.i_configuration(obj.idx)}
                case 'loiter'
                    % determine flight time
                    indicator = obj.time(obj.idx) - obj.time(max(obj.idx-i,1));
                case 'climb'
                    % determine altitude difference
                    indicator = obj.altitude(obj.idx) - ...
                        obj.altitude(max(obj.idx-i,1));
                case 'cruise'
                    % determine distance difference
                    indicator = obj.distance(obj.idx) - ...
                        obj.distance(max(obj.idx-i,1));
            end
            % increment the current configuration
            if indicator >= obj.params.paramForConfig(obj.i_configuration(obj.idx-1)) 
                if obj.i_configuration(obj.idx-1) < length(obj.params.configuration)
                    % increment
                    obj.i_configuration(obj.idx:end) = obj.i_configuration(obj.idx-1) + 1;
                else
                    % indicate that mission is completed
                    obj.i_configuration(obj.idx:end) = 0;
                end
            end
        end

     end
 end
