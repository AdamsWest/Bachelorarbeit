classdef EnvironmentClass
    properties
        time;
        idx = 1;
        rho;
        p;
        T;
        a;
        u_Wg;
        g;
        params;
    end
        methods
         
        function obj = initEnvironmentProperties(obj, n)
            initVar = zeros(n,1);
            obj.time = initVar;
            obj.rho = obj.params.rho_0 * ones(n,1);
            obj.p = initVar;
            obj.T = initVar;
            obj.a = initVar;
            obj.u_Wg = obj.params.u_Wg * ones(n,1);
            obj.g = obj.params.g * ones(n,1);
        end
        
        function obj = pickEnvironmentState(obj, idx)
            obj.time = obj.time(idx);
            obj.rho = obj.rho(idx);
            obj.p = obj.p(idx);
            obj.T = obj.T(idx);
            obj.a = obj.a(idx);
            obj.u_Wg = obj.u_Wg(idx);
            obj.g = obj.g(idx);
        end
         
        function objStore = storeEnvironmentState(objStore, obj, idx)
            objStore.time(idx) = obj.time;
            objStore.rho(idx) = obj.rho;
            objStore.p(idx) = obj.p;
            objStore.T(idx) = obj.T;
            objStore.a(idx) = obj.a;
            objStore.u_Wg(idx) = obj.u_Wg;
            objStore.g(idx) = obj.g;
        end
         
        function obj = incrementEnvironment( obj, dt )
            obj.idx = obj.idx + 1;
            obj.time(obj.idx) = obj.time(max(obj.idx-1,1)) + dt;
        end
        
        function obj = stratopauseConditions(obj)
            % Tropopause
            % Temperature at 11000m altitude 
            obj.params.T_11 = obj.params.T_0 - ...
                0.0065 *(11000-obj.params.H_0);              
            % Air densitiy at 11000m altitude 
            obj.params.rho_11 = obj.params.rho_0 * ...
                (1 - 0.0065*(11000/obj.params.T_0))^4.256;
            % Pressure at 11000m altitude 
            obj.params.p_11 = obj.params.p_0 * ...
                (1 - 0.0065*(11000/obj.params.T_0))^5.256; 
        end
     end
 end
