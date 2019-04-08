classdef MotorClass
     properties
         time;
         idx = 1;
        U;
        I;
        tau;
        Omega;
        eta;
        workload;
        params;
     end
     methods
         
        function obj = initMotorProperties(obj, n)
            initVar = zeros(n,1);
            obj.time = initVar;
            obj.U = initVar;
            obj.I = initVar;
            obj.tau = initVar;
            obj.Omega = initVar;
            obj.eta = initVar;
            obj.workload = initVar;
        end
         
        function obj = pickMotorState(obj, idx)
            obj.time = obj.time(idx);
             obj.U = obj.U(idx);
             obj.I = obj.I(idx);
             obj.tau = obj.tau(idx);
             obj.Omega = obj.Omega(idx);
             obj.eta = obj.eta(idx);
             obj.workload = obj.workload(idx);
        end
         
        function objStore = storeMotorState(objStore, obj, idx)
             objStore.U(idx) = obj.U;
             objStore.I(idx) = obj.I;
             objStore.tau(idx) = obj.tau;
             objStore.Omega(idx) = obj.Omega;
             objStore.eta(idx) = obj.eta;
             objStore.workload(idx) = obj.workload;
        end
         
        function obj = incrementMotor( obj, dt )
            obj.idx = obj.idx + 1;
            obj.time(obj.idx) = obj.time(max(obj.idx-1,1)) + dt;
        end
         
         function obj = getMotorWorkload(obj)
            obj.workload(1) = obj.I / obj.params.I_max;
            obj.workload(obj.workload < 0) = 2;
         end
         
         function obj = setSpecificMotorParams(obj,Data)
            % call Motordata function
            [ K_V, I_0, R_i, m_Mot, S_max, I_max ] = ...
                motorLoadParams( Data, obj.params.name );
            % convert unit to 1/(V*s)
            K_V = K_V*2*pi/60;
            % set parameters
            obj.params.K_V = K_V;
            obj.params.I_0 = I_0;
            obj.params.R_i = R_i;
            obj.params.m = m_Mot;
            obj.params.S_max = S_max;
            obj.params.I_max = I_max;
         end
            
     end
end
