classdef BatteryClass
     properties
        time;
        idx = 1;
        U;
        C_rate;
        SoC;
        workload;
        params;
     end
     methods
         
        function obj = initBatteryProperties(obj, n)
            initVar = zeros(n,1);
            obj.time = initVar;
            obj.U = obj.params.U_cell * obj.params.n_cell * ones(n,1);
            obj.C_rate = initVar;
            obj.SoC = ones(n,1);
            obj.workload = [initVar, initVar];
        end
         
        function obj = pickBatteryState(obj, idx)
            obj.time = obj.time(idx);
            obj.U = obj.U(idx);
            obj.C_rate = obj.C_rate(idx);
            obj.SoC = obj.SoC(idx);
            obj.workload = obj.workload(idx);
        end
         
        function objStore = storeBatteryState(objStore, obj, idx)
            objStore.time(idx) = obj.time;
            objStore.U(idx) = obj.U;
            objStore.C_rate(idx) = obj.C_rate;
            objStore.SoC(idx) = obj.SoC;
            objStore.workload(idx,:) = obj.workload;
        end
         
        function obj = incrementBattery( obj, dt )
            obj.idx = obj.idx + 1;
            obj.time(obj.idx) = obj.time(max(obj.idx-1,1)) + dt;
        end
         
        function obj = batteryVoltage2(obj,I_bat)
            i_int = 1 - obj.SoC;
            Eo = obj.params.dynVolt.Eo;
            R = obj.params.dynVolt.R;
            K = obj.params.dynVolt.K;
            Q = obj.params.dynVolt.Q;
            A = obj.params.dynVolt.A;
            B = obj.params.dynVolt.B;
            n_cell = obj.params.n_cell;
            U_cell = Eo - R*I_bat - K * Q / (Q - i_int) * ...
                (i_int + I_bat*0) + A * exp(-B*i_int);
            obj.U = n_cell * U_cell;
        end
        
        function obj = computeBatteryMassOrCapacity(obj)
            obj.params.U_nom = obj.params.U_cell * obj.params.n_cell;
            if obj.params.massKnown == 0
                obj.params.m = obj.params.C * obj.params.U_nom / obj.params.E_density;
            elseif obj.params.massKnown == 1
                obj.params.C = obj.params.E_density * obj.params.m / obj.params.U_nom;
            end
        end
        
        function obj = integrateCapacity(obj)
            dt = obj.time( obj.idx ) - obj.time( max( obj.idx - 1, 1 ) );
            C = obj.params.C;
            P = obj.params.P_Peukert;
%             obj.C_rate = I_bat / C *3600;
            C_Peukert = C * (1/obj.C_rate(obj.idx-1))^(P-1);
            Delta_C = obj.C_rate(obj.idx-1)/3600 * C* dt * C / C_Peukert;
            obj.SoC(obj.idx) = obj.SoC( max( obj.idx - 1, 1 ) ) - ...
                Delta_C / C;
            if obj.SoC(obj.idx) <= 0
                disp('battery empty')
            end
        end
        
        function obj = getBatteryWorkload(obj)
            obj.workload(1) = obj.C_rate / obj.params.C_rate_max;
            obj.workload(2) = obj.params.U_cell_min / ...
                ( obj.U / obj.params.n_cell );
            obj.workload(obj.workload < 0) = 2;
        end
        
     end
 end
