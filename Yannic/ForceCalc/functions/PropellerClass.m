classdef PropellerClass
     properties
         time;
         idx = 1;
        Omega;
        V;
        T;
        P;
        tau;
        M_tip;
        eta;
        mapRho;
        workload;
        params;
     end
     methods
         
         function obj = initPropellerProperties(obj, n)
             initVar = zeros(n,1);
             obj.time = initVar;
             obj.Omega = initVar;
             obj.V = initVar;
             obj.T = initVar;
             obj.P = initVar;
             obj.tau = initVar;
             obj.M_tip = initVar;
             obj.eta = initVar;
             obj.workload = [initVar, initVar];
         end
         
        function obj = pickPropellerState(obj, idx)
            obj.time = obj.time(idx);
            obj.Omega = obj.Omega(idx);
            obj.V = obj.V(idx);
            obj.T = obj.T(idx);
            obj.P = obj.P(idx);
            obj.tau = obj.tau(idx);
            obj.M_tip = obj.M_tip(idx);
            obj.eta = obj.eta(idx);
            obj.workload = obj.workload(idx);
%             obj.mapRho = obj.mapRho;
%             obj.workload = obj.workload;
%             obj.params = obj.params;
        end
         
        function objStore = storePropellerState(objStore, obj, idx)
             objStore.Omega(idx) = obj.Omega;
             objStore.V(idx) = obj.V;
             objStore.T(idx) = obj.T;
             objStore.P(idx) = obj.P;
             objStore.tau(idx) = obj.tau;
             objStore.M_tip(idx) = obj.M_tip;
             objStore.eta(idx) = obj.eta;
             objStore.workload(idx,:) = obj.workload;
%              obj.mapRho = obj.mapRho;
%              obj.workload = obj.workload;
%              obj.params = obj.params;
        end
         
        function obj = incrementPropeller( obj, dt )
            obj.idx = obj.idx + 1;
            obj.time(obj.idx) = obj.time(max(obj.idx-1,1)) + dt;
        end
         
         function obj = setPropellerMapRho(obj,rho,rho_0)
             obj.mapRho.T = obj.params.map.T_0 * rho / rho_0;
             obj.mapRho.P = obj.params.map.P_0 * rho / rho_0;
             obj.mapRho.tau = obj.params.map.tau_0 * rho / rho_0;
         end
         
         function obj = computeTipMachNumber(obj,a)
             V_tip = obj.Omega * obj.params.r;
             obj.M_tip = V_tip / a;
         end     
         
         function obj = setPropGeometicParams(obj)
         	name = obj.params.name;
            obj.params.D = str2double(name(1:strfind(name,'x')-1));      % Durchmesser extrahieren
            P_75 = name(strfind(name,'x')+1:end);             % Pitch extrahieren
            while isnan(str2double(P_75)) == 1
                P_75(end) = [];
            end
            obj.params.P_75 = str2double(P_75);                    % Pitch festlegen
            obj.params.r = obj.params.D * 0.0254 / 2;                         % Propellerradius in Meter
            obj.params.F = pi * obj.params.r^2;                               % Fl�che eines Propellers in Quadratmeter
            obj.params.Theta_75 = atan( 4*obj.params.P_75 / (3*pi * obj.params.D) );     % geometrischer Anstellwinkel des Propellers bei 75% des Radius
         end
         
        function obj = getPropellerWorkload(obj)
            obj.workload(1) = obj.Omega / obj.params.Omega_max;
            obj.workload(2) = obj.M_tip;
            obj.workload(obj.workload < 0) = 2;
        end
         
     end
 end
