classdef EscClass
     properties
         time;
         idx = 1;
        PWM;
        eta;
        I_bat;
        params;
     end
     methods
         
          function obj = initEscProperties(obj, n)
             initVar = zeros(n,1);
             obj.time = initVar;
             obj.PWM = initVar;
             obj.eta = initVar;
             obj.I_bat = initVar;
          end
         
          function obj = pickEscState(obj, idx)
              obj.time = obj.time(idx);
             obj.PWM = obj.PWM(idx);
             obj.eta = obj.eta(idx);
             obj.I_bat = obj.I_bat(idx);
          end
          
          function objStore = storeEscState(objStore, obj, idx)
              objStore.time(idx) = obj.time;
             objStore.PWM(idx) = obj.PWM;
             objStore.eta(idx) = obj.eta;
             objStore.I_bat(idx) = obj.I_bat;
          end
         
         function obj = incrementEsc( obj, dt )
            obj.idx = obj.idx + 1;
            obj.time(obj.idx) = obj.time(max(obj.idx-1,1)) + dt;
        end
                   
     end
end
