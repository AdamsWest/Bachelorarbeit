%% create propeller map with user defined function Propeller_map
load('DATA_APC.mat');
copter.prop_name = '10x5E';                             % propeller name within the APC data base
[prop_map.RPM, prop_map.V, prop_map.T, prop_map.P,prop_map.Tau] = Propeller_map(DATA_APC,copter.prop_name);
%for m = 1:length(prop_map.V)                                            % 46 is hard coded?
%    prop_map.tau(:,m) = prop_map.P(:,m)./(prop_map.RPM*2*pi/60);
%end