%% aerodynamic parameters
%
% empty mass (without propulsion components, the total mass will be 
% computed automatically), kg
aeroParams.m_empty = 0.354;
% payload, kg
aeroParams.payload = 0;

% type of aircraft
%   (0) : multicopter
%   (1) : airplane
aeroParams.type = 1;

% multicopter parameters (only used if type = 0)
%
% upper front face, m^2
aeroParams.A_copter = 0.15*0.05 + 0.12*0.02*4;
% sideways front face, m^2
aeroParams.A_copter_sideways = 1.5 * aeroParams.A_copter;
% drag coefficient (like a form factor), -
aeroParams.c_W_copter_upper = 1;
% sideways drag coefficient
aeroParams.c_W_copter_sideways = 1 * aeroParams.A_copter_sideways / ...
    aeroParams.A_copter;
% maximum lift coefficient (reached at +/- 45 degrees angle of attack), -
aeroParams.c_A_copter_max = 0.3;    
    

% airplane parameters (only used if type = 1)
% 
% design glide ratio, -
aeroParams.E_star = 4;
% design airspeed, m/s
aeroParams.V_star = 100/3.6;
% design air density, kg/m^3
aeroParams.rho_star = 1.225;


%% propeller parameters
%
% name of the propeller (diameter, inch x pitch, inch)
propParams.name = '9x7';
% number of propellers, -
propParams.n = 1;


%% motor parameters according to
% Drela, M. (2007). First-order DC electric motor model. Massachusetts
% Institute of Techno* sin(mission.gamma)logy.
%
% define if custom parameters or a specific motor should be used
%   (0) : use specific motor (the following lines will not be considered)
%   (1) : use custom parameter (to be defined in the following lines)
motorParams.isCustom = 0;
% define the name of a specific AXI motor within the data base
motorParams.name = 'AXI 2814/12 GOLD LINE';
% velocity constant, rpm/V
motorParams.K_V = 1390;
% zero-torque current, A
motorParams.I_0 = 1.8;
% resistance, Ohm
motorParams.R_i = 0.053;
% mass, kg
motorParams.m = 0.106;
% maximum allowed number of LiPo cells, -
motorParams.S_max = 3;
% maximum allowed current, A
motorParams.I_max = 35;


%% battery parameters
%
% energy density of the battery, J/kg
batParams.E_density = 890540;
% battery mass, kg
batParams.m = 0.56;
% battery capacity, A*s
batParams.C = 1;
% nominal voltage per cell, V
batParams.U_cell = 3.9;
% number of battery cells in a row, -
batParams.n_cell = 6;
% minimal allowed voltage per cell, V
batParams.U_cell_min = 2.85;
% Peukert's constant (estimation), -
batParams.P_Peukert = 1.00;
% maximum allowed C-rate, 1/h
batParams.C_rate_max = 30;
% define if the battery voltage is computed dynamically (only for LiPos)
%   (0) : fixed battery voltage (nominal voltage is used)
%   (1) : voltage depends on the state of charge and discharge rate
batParams.isDynVolt = 1;
% define if the battery mass is known
%   (0) : the mass is computed from the energy density, voltage and
%         capacity
%   (1) : the capacity is computed from the energy density, voltage and
%         mass
%   (2) : neither mass nor capacity will be computed
batParams.massKnown = 1;
    

%% environment parameters
%
% earth's gravity, m/s^2
envParams.g = 9.81;
% altitude of the airfield (relative to mean sea level), m
envParams.H_0 = 0;
% air temperature at the airfield, K
envParams.T_0 = 288.15;
% air pressure at the airfield, Pa
envParams.p_0 = 101325;
% air density at the airfield, kg/m^3
envParams.rho_0 = 1.225;
% adiabatic index of the air, -
envParams.kappa = 1.4;
% sidewind, m/s
envParams.u_Wg = 10;

    
%% mission parameters
%
% set the list of configurations in which the aircraft should operate
%   loiter: flight with maximum endurance
%   climb:  flight with most efficient climb
%   cruise: flight with highest range
missionParams.configuration = { 'climb', ...
                                'climb', ...
                                'climb' }; 
% set the required parameter for each mission configuration
%   loiter - parameter: time, s
%   climb - parameter:  vertical distance, m
%   cruise - parameter: longitudinal distance, m
missionParams.paramForConfig = [ 1000, ...
                                    15000, ...
                                    10000 ];
% set whether the computation of the mission should be dynamic or not
%   (0) : static computation of each mission element
%   (1) : dynamic quasi-static computation of the whole mission
missionParams.isDynamic = 1;
% set the time step (only used if isDynamic = 1), s
missionParams.dt = 10;
% set the maximum flight time (only used if isDynamic = 1), s
missionParams.maxTime = 24*3600;
% set how many times parameters will be changed during the optimization
missionParams.gamma = { 0:pi/2/90:pi/2; ...
                        0:pi/2/90:pi/2; ...
                        0:1:0 };
missionParams.V_K = { 0:5:0; ...
                        50:1:50; ...
                        50:1:50 };
                        
