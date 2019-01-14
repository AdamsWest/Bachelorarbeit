clc
clear  
close all

%% Dateinamen eingeben %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dateiname = 'Vergleich (Kennfeld)';

load('DATA_APC.mat');
load('Elektromodellflug');


%% Flugsystem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Handelt es sich bei dem zu untersuchenden Flugobjekt um einen Multicopter
% (1) oder um ein Flächenflugzeug (0)?

Abfrage_Flugsystem = 0;


%% Initialisierung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figures definieren
figure_C_Rest_V = figure;
figure_omega = figure;
figure_I_mot = figure;
figure_U_mot = figure;
figure_I_Bat = figure;
figure_PWM = figure;
figure_eta = figure;
figure_gamma = figure;


%% allgemeine Parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Motor
R_i = 0.123;            % Innenwiderstand in Ohm
K_V = 1400*2*pi/60;     % K_V Wert in 1/(V*s)
I_0 = 0.52;             % Leerlaufstrom in Ampere
I_max = 17;             % Max Continuous Current
m_Mot = 0.0365;         % Motorgewicht in kg

% Propeller
prop_name = '7x3.8';    % Propellerbezeichnung
n_Prop = 4;             % Anzahl der Propeller
%D = 14;                % Propellerdurchmesser in inch
%P_75 = 8;              % Propellersteigung bei 75% des Radius in inch
c_d0 = 0.05;            % Schaetzung des mittleren Nullwiderstandbeiwerts
a_alpha = 5;            % Anstieg des Auftriebsbeiwerts ueber dem Anstellwinkel (Profil), Schaetzung
alpha_stall = 10;       % Anstellwinkel, bei dem die Strömung abreisst in Grad, Schaetzung

% Batterie
E_Dichte = 750000;      % Energiedichte des LiPos in J/kg
N_Bat_cell = 4;         % Anzahl der Batteriezellen in Reihe
N_Bat_cell_p = 3;       % Anzahl der Batteriezellen parallel
C_Bat_cell = 3.160;     % Kapazität einer Zelle in Ah
U_Bat_cell = 3.9;       % nominale Spannung pro Batteriezelle
U_Bat_cell_min = 2.85;  % minimale Spannung pro Batteriezelle
P_Bat_Peukert = 1.05;   % Peukert-Konstante (Schaetzung)    
C_Rate_max = 30;        % maximale C-Rate bezogen auf eine nominale Entladezeit von 1 Stunde
m_Bat = 0.56;           % Batteriemasse in kg

% Missionsparameter
m_nutz = 0.25;          % Nutzlast in kg           


%% Parameter Multicopter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gesamtsystem
m_copter = 0.354;                       % Multicopter Leermasse in kg
A_copter = 0.15*0.05 + 0.12*0.02*4;     % obere Stirnflaeche des Multicopter in m^2
A_copter_seitlich = 1.5 * A_copter;     % seitliche Stirnflaeche des Multicopter in m^2
c_W_copter_oben = 1;                    % Widerstandsbeiwert des Multicopters 
c_W_copter_seitlich = 1 * A_copter_seitlich / A_copter;         % seitlicher Widerstandsbeiwert  des Multicopters
c_A_copter_max = 0.3;                   % maximaler Auftriebsbeiwert des Multicopters (bei +/-45° Anstellwinkel)


%% Parameter Flächenflugzeug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_flugzeug = 0.354;             % Flächenflugzeug Leermasse in kg
E = 3;                 		% Gleitzahl
E_stern = 4;			% Auslegungsgleitzahl
V_stern = 100/3.6;		% Auslegungsgeschwindigkeit in m/s
rho_stern = 1.225;		% Auslegungshöhe repräsentiert durch die Dichte (Bodennähe) in kg/m^3

% Diskretisierung des Bahnneigungswinkels zur Ermittlung des optimalen Steigwinkels 

gamma_min = 5;			% kleinster Bahnneigungswinkel
gamma_Delta = 5;		% Schrittweite Batteriemasse
gamma_max = 90;			% größter Bahnneigungswinkel


%% Flugparameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bahngeschwindigkeit
V_Kg = 10;                                  % Steiggeschwindigkeitin m/s
gamma_copter = 90 * pi/180;                 % Bahnanstellwinkel für den Multicopter    
 

%% Umgebungsparameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 9.81;                                   % Erdbeschleunigung in m/s^2

H_0 = 0;                                    % Höhe des Abflugplatzes über Normalnull in m
Delta_H = 50;                               % Inkrementweite in m 
H_max = 20000;                              % Maximalhöhe in m

T_0 = 288.15;                               % Temperatur in K am Flugplatz
p_0 = 101325;                               % Druck am Abflugplatz in Pa
rho_0 = 1.225;                              % Dichte am Startort in kg/m^3
kappa = 1.4;                                % Adiabatenexponent

u_Wg = 10;                                  % Seitenwindgeschwindigkeit in m/s


%% Aufruf des Hauptskripts: Leistungsberechnung starten %%%%%%%%%%%%%%%%%%%

run('Leistungsberechnung');