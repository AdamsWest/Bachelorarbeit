clc
clear  
close all

%% Dateinamen eingeben %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dateiname = 'Vergleich (Kennfeld)';

load('DATA_APC.mat');


%% Flugsystem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Handelt es sich bei dem zu untersuchenden Flugobjekt um einen Multicopter
% (1) oder um ein Flächenflugzeug (0)?

Abfrage_Flugsystem = 1;


%% Initialisierung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figures definieren
figure_C_Rest_V = figure;
figure_omega = figure;
figure_I_mot = figure;
figure_U_mot = figure;
figure_I_Bat = figure;
figure_PWM = figure;
figure_Pitch = figure;

%% allgemeine Parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Motor
R_i = 0.123;            % Innenwiderstand in Ohm
K_V = 1400*2*pi/60;     % K_V Wert in 1/(V*s)
I_0 = 0.52;             % Leerlaufstrom in Ampere
I_max = 40;             % Max Current
m_Mot = 0.0365;         % Motorgewicht in kg

% Propeller
%prop_name = '9x4';    % Propellerbezeichnung
n_Prop = 4;             % Anzahl der Propeller
D = 9;                  % Propellerdurchmesser in inch
%P_75 = 8;             % Propellersteigung bei 75% des Radius in inch
c_d0 = 0.05;            % Schaetzung des mittleren Nullwiderstandbeiwerts
a = 5;                  % Anstieg des Auftriebsbeiwerts ueber dem Anstellwinkel (Profil), Schaetzung
alpha_stall = 10;       % Anstellwinkel, bei dem die Strömung abreisst in Grad, Schaetzung

% Batterie
E_Dichte = 750000;      % Energiedichte des LiPos in J/kg
N_Bat_cell = 4;         % Anzahl der Batteriezellen
U_Bat_cell = 3.7;       % nominale Spannung pro Batteriezelle
U_Bat_cell_min = 3.4;   % minimale Spannung pro Batteriezelle
P_Bat = 1.05;           % Peukert-Konstante (Schaetzung)    
C_Rate_max = 50;        % maximale C-Rate bezogen auf eine nominale Entladezeit von 1 Stunde
m_Bat = 1;           % Batteriemasse in kg


%% Parameter Multicopter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gesamtsystem
m_copter = 0.5;                        % Quadrocopter Leermasse in kg
A_copter = 0.1*0.1  + 0.1*0.02*4;     % obere Stirnflaeche des Quadrocopters in m^2
%A_copter = 2 * A_copter;
A_copter_seitlich = 1.5 * A_copter;     % seitliche Stirnflaeche des Quadrocopters in m^2
c_W_copter_oben = 1;                    % Widerstandsbeiwert des Quadrocopters 
c_W_copter_seitlich = 1 * A_copter_seitlich / A_copter;         % seitlicher Widerstandsbeiwert  des Quadrocopters
c_A_copter_max = 0.3;                   % maximaler Auftriebsbeiwert des Quadrocopters (bei +/-45° Anstellwinkel)


%% Parameter Flächenflugzeug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon = 1/38;          % reziproke Gleitzahl
m_flugzeug = 15;          % Flächenflugzeug Leermasse in kg


%% Flugparameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bahngeschwindigkeit
V_Kg = 10;                                  % Steiggeschwindigkeitin m/s
%V_Profil = [11 13 14.5 12 10 7 2 3];        % Geschwindigkeitsprofil für den Steigflug (Russland)
gamma = 90 * pi/180;                        % Bahnanstellwinkel für den Multicopter    


%% Umgebungsparameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 9.81;                                   % Erdbeschleunigung in m/s^2

H_0 = 0;                                    % Höhe des Abflugplatzes über Normalnull in m
Delta_H = 100;                              % Inkrementweite in m 
H_max = 20000;                              % Maximalhöhe in m

Temp_0 = 288.15;                               % Temperatur in K am Flugplatz

rho_0 = 1.225;                              % Dichte am Startort in kg/m^3

u_Wg = 10;                                  % Seitenwindgeschwindigkeit in m/s


%% Aufruf des Hauptskripts: Leistungsberechnung starten %%%%%%%%%%%%%%%%%%%

run('Leistungsberechnung_vprop');