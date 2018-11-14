clc
clear 
clf 
close all

%% Dateinamen eingeben

Dateiname = 'Vergleich (Kennfeld)';

load('DATA_APC.mat');


%% Mission

% Bei jeder Spalte handelt es sich um eine Flugphase, es sind beliebig
% viele Flugphasen moeglich - sofern diese konsistent angegeben werden

% Bei jeder Flugphase kann eine andere Windgeschwindigkeit angegeben
% werden. Bei Hin- und Rueckflug kann beispielsweise das Vorzeichen
% vertauscht werden.

% Handelt es sich bei dem zu untersuchenden Flugobjekt um einen Multicopter
% (1) oder um ein Flächenflugzeug (0)?

Abfrage_Flugsystem = 0;

% Windgeschwindigkeiten in x-Richtung in m/s, 
% V_Wg =(u_Wg, v_Wg, w_Wg)^-1
% v_Wg und w_Wg in dieser Betrachtung irrelevant
% Falls kein Wind vorhanden 0 eintragen
u_Wg = 10;

% Handelt es sich bei der Flugphase um hovern (0) oder um eine
% Bahngeschwindigkeit (1) ?
%Abfrage_hovern = [0;1;0]; 

% Je nach dem was bei hovern eingetragen wurde, werden im Folgenden nur die
% Eintraege unter Hovern bzw. unter Bahngeschwindigkeit beruecksichtigt

% Hovern
%t_hover = [0;0;0];

% Bahngeschwindigkeit
Abfrage_V_variabel1 = [1;0;0];              % wenn bei Abfrage_V_variabel1 oder Abfrage_V_variabel2 eine 1 steht, ist der Eintrag bei V_Kg belanglos
Abfrage_V_variabel2 = [0;0;1];
Bahngeschwindigkeit = [0;0;0];  
gamma = 90 * pi/180;                        % Bahnanstellwinkel für den Multicopter           
Strecke = [1000,0,1000];




%% Initialisierung

% figures definieren
figure_C_Rest_V = figure;
figure_omega = figure;
figure_I_mot = figure;
figure_U_mot = figure;
figure_I_Bat = figure;
figure_PWM = figure;

%% allgemeine Parameter

% Motor
R_i = 0.123;            % Innenwiderstand in Ohm
K_V = 1400*2*pi/60;
I_0 = 0.52;             % Leerlaufstrom in Ampere
I_max = 40;             % Max Current

% Propeller
prop_name = '7x3.8';    % Propellerbezeichnung
n_Prop = 4;             % Anzahl der Propeller
D = 7;                  % Propellerdurchmesser in inch
P_75 = 3.8;             % Propellersteigung bei 75% des Radius in inch
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
m_Bat = 0.55;           % Batteriemasse in kg

%% Parameter Multicopter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gesamtsystem
m_copter = 1.06;                        % Quadrocopter Leermasse in kg
A_copter = 0.1*0.1  + 0.1*0.02*4;     % obere Stirnflaeche des Quadrocopters in m^2
A_copter = 2 * A_copter;
A_copter_seitlich = 1.5 * A_copter;     % seitliche Stirnflaeche des Quadrocopters in m^2
c_W_copter_oben = 1;                    % Widerstandsbeiwert des Quadrocopters 
c_W_copter_seitlich = 1 * A_copter_seitlich / A_copter;         % seitlicher Widerstandsbeiwert  des Quadrocopters
c_A_copter_max = 0.3;                   % maximaler Auftriebsbeiwert des Quadrocopters (bei +/-45° Anstellwinkel)


%% Parameter Flächenflugzeug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon = 1/38;          % reziproke Gleitzahl
m_flugzeug = 15;          % Flächenflugzeug Leermasse in kg


% Flugmission %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%H = 1000;               % Flughoehe in Meter
%t_Hover = 1*60;         % Hover-Zeit in Sekunden
%V_Sink = 2:2:6;         % Sinkgeschwindigkeit (positiv) in m/s

% Diskretisierung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Diskretisierung der Batteriemasse
%m_Bat_min = 4;          % kleinste Batteriemasse
%m_Bat_Delta = 1.5;       % Schrittweite Batteriemasse
%m_Bat_max = 12;         % groesste Batteriemasse

% Diskretisierung der 1. zu optimierenden Bahngeschwindigkeit
%V_Kg_1_min = 0;           % kleinste Geschwindigkeit
%V_Kg_1_Delta = 0.25;          % Schrittweite Geschwindigkeit
%V_Kg_1_max = 15;          % groesste Geschwindigkeit

% Diskretisierung der 2. zu optimierenden Bahngeschwindigkeit
%V_Kg = 10;                                  % Steiggeschwindigkeitin m/s
V_Profil = [11 13 14.5 12 10 7 2 3];        % Geschwindigkeitsprofil für den Steigflug (Russland)
%V_Kg_2_min = 2;
%V_Kg_2_max = 6;
%V_Kg_2_Delta = (V_Kg_2_max - V_Kg_2_min) / 2;


%% Umgebungsparameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 9.81;                                   % Erdbeschleunigung in m/s^2

H_0 = 0;                                    % Höhe des Abflugplatzes über Normalnull in m
Delta_H = 100;                              % Inkrementweite in m 
H_max = 15000;                              % Maximalhöhe in m

T_0 = 288.15;                               % Temperatur in K am Flugplatz

rho_0 = 1.225;                              % Dichte am Startort in kg/m^3



%% Aufruf des Hauptskripts: Leistungsberechnung starten

run('Leistungsberechnung');