clc
clear  
close all

%% Dateinamen eingeben %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dateiname = 'Flächenflugzeug';                        % Hier Dateiname
% festlegn oder unten definieren lassen

load('DATA_APC.mat');
load('Elektromodellflug');
load('axi_motor_db.mat');


%% Flugsystem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Handelt es sich bei dem zu untersuchenden Flugobjekt um einen Multicopter
% (1) oder um ein Flächenflugzeug (0)?

% Abfrage_Flugsystem = 1;


%% Initialisierung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figures definieren
figure_ges = figure;


%% allgemeine Parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Motor
motor_name = axi_motor_db{21,1}; % Motorname
[K_V, I_0, R_i, m_Mot, S_max, I_max] = Motordata('axi_motor_db',motor_name);
K_V = K_V*2*pi/60;          % Umrechnung in 1/(V*s)
% R_i = 0.123;            % Innenwiderstand in Ohm
% K_V = 1400*2*pi/60;     % K_V Wert in 1/(V*s)
% I_0 = 0.56;             % Leerlaufstrom in Ampere
% I_max = 17;             % Max Continuous Current
% m_Mot = 0.0365;         % Motorgewicht in kg

% Propeller
prop_name = '11x4';    % Propellerbezeichnung
n_Prop = 4;             % Anzahl der Propeller
D = 10;                % Propellerdurchmesser in inch
%P_75 = 8;              % Propellersteigung bei 75% des Radius in inch
c_d0 = 0.05;            % Schaetzung des mittleren Nullwiderstandbeiwerts
a_alpha = 5;            % Anstieg des Auftriebsbeiwerts ueber dem Anstellwinkel (Profil), Schaetzung
alpha_stall = 10;       % Anstellwinkel, bei dem die Strömung abreisst in Grad, Schaetzung

% Batterie
E_Dichte = 890540;      % Energiedichte des LiPos in J/kg
N_Bat_cell = 5;         % Anzahl der Batteriezellen in Reihe
N_Bat_cell_p = 3;       % Anzahl der Batteriezellen parallel
C_Bat_cell = 3.12;     % Kapazität einer Zelle in Ah
U_Bat_cell = 4;       % nominale Spannung pro Batteriezelle
U_Bat_cell_min = 2.85;  % minimale Spannung pro Batteriezelle
P_Bat_Peukert = 1.00;   % Peukert-Konstante (Schaetzung)    
C_Rate_max = 30;        % maximale C-Rate bezogen auf eine nominale Entladezeit von 1 Stunde
% m_Bat = 0.56;           % Batteriemasse in kg

% Batteriemassendiskreitisierung
% m_Bat_min = 0.5;
% m_Bat_Delta = 0.02;
% m_Bat_max = 0.72;

% Missionsparameter
m_nutz = 0.0;          % Nutzlast in kg           

m_ges = (n_Prop * m_Mot)/(4*0.0365/1.06);     % Gesamtmasse
% m_Bat = m_ges * (0.56/1.06);              % Batteriemasse
m_copter = m_ges * (0.354/1.06);           % Leermasse
m_Bat = 2/3*m_ges;
% m_Bat = (m_Mot*n_Prop+m_copter)*2;
m_Bat = m_Bat - 0.25;


%% Parameter Multicopter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gesamtsystem
% m_copter = 0.354;                       % Multicopter Leermasse in kg
% A_copter = 0.15*0.05 + 0.12*0.02*4;     % obere Stirnflaeche des Multicopter in m^2
A_copter = 0.15*0.05 + (D/2*0.0254)*1.2*0.02* n_Prop;     % obere Stirnflaeche des Multicopter in m^2
A_copter_seitlich = 1.5 * A_copter;     % seitliche Stirnflaeche des Multicopter in m^2
c_W_copter_oben = 0.5;                    % Widerstandsbeiwert des Multicopters 
c_W_copter_seitlich = 1 * A_copter_seitlich / A_copter;         % seitlicher Widerstandsbeiwert  des Multicopters
c_A_copter_max = 0.3;                   % maximaler Auftriebsbeiwert des Multicopters (bei +/-45° Anstellwinkel)

% Diskretisierung der Steiggeschwindigkeit
V_Kg_min = 1;			% kleinster Bahnneigungswinkel
V_Kg_Delta = 1;		% Schrittweite Batteriemasse
V_Kg_max = 50;			% größter Bahnneigungswinkel

%% Flugparameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bahngeschwindigkeit
% V_Kg = 10;                                  % Steiggeschwindigkeitin m/s
gamma_copter = 90 * pi/180;                 % Bahnanstellwinkel für den Multicopter    
 

%% Umgebungsparameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 9.81;                                   % Erdbeschleunigung in m/s^2

H_0 = 0;                                    % Höhe des Abflugplatzes über Normalnull in m
Delta_H = 100;                              % Inkrementweite in m 
H_max = 18500;                              % Maximalhöhe in m

T_0 = 288.15;                               % Temperatur in K am Flugplatz
p_0 = 101325;                               % Druck am Abflugplatz in Pa
rho_0 = 1.225;                              % Dichte am Startort in kg/m^3
kappa = 1.4;                                % Adiabatenexponent

u_Wg = 100/3.6;                                  % Seitenwindgeschwindigkeit in m/s


%% Festlegung des Dateinamen
    
% Dateiname = ['Multicopter, m_Mot = ' num2str(m_Mot) ', n_Prop = ' num2str(n_Prop) ', K_V = ' num2str(K_V*60/(2*pi)) ', Prop = ' prop_name ...
%     ', n_Bat_cell = ' num2str(N_Bat_cell) ', c_W = ' num2str(c_W_copter_oben) ', u_Wg = ' num2str(u_Wg) 'ms.pdf'];

Dateiname = 'aeromet_randbedingungen';
%% Aufruf des Hauptskripts: Leistungsberechnung starten %%%%%%%%%%%%%%%%%%%

run('Leistungsberechnung'); % _Var_m_Bat