clc
clear all

%% Dateinamen eingeben

Dateiname = 'Vergleich (ursprl.)';


%% Mission

% Bei jeder Spalte handelt es sich um eine Flugphase, es sind beliebig
% viele Flugphasen moeglich - sofern diese konsistent angegeben werden

% Bei jeder Flugphase kann eine andere Windgeschwindigkeit angegeben
% werden. Bei Hin- und Rueckflug kann beispielsweise das Vorzeichen
% vertauscht werden.

% Seitenwindgeschwindigkeit
u_Wg = [20;20;20];

% Handelt es sich bei der Flugphase um hovern (0) oder um eine
% Bahngeschwindigkeit (1) ?
Abfrage_hovern = [0;1;0]; 

% Je nach dem was bei hovern eingetragen wurde, werden im Folgenden nur die
% Eintraege unter Hovern bzw. unter Bahngeschwindigkeit beruecksichtigt

% Hovern
t_hover = [0;60;0];

% Bahngeschwindigkeit
Abfrage_V_variabel1 = [1;0;0];           % wenn bei Abfrage_V_variabel1 oder Abfrage_V_variabel2 eine 1 steht, ist der Eintrag bei V_Kg belanglos
Abfrage_V_variabel2 = [0;0;1];
Bahngeschwindigkeit = [0;0;0];
gamma = [90;0;-90] * pi/180;
Strecke = [1000,0,1000];




%% Initialisierung

% figures definieren
figure_C_Rest_V = figure;
figure_omega = figure;
figure_I_mot = figure;
figure_U_mot = figure;


% Parameter Multicopter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Motor
R_i = 0.057;            % Innenwiderstand in Ohm
K_V = 120*2*pi/60;
I_0 = 0.7;              % Leerlaufstrom in Ampere
I_max = 80;

% Propeller
n_Prop = 4;             % Anzahl der Propeller
D = 26;                 % Propellerdurchmesser in inch
P_75 = 8.5;             % Propellersteigung bei 75% des Radius in inch
c_d0 = 0.05;            % Schaetzung des mittleren Nullwiderstandbeiwerts
a = 5;                  % Anstieg des Auftriebsbeiwerts ueber dem Anstellwinkel (Profil), Schaetzung
alpha_stall = 10;       % Anstellwinkel, bei dem die Strömung abreisst in Grad, Schaetzung

% Antrieb (Testdaten Motor-Propeller-Kombination): verfuegbar auf rctigermotor.com
Throttle_static = [0;50;65;75;85;100];                      % Testdaten fuer Throttle (statisch)
Amps_static = [I_0;10.1;18.3;25;33.5;47.4];                 % Testdaten fuer Strom in Ampere (statisch)
Thrust_static = [0;4590;6700;8110;9640;12420]/1000*9.81;    % Testdaten fuer Schub in Newton (statisch)
omega_static = [0;2860;3600;3900;4300;4600]/60*2*pi;        % Testdaten fuer Winkelgeschw. in rad/s (statisch)

% Batterie
E_Dichte = 444000;      % Energiedichte des LiPos in J/kg
N_Bat_cell = 12;        % Anzahl der Batteriezellen
U_Bat_cell = 3.7;       % nominale Spannung pro Batteriezelle
U_Bat_cell_min = 3.4;   % minimale Spannung pro Batteriezelle
P_Bat = 1.05;           % Peukert-Konstante (Schaetzung)    
C_Rate_max = 50;        % maximale C-Rate bezogen auf eine nominale Entladezeit von 1 Stunde

% Gesamtsystem
m_copter = 15;                      % Quadrocopter Leermasse in kg
A_copter = 0.4*0.4 + 0.043*0.3*4;   % obere Stirnflaeche des Quadrocopters in m^2
A_copter_seitlich = 1.5 * A_copter; % seitliche Stirnflaeche des Quadrocopters in m^2
c_W_copter_oben = 1;                % Widerstandsbeiwert des Quadrocopters 
c_W_copter_seitlich = 1 * A_copter_seitlich / A_copter;         % seitlicher Widerstandsbeiwert  des Quadrocopters
c_A_copter_max = 0.3;               % maximaler Auftriebsbeiwert des Quadrocopters (bei +/-45° Anstellwinkel)


% Flugmission %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%H = 1000;               % Flughoehe in Meter
%t_Hover = 1*60;         % Hover-Zeit in Sekunden
%V_Sink = 2:2:6;         % Sinkgeschwindigkeit (positiv) in m/s

% Diskretisierung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Diskretisierung der Batteriemasse
m_Bat_min = 4;          % kleinste Batteriemasse
m_Bat_Delta = 1.5;       % Schrittweite Batteriemasse
m_Bat_max = 12;         % groesste Batteriemasse

% Diskretisierung der 1. zu optimierenden Bahngeschwindigkeit
V_Kg_1_min = 0;           % kleinste Geschwindigkeit
V_Kg_1_Delta = 0.25;          % Schrittweite Geschwindigkeit
V_Kg_1_max = 15;          % groesste Geschwindigkeit

% Diskretisierung der 2. zu optimierenden Bahngeschwindigkeit
V_Kg_2_min = 2;
V_Kg_2_max = 6;
V_Kg_2_Delta = (V_Kg_2_max - V_Kg_2_min) / 2;

% Umgebungsparameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho = 1.1;              % Luftdichte in kg/m^3
g = 9.81;               % Erdbeschleunigung in m/s^2



%% Aufruf des Hauptskripts: Leistungsberechnung starten

run('MulticopterLeistungsberechnung');