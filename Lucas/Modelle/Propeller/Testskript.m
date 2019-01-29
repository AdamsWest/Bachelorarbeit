% Propeller mit dem gleichen Durchmesser auslesen
clear
clc
% Diameter:

D = 9; %in [in]
Thrust = 70; % Schub in N
V_A = 35;
alpha = 0;

load('DATA_APC.mat'); 
DATA = evalin('base','DATA_APC');

[l,w] = size(DATA_APC);


for m = l:-1:1
    a = DATA{m,1};                                  % Speichern des Propellernames unter der Variablen a
    diameter = str2double(a(1:strfind(a,'x')-1));   % Extrahieren des Durchmessers aus dem Propellernamen
    if diameter ~= D                                % Entfernen der Zeile mit einem anderen Propellerdurchmesser
        DATA(m,:) = [];
    end
end