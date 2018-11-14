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

%% 
% Aufsuchen und finden aller Propellernamen, die den gleichen Durchmesser
% wie der gesuchte haben. Dabei werden alle anderen Propeller entfernt
m = 1;
while m <= l
    a = DATA{m,1};                                  % Speichern des Propellernames unter der Variablen a
    diameter = str2double(a(1:strfind(a,'x')-1));   % Extrahieren des Durchmessers aus dem Propellernamen
    if diameter ~= D                                % Entfernen der Zeile mit einem anderen Propellerdurchmesser
        DATA(m,:) = [];
        l = l - 1;
        m =  m - 1;
    end 
    m = m + 1;                                      % Erhöhen der Zählervariablen
end

%% 
% Entnahme jedes einzelnen Propellers
counter = 0;
ind = 1;
len = length(DATA);
d = 1;
while ind < len
    
    prop_name = DATA{ind,1};
    pitch = prop_name(strfind(prop_name,'x')+1:end);
    ind = find(strcmp(DATA(:,1),prop_name));
    
    
    ind_unten = min(ind);
    ind_oben = max(ind);
    ind = max(ind) + 1;
    
    
    % nicht passende Propeller werden entfernt (-E, 3- and 4-blade, -SF, etc.)
    if isnan(str2double(pitch)) == 1
        
        for d = ind_oben:-1:ind_unten
            DATA(d,:) = [];
            len = len - 1;
        end
        
    else
        
        [RPM, V, T, P, Tau] = Propeller_map(DATA,prop_name);
        
        
        counter = counter + 1;
        
        % Die Daten für RPM, V, T, usw. werden fortlaufenden Vektoren
        % zugewiesen
        assignin ('base',['prop_name' num2str(counter)], prop_name);
        assignin ('base',['pitch_' num2str(counter)], pitch);
        assignin ('base',['RPM_' num2str(counter)], RPM);
        assignin ('base',['T_' num2str(counter)], T);
        assignin ('base',['P_' num2str(counter)], P);
        assignin ('base',['Tau_' num2str(counter)], Tau);
    end
end

Omega = zeros(counter,1);
tau = zeros(counter,1);

%% Bestimmung des Pitches mit der geringsten Leistung

for i = 1:counter
    RPM = evalin('base',['RPM_' num2str(i)]);
    T = evalin('base',['T_' num2str(i)]);
    Tau = evalin('base',['Tau_' num2str(i)]);

    [Omega(i),tau(i)] = Propeller(V_A, alpha, Thrust, RPM, V, T, Tau);
end 
    
Tau_min = min(tau);
mintorqueind = find(tau == Tau_min);
mintorquepitch = evalin('base',['pitch_' num2str(mintorqueind)]);
    
x = 1:counter;
plot(x,tau,'rx')


