


% Drehzahl und Drehmoment bestimmen

%[Omega(i),tau(i)] = Propeller(V_A, alpha(i), Thrust(i), RPM_map, V_map, T_map, TAU_map);

omega = zeros(counter,1);
TAU = zeros(counter,1);

%% Bestimmung des Pitches mit der geringsten Leistung
% Theoretisch in gleicher While-Schleife möglich wie in
pitch_all = zeros(counter,1);
for n = 1:counter
    RPM = evalin('base',['RPM_' num2str(n)]);
    T = evalin('base',['T_' num2str(n)]) * rho_2/rho_1;
    Tau = evalin('base',['Tau_' num2str(n)]) * rho_2/rho_1;
    
    a = str2double(evalin('base',['pitch_' num2str(n)]));
    pitch_all(n) = a;
    
    
    if Thrust(i) > max(max(T))            % !!!!!!!!!!!!! Abändern                      % wenn Schub zu gross (Ergebnis verwerfen)
        Omega(i) = NaN;
        I_mot(i) = NaN;
        C_Rate(i) = NaN;
        C_Rest_V(i) = NaN;
        I_bat_ges(i) = NaN;
        PWM(i) = NaN;
        Pitch(i) = NaN;
    else
        
        [omega(n),TAU(n)] = Propeller(V_A, alpha(n), Thrust(n), RPM, V, T, Tau);
        
        
        % Motorzustand berechnen
        
        [U_mot(i),I_mot(i)] = Motor(tau(i),K_V,I_0,R_i,Omega(i));
        
        
        % Zustand der Motorregler berechnen
        
        [PWM(i),eta_PWM] = ESC(U_mot(i),U_Bat_nom);
        
        
        % Batteriezustand berechnen
        
        [I_Bat,C_Rate(i),Delta_C_Bat,C_Rest_V(i)] = Batterie(PWM(i),eta_PWM,I_mot(i),n_Prop,C_Bat,P_Bat,Delta_C_Bat,t_Flug);
        I_bat_ges(i) = I_Bat;
        
        
    end
    
    
    
    
end

tau(i) = min(TAU);
mintorqueind = find(TAU == tau(i));
Omega(i) = omega(mintorqueind);
Pitch(i) = pitch_all(mintorqueind);
%mintorquepitch = evalin('base',['pitch_' num2str(mintorqueind)]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialisierungen
omega = zeros(counter,1);
TAU = zeros(counter,1);
pitch_all = zeros(counter,1);

for n = 1:counter
    
    RPM = evalin('base',['RPM_' num2str(n)]);
    T = evalin('base',['T_' num2str(n)]) * rho_2/rho_1;
    Tau = evalin('base',['Tau_' num2str(n)]) * rho_2/rho_1;
    
    a = str2double(evalin('base',['pitch_' num2str(n)]));
    pitch_all(n) = a;
    
    if Thrust(i) > max(max(T))
        omega(n) = NaN;
        TAU(n) = NaN;
        
    else
        
        [omega(n),TAU(n)] = Propeller(V_A, alpha(n), Thrust(n), RPM, V, T, Tau);
        
    end
    
end

tau(i) = min(TAU);
mintorqueind = find(TAU == tau(i));
Omega(i) = omega(mintorqueind);
Pitch(i) = pitch_all(mintorqueind);

if isnan(mean(omega)) == 1
    
    Omega(i) = NaN;
    I_mot(i) = NaN;
    C_Rate(i) = NaN;
    C_Rest_V(i) = NaN;
    I_bat_ges(i) = NaN;
    PWM(i) = NaN;
    Pitch(i) = NaN;
    
else
    
    tau(i) = min(TAU);
mintorqueind = find(TAU == tau(i));
Omega(i) = omega(mintorqueind);
Pitch(i) = pitch_all(mintorqueind);

end
    
    
