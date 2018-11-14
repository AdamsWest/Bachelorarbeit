function [I_bat,C_rate,Delta_C_bat,C_Rest_V] = Batterie(PWM,eta_PWM,I_mot,n_Prop,C_bat,P_bat,Delta_C_bat,t_Flug)

I_bat = PWM * I_mot / eta_PWM * n_Prop;
C_rate = I_bat / (C_bat/3600);                              % C-Rate bezogen auf eine nominale Entladezeit von 1 Stunde
C_bat_Peukert = C_bat * (1/C_rate)^(P_bat-1);               % nutzbare Kapazitaet nach Peukert
Delta_C_bat = I_bat * t_Flug + Delta_C_bat;                 % entnommene Ladung
C_Rest_V = (C_bat_Peukert - Delta_C_bat) / C_bat_Peukert;   % Ladezustand als Verhaltnis

end