function [I_mot,Omega,tau] = Motor_umgk(PWM, U_mot, P_map,I_0,K_V,R_i)

% MOTOR berechnet die Drehzahl des Motors und den Motorstrom


%          [I_mot,Omega] = Motor_umgk(PWM, U_mot, P_map,I_0,K_V) bekommt
%          die Motorparameter übergeben und berechnet mittels der
%          Motorspannung den Motorstrom und und die Drehzahl


% Annahme PWM * P_max = U_mot * I_mot

I_mot = (PWM * max(max(P_map)))/U_mot;

tau = (I_mot - I_0)/K_V;

Omega = (U_mot - R_i*I_mot)*K_V;

end

