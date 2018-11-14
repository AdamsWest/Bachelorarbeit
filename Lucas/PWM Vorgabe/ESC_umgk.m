function [U_mot,eta_PWM] = ESC_umgk(PWM,U_bat_nom)

% ESC_UMGK bestimmt mittels der Thrust Position den Motorstrom

%          [U_mot,eta_PWM] = ESC_umgk(PWM,U_bat_nom) bekommt die Thrust
%          lever Postion gegeben und ermittelt mithilfe der
%          Batteriespannung die Motorspannung und den Wirkungsgrad des
%          Reglers.

U_mot = PWM * U_bat_nom;

if PWM > 0 && PWM < 0.5 
    eta_PWM = 0.7 * PWM + 0.50;
elseif PWM >= 0.5 && PWM <= 1
	eta_PWM = 0.2 * PWM + 0.75;
else
	eta_PWM = 1;
end

end

