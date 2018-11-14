function [PWM,eta_PWM] = ESC(U_mot,U_bat_nom)

PWM = U_mot / U_bat_nom;
            
if PWM > 0 && PWM < 0.5 
    eta_PWM = 0.7 * PWM + 0.50;
elseif PWM >= 0.5 && PWM <= 1
	eta_PWM = 0.2 * PWM + 0.75;
else
	eta_PWM = 1;
end

end