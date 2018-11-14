function [U_mot,I_mot] = Motor(tau,K_V,I_0,R_i,Omega)

% MOTOR berechnet den Motorstrom und -spannung 

%   [U_mot,I_mot] = Motor(tau,K_V,I_0,R_i,) berechnet anhand eines 
%   einfachen Motormodells nach Drela (2007) den Motorstrom und -spannung 

I_mot = tau*K_V + I_0;                      % Motorstrom
U_mot = Omega/(K_V) + R_i*I_mot;            % Motorspannung

end

