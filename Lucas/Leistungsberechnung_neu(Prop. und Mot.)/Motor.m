function [U_mot,I_mot] = Motor(tau,K_V,I_0,R_i,Omega_neu)

% MOTOR berechnet den neuen Motorstrom und -spannung nach einem einfachen
% Motormodell (Drela 2007)

%   [U_mot,I_mot] = Motor(tau,K_V,I_0,R_i,) berechnet anhand eines 
%   einfachen Motormodells nach Drela (2007) den Motorstrom und -spannung 

I_mot = tau*K_V + I_0;                       % Motorstrom
U_mot = Omega_neu/(K_V) + R_i*I_mot;         % Motorspannung

end

