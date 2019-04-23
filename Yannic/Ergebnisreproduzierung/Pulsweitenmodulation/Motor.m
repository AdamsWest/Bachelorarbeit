function [U_mot,I_mot,eta_mot] = Motor(tau,K_V,I_0,R_i,Omega)  

% MOTOR berechnet den Motorstrom und -spannung sowie den Wirkungsgrad
%
%   Motor berechnet anhand eines einfachen Motormodells nach Drela (2007) 
%   den Motorstrom und -spannung
%
% Syntax:  [U_mot,I_mot] = Motor(tau,K_V,I_0,R_i,Omega)
%
% Inputs:
%    tau	Drehmoment des Propellers
%    K_V 	Motorkonstante
%    I_0	Leerlaufstrom
%    R_i	Innenwiderstand
%    Omega	Drehzahl des Propellers
%
% Outputs:
%    U_mot	Motorspannung
%    I_mot	Motorstrom
%    eta_mot    Motorwirkungsgrad
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    [ values, derivatives ] = Untitled( x, myStruct.y, 2 )
%    [ values, ~ ] = Untitled( [ 1:.1:100 ] , myStruct.y, 1 )
%
% See also: PROPELLER,  ESC,  BATTERIE

%   Copyright 2018 TU-Braunschweig
% ******************************************************************************


I_mot = tau*K_V + I_0;                      % Motorstrom
U_mot = Omega/(K_V) + R_i*I_mot;            % Motorspannung
eta_mot = (tau*Omega)/(U_mot*I_mot);        % Motorwirkungsgrad

end

