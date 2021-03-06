function [Thrust] = Vertikalflug(m,g,E_stern,V_stern,rho_stern,V_A,rho)
% VERTIKALFLUG berechnet den Schub f�r ein l�chenflugzeug im vertikalen
% Steiglfug
%
%   Vertikalflug(m,g,E_stern,V_stern,rho_stern,V_A,rho) errechnet f�r aus
%   einem gegebenen Auslegungszustand und f�r eine gegebene
%   Steiggeschwindigkeit den ben�tigten Schub.
%
%
%   Dabei bleiben Windeinfl�sse unber�cksichtig, sodass gilt
%   alpha=sigma=beta= 0. Der Wind hat keinen Einfluss auf das Steigverm�gen
%   eines Flugzeugs sondern nur auf die zur�ckgelegte Strecke �ber Grund
%   und auf den Steigwinkel, der dennoch ohne Belang in dieser Betrachtung
%   ist.
%   (vgl. Scheiderer, J., �/Angewandte Flugleistung: Eine Einf�hrung in die
%   operationelle Flugleistung vom Start bis zur Landung/�, Springer, 2008,
%   S. 241f.)
%
% Syntax:  [Thrust] = Vertikalflug(m,g,E_stern,V_stern,rho_stern,V_A,rho)
%
%
% Inputs:
%    m		Masse des Fl�chenflugzeugs
%    g		Erdbeschleunigung
%    E_stern	Auslegungsgleitzahl
%    V_stern  	Auslegungsgeschwindigkeit
%    rho_stern	Auslegungsh�he bzw. repr�sentiert durch die Dichte
%    V_A    vertikale Steiggeschwindigkeit
%    rho	aktuelle H�he bzw. Dichte
%
% Outputs:
%    Thrust	ben�tigter Schub
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    [ values, derivatives ] = Untitled( x, myStruct.y, 2 )
%    [ values, ~ ] = Untitled( [ 1:.1:100 ] , myStruct.y, 1 )
%
% See also: MULTICOPTERAERODYNAMIK, FLAECHENFLUGZEUGAERODYNAMIK
%
%   Copyright 2018 TU-Braunschweig
% ******************************************************************************
    
    A_stern = m*g;					% Berechnung des Auftriebs im Auslegungszustand, Horizontalflug mit gamma = 0
    
    W_stern = A_stern / E_stern;    % �ber die Auslegungsgleitzahl berechnet sich der Widerstand
    
    W_0_stern = 0.5 * W_stern;		% Der Nullwiderstand ergibt sich f�r einen optimalen Flug als die H�lfte des Gesamtwiderstands
    
    W_0 = W_0_stern * (V_A^2*rho/2)/(V_stern^2*rho_stern/2);	% Skalierung des Nullwiderstandsbeiwertes aus der
    
    Thrust = m*g + W_0;             % Schub besteht aus der Kompensation der Gewichtskraft und des Nullwiderstandes
    
end

