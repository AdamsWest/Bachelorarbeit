function [Thrust] = Vertikalflug(m,g,E_stern,v_ver_max,rho_stern,rho)
% FLAECHENFLUGZEUGAERODYNAMIK   berechnet auf Basis eines einfachen
% aerodynamischen Modells den ben�tigten Schub f�r ein Fl�chenflugzeug im
% station�ren Steigflug
%
%   FlaechenflugzeugAerodynamik(m,g,E_stern,V_stern,rho_stern,E,gamma,rho)
%   bestimmt aus dem Auslegungszustand des Flugzeuges im Horizontalflug
%   (E_stern, V_stern und rho_stern) den Schub f�r den Steigflug (gamma > 0).
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
% Syntax:  [Thrust,V_A,Flugzustand_Flaechenflzg] =
%	FlaechenflugzeugAerodynamik(m,g,E_stern,V_stern,rho_stern,E,gamma,rho)
%
%
% Inputs:
%    m		Masse des Fl�chenflugzeugs
%    g		Erdbeschleunigung
%    E_stern	Auslegungsgleitzahl
%    V_stern  	Auslegungsgeschwindigkeit
%    rho_stern	Auslegungsh�he bzw. repr�sentiert durch die Dichte
%    E		Gleitzahl
%    gamma	Bahnneigungswinkel
%    rho	aktuelle H�he bzw. Dichte
%
% Outputs:
%    Thrust	ben�tigter Schub
%    V_A	Fluggeschwindigkeit
%    Flugzustand_Flaechenflzg 	Kontrollvariable des Flugzustandes
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    [ values, derivatives ] = Untitled( x, myStruct.y, 2 )
%    [ values, ~ ] = Untitled( [ 1:.1:100 ] , myStruct.y, 1 )
%
% See also: MULTICOPTERAERODYNAMIK,
%
%   Copyright 2018 TU-Braunschweig
% ******************************************************************************


    
    V_H(i) = v_vert_variabel;
    
    A_stern = m*g;					% Berechnung des Auftriebs im Auslegungszustand, Horizontalflug mit gamma = 0
    
    W_stern = A_stern / E_stern;			% �ber die Auslegungsgleitzahl berechnet sich der Widerstand
    
    W_0_stern = 0.5 * W_stern;			% Der Nullwiderstand ergibt sich f�r einen optimalen Flug als die H�lfte des Gesamtwiderstands
    
    W_0 = W_0_stern * (V_vert_opt^2*rho/2)/(V_stern^2*rho_stern/2);	% Skalierung des Nullwiderstandsbeiwertes aus der
    
    Thrust = m*g + W_0;
    
    
end
end

