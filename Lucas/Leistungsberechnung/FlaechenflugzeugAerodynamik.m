function [Thrust,V_A] = FlaechenflugzeugAerodynamik(m,g,epsilon,V_Kg)
% FLAECHENFLUGZEUGAERODYNAMIK berechnet auf Basis eines einfachen
% aerodynamischen Modells den ben�tigten Schub f�r ein Fl�chenflugzeug im
% station�ren Steigflug


%   FlaechenflugzeugAerodynamik(m,g,epsilon) bestimmt mit der reziproken
%   Gleitzahl epsilon und der Masse m sowie zus�tzlich der Erdbeschleunigung g 
%   den Schub. Dabei bleiben Windeinfl�sse unber�cksichtig, sodass gilt
%   alpha=sigma=beta= 0. Der Wind hat keinen Einfluss auf das Steigverm�gen
%   eines Flugzeugs sondern nur auf die zur�ckgelegte Strecke �ber Grund
%   und auf den Steigwinkel, der dennoch ohne Belang in dieser Betrachtung
%   ist.
%   (vgl. Scheiderer, J., �/Angewandte Flugleistung: Eine Einf�hrung in die 
%   operationelle Flugleistung vom Start bis zur Landung/�, Springer, 2008,
%   S. 241f.)

gamma = atan(epsilon);

Thrust = m*g * ( sin(gamma) + epsilon * cos(gamma));

V_A = V_Kg;

end

