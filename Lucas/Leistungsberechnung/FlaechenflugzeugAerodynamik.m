function [Thrust,V_A] = FlaechenflugzeugAerodynamik(m,g,epsilon,V_Kg)
% FLAECHENFLUGZEUGAERODYNAMIK berechnet auf Basis eines einfachen
% aerodynamischen Modells den benötigten Schub für ein Flächenflugzeug im
% stationären Steigflug


%   FlaechenflugzeugAerodynamik(m,g,epsilon) bestimmt mit der reziproken
%   Gleitzahl epsilon und der Masse m sowie zusätzlich der Erdbeschleunigung g 
%   den Schub. Dabei bleiben Windeinflüsse unberücksichtig, sodass gilt
%   alpha=sigma=beta= 0. Der Wind hat keinen Einfluss auf das Steigvermögen
%   eines Flugzeugs sondern nur auf die zurückgelegte Strecke über Grund
%   und auf den Steigwinkel, der dennoch ohne Belang in dieser Betrachtung
%   ist.
%   (vgl. Scheiderer, J., „/Angewandte Flugleistung: Eine Einführung in die 
%   operationelle Flugleistung vom Start bis zur Landung/“, Springer, 2008,
%   S. 241f.)

gamma = atan(epsilon);

Thrust = m*g * ( sin(gamma) + epsilon * cos(gamma));

V_A = V_Kg;

end

