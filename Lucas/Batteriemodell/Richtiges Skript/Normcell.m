function [Q,Qnom,Qexp,Vfull,Vexp,Vnom,i,R] = Normcell(DATA)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Initialisierungen
sum_1 = 0;
sum_2 = 0;
sum_3 = 0;
sum_4 = 0;
sum_5 = 0;
sum_6 = 0;
% sum_7 = 0;
sum_8 = 0;

for n = 1:length(DATA)
    
    % Normierung von Q, Q_nom und Q_exp mit der Kapazität in As
    
    DATA{n,3}(1) = DATA{n,3}(1) * 1000 / (DATA{n,5}*3600);
    DATA{n,3}(2) = DATA{n,3}(2) * 1000 / (DATA{n,5}*3600);
    DATA{n,3}(3) = DATA{n,3}(3) * 1000 / (DATA{n,5}*3600);

    % Normzelle: (arithmetischer Mittelwert)
    
    sum_1 = sum_1 + DATA{n,3}(1);
    sum_2 = sum_2 + DATA{n,3}(2);
    sum_3 = sum_3 + DATA{n,3}(3);
    sum_4 = sum_4 + DATA{n,3}(4);
    sum_5 = sum_5 + DATA{n,3}(5);
    sum_6 = sum_6 + DATA{n,3}(6);
    % sum_7 = sum_7 + Elektromodellflug_norm{i,3}(7);
    sum_8 = sum_8 + DATA{n,3}(8);
    

end


% arithmetischer Mittelwert über alle Batterien
% durchschnittliche Kapazität pro Zelle:

Q = sum_1 / length(DATA);
Qnom = sum_2 / length(DATA);
Qexp = sum_3 / length(DATA);
Vfull = sum_4 / length(DATA);
Vexp = sum_5 / length(DATA);
Vnom = sum_6 / length(DATA);
i = 1/100;
R = sum_8 / length(DATA);

% Abgeschätzter Innenwiderstand genauer Regression wäre schön an dieser
% Stelle
R = 0.0160;


end

