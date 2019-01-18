function [K_V, I_0, R_i, m_Mot, S_max, I_max] = Motordata(filename, motor_name)

%MOTORDATA Kennzahlenentnahme
%
%   [K_V, I_0, R_i, m_Mot, S_max, I_max] = Motordata(filename, motor_name)
%   entnimmt Kennzahlen eines Motors einer angegebenen Datenbank. Die Kennzahlen 
%   sind: K_V, I_0,R_i, m_Mot, S_max und I_max

%   Lucas Schreer, 2018


%%%%%%%%%%%%%%%%%%%%%%%% Beginn des Programms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DATAp = evalin('base',filename);                % Speichern aller in der Datenbank enthaltenen Daten unter DATAp
ind = find(strcmp(DATAp(:,1),motor_name));      % Suchen und Speichern des gesuchten Motors sowie der zu ihm gehörigen Spalten
DATAp = DATAp{ind,2};                           % Ablegen des Datenvektors unter DATAp


% Da die Datenbanken eine unterschiedliche Menge an Informationen über den
% jeweiligen Motor enthalten, werden im folgenden die entnommenen
% Datenvektoren auf eine genormte Länge zurückgeführt.

switch filename
    case 'axi_motor_db'
        DAT = DATAp([1 2 3 4 6 7]);
        
    case 'hacker_motor_db'
        DAT = [DATAp(1:3) NaN NaN DATAp(3)];
        
    case 'motocalc_db'
        DAT = [DATAp(1:3) NaN NaN DATAp(3)];
        
    case 'tetacalc_motor_db'
        DAT = [DATAp(1:4) NaN DATAp(5)];
        
    otherwise
        error('databank not found');
        
end

% Belegung der Variablen

K_V = DAT(1);
I_0 = DAT(2);
R_i = DAT(3);
m_Mot = DAT(4);
S_max = DAT(5);
I_max = DAT(6);    
          

end



