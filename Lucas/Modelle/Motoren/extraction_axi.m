function [K_V, I_0, R_i, m_Mot, S_max, I_max, ges] = extraction_axi(axi_motor_db)
% EXTRACTION Fasst alle Kennzahlen jedes einzelnen Motors in einem
% gemeinsamen Vektor zusammen

%   function [K_V, I_0, R_i, m_Mot, S_max, I_max] = extraction(axi_motor_db)


DATAp = axi_motor_db(:,2);

% Initialisierungen

K_V = zeros(length(DATAp),1); 
I_0 = zeros(length(DATAp),1); 
R_i = zeros(length(DATAp),1); 
m_Mot = zeros(length(DATAp),1);
S_max = zeros(length(DATAp),1);
I_max = zeros(length(DATAp),1);


% Belegung der Einträge der mit den Kennzahlen der Motoren

for i = 1:1:length(DATAp)
    el = DATAp{i};
    
    for j = 1:1:length(el)
        K_V(i) = el(1);
        I_0(i) = el(2);
        R_i(i) = el(3);
        m_Mot(i) = el(4);
        S_max(i) = el(6);
        I_max(i) = el(7);
    end
end

ges = [K_V I_0 R_i m_Mot S_max I_max];          % Alle Infos in einer Matrix

end

