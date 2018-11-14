function [z_val] = Interpolation_lin(X,Y,Z,x_val,y_val)


%INTERPOLATION_LIN interpoliert einen Wert z_val für die Werte x_val und 
%y_val für die gegeben X-Vektor, die Matrix Y und den Vektor Z. 


%   function [z_val] = Interpolation_lin(X,Y,Z,x_val,y_val) ...
%   Beschreibung folgt.

% Initialisierung
x_unten = 0;
x_oben = 0;

%Suchen der Position des x_val im X-Vektor                              %Zähler initialisieren
for i = 1:1:length(X)                    
    if x_val <= X(i)                %Laufe Vektor solange durch bis Position von x_val gefunden
        x_oben = X(i);              %Speichere die Position der Werte größer und kleiner als x_val
        x_unten = X(i-1);
        break
    end
end

%Finden und Speicher der Indizes der Werte größer und kleiner als x_val
ind_1 = find(X==x_unten);
ind_2 = find(X==x_oben);

y_1 = Y(:,ind_1);
y_2 = Y(:,ind_2);

Y_inter = zeros(length(Z),1);       %Initialisierung des ersten interpolierten Vektors


%1. Interpolation

for j = 1:1:length(Z)
    Y_inter(j) = y_1(j) + (y_2(j) - y_1(j))*(x_oben-x_val)/(x_oben-x_unten);
    
end

%Initialisierung
y_unten = 0;
y_oben = 0;

%Suchen der Position des y_val im Y-Vektor, analog zur 1. Suche

for k= 1:1:length(Y_inter)
    if y_val <= Y_inter(k)
        y_unten = Y_inter(k-1);
        y_oben = Y_inter(k);
        break
    end
end

%Berechnung des gesuchten Wertes z_val

%Bestimmung Indizes

ind_1 = find(Y_inter == y_unten);
ind_2 = find(Y_inter == y_oben);

z_val = Z(ind_1) + (Z(ind_2) - Z(ind_1))* (y_oben - y_val)/(y_oben - y_unten);

end

