function [z_val] = Interpolation(X,Y,Z,x_val,y_val)

%INTERPOLATION doppelte Interpolation für Kennfelder
%X = V
%Y = T
%Z = RPM
%x_val = v_interp
%y_val = t_interp

%   [z_val] = Interpolation(X,Y,Z,x_val,y_val) gibt für einen x- und y-Wert
%   den gewünschten z-Wert aus einem Kennfeld wieder



Y_interp = zeros(length(Z),1);                      %Initialisierung eines Ergebnisvektors der ersten interpolierten Werte
y_row = 0;                                          %Initialisierung der Spalten


% 1. Interpolation

for i = 1:length(Z)                                 %Berechnung des interpolierten Wertes für x_val in jeder Zeile
    y_row = Y(i,:);
    Y_interp(i) = interp1(X,y_row,x_val,'pchip');
end

% 2. Interpolation

z_val = interp1(Y_interp,Z,y_val);

