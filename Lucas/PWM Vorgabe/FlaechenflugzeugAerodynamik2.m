function [Thrust,V_A] = FlaechenflugzeugAerodynamik2(u_Wg,v_Wg,w_Wg,V_Kg,epsilon,S,rho,m,g)

% AERODYNAMIK_FLAECHENFLUGZEUG2 berechnet den Schub für ein Flaechenflugzeug auf statische Weise

%  [Thrust,V_A] =
%  Aerodynamik_Flaechenflugzeug(u_Wg,v_Wg,w_Wg,V_Kg,epsilon,S,rho,m,g)
%  berechnet über die Windgeschwindigkeit und die Bahngeschwindigkeit die
%  am Flugzeug angreifenden Kräfte. Diese werden mittels
%  Koordinatentransformation in das geodätische überführt und anschließend
%  der benötigte Schub berechnet.




% nach G. Brüning, X.Hafer, "Flugleistungen - Grundlagen, Flugzustände,
% Flugabschnitte" 1983, Springer-Verlag, Berlin Heidelberg New York: 
% S. 77
% "Da der Winkel (alpha + sigma) für konventionelle Flugzeuge sehr klein
% ist [...] kann man mit hinreichender Genauigkeit cos(alpha + sigma) = 1,
% F*sin(alpha + sigma) << A; mg*cos(gamma) setzen."



% Windgeschwindigkeiten null setzen, wenn sie nicht vorhanden sind

if u_Wg == 0
	u_Wg = 0.01;
end

if u_Wg == 0
	u_Wg = 0.01;
end

if u_Wg == 0
	u_Wg = 0.01;
end

gamma = atan(epsilon);

u_Kg = cos(gamma) * V_Kg;
w_Kg = -sin(gamma) * V_Kg;

V_A = sqrt((u_Kg + u_Wg)^2 + (v_Wg)^2 + (w_Kg + u_Wg)^2);           % Absolute Fluggeschwindigkeit
gamma_a = atan( (w_Kg + w_Wg) / sqrt((u_Kg + u_Wg)^2 + v_Wg));		% Windanstellwinkel
chi_a = atan(v_Wg / (u_Kg  + u_Wg));                                % Flugwindazimut

A = cos(gamma) * m * g;                                             % aus GGW der Kräfte
c_A = A / (rho/2 * S * V_A^2);
c_W = epsilon * c_A;
W = c_W * rho/2 * V_A^2 * S;                                
X_g = - W * cos(gamma_a) * cos(chi_a) + A * sin(gamma_a);           % Koordinatentransformation
Y_g = W * sin(chi_a);
Z_g = - W *sin(gamma_a) * cos(chi_a) - A * cos(gamma_a) + m*g;

Thrust = sqrt(X_g^2 + Y_g^2 + Z_g^2);                               % Schubberechnung

end