function [Eo, A, K] = Batterie_parameter(points)

% the experimental points on the discharge curve
Q = points(1);
Qnom = points(2);
Qexp = points(3);
Vfull = points(4);
Vexp = points(5);
Vnom = points(6);
i = points(7);
R = points(8);

% equation to calculate the 3 parameters Eo, A, K of the discharge curve
% according to Trembley
M_A = [1, 1, 0 ; 1, exp(-3), -Q/(Q-Qexp)*(Qexp+i) ; 1, exp(-3*Qnom/Qexp), -Q/(Q-Qnom)*(Qnom + i)];
b = [Vfull + R*i ; Vexp + R*i ; Vnom + R*i];
x = M_A\b;

Eo = x(1);
A = x(2);
K = x(3);

end