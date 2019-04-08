function [ rho, a, p, T ] = atmosphere( H, params )
% atmosphere computes the state of the air assuming norm conditions
%
%   The state of the air, density (rho), speed of sound (a), pressure (p),
%   and temperature (T), are computed according to the international
%   standard atmosphere.
%
% Syntax:    [ rho, a, p, T ] = atmosphere( H, params )
%
% Inputs:
%   H           altitude above mean sea level (scalar), m
%   params      a struct containing several parameters
%
% Outputs:
%   rho         air density (scalar), in kg/m^3
%   a           air speed of sound (scalar), m/s
%   p           air pressure (scalar), Pa
%   T           air temperature (scalar), K
%
%
% See also: aerodynamicsMulticopter
%
%   Copyright 2019 TU Braunschweig
% *************************************************************************

% get parameter (conditions at the ground and at the stratopause
T_0 = params.T_0;
T_11 = params.T_11;
rho_0 = params.rho_0;
rho_11 = params.rho_11;
p_0 = params.p_0;
p_11 = params.p_11;
g = params.g;
kappa = params.kappa;

% computation of the physical quantaties depending on the altitude
if H <= 11000
    % troposhere
    rho = rho_0 *(1-0.0065*(H/T_0))^4.256;
    p = p_0 *(1-0.0065*(H/T_0))^5.256;
    T = T_0 - 0.0065 * H;
else
    % stratosphere
    rho = rho_11 * exp(-g/(287*T_11)*(H-11000));
    p = p_11 * exp(-g/(287.1*T_11)*(H-11000));
    T = T_11;
end

% speed of sound, m/s
a = sqrt(kappa*p/rho);

end