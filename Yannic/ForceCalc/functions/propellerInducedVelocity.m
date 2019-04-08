function [ vi, mu_z, vi0 ] = ...
    propellerInducedVelocity( V_A, alpha, m, g, rho, F, n_Prop )
%propellerInducedVelocity computes the induced velocity of a propeller.
%
%   The computation is based on the blade element theory according to
%   [1, S.153].
%
% Literature: 
%   [1] Berend Gerdes Wall. Grundlagen der Hubschrauber-Aerodynamik. 
%       Berlin, Heidelberg: Springer Berlin Heidelberg, 2015. 
%       isbn: 978-3-662-44399-6. 
%
% Syntax [ vi, mu_z, vi0 ] = ...
%    propellerInducedVelocity( V_A, alpha, m, g, rho, F, n_Prop )
%
% Inputs:
%   V_A         absolute velocity of the propeller relative to the air
%               (scalar), in m/s
%   alpha       angle between the rotor plane and the velocity vector 
%               relative to the air (scalar), in rad
%   m           mass of the aircraft (scalar), in kg
%   rho         air density (scalar), in kg/m^3
%   F           surface of the propeller (scalar), in m^2
%   n_Prop      number of propellers (scalar), -
%
% Outputs:
%   vi          induced velocity of the propeller (scalar), in m/s
%   mu_z        velocity of the propeller relative to the air perpendicular
%               to the propeller plane (scalar), in m/s
%   vi0         induced velocity for hover (scalar), in m/s
%
%
%   Copyright 2019 TU Braunschweig
% *************************************************************************

% computation of the induced velocity for hover
vi0 = sqrt(m*g / ( 2*rho*F*n_Prop ) );
% computation of the velocity of the propeller relative to the air
% perpendicular to the propeller plane
mu_z = -V_A*sin(alpha);
% computation of the velocity of the propeller relative to the air
% tangential to the propeller plane
mu = V_A*cos(alpha);
% initial values for the iteration
v = vi0;
krit = 1;
% iteration
while krit > 0.0005
    f = v - mu_z - vi0^2 / sqrt(mu^2 + v^2);
    fs = 1 + v * vi0^2 / (mu^2 + v^2)^(3/2);
    v_i_neu = v - f/fs;
    krit = abs(v_i_neu - v) / v_i_neu;
    v = v_i_neu;
end
% computation of the ratio vi/vi0
vi_vi0 = (v - mu_z) / vi0;
% computation of the induced velocity
vi = vi0 * vi_vi0;
                
end

