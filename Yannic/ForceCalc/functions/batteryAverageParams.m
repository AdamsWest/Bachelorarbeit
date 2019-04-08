function [ SoC_full, SoC_nom, SoC_exp, V_full, V_exp, V_nom, C_rate ] = ...
    batteryAverageParams( DATA )
% batteryAverageParams   creates average values for the standardised 
%   battery cell
%  
%   For the standardised battery cell the values of the state of charge
%   (SoC) for all three measuring points (full, nom and exp) as well as the
%   voltage (full, exp and nom) are beeing averaged with arithmetic mean.
%   The measuring points were determined at a discharge rate of 10 1/h [1].
%
% Literature: 
%   [1] Gerd Giese. Elektromodellflug - Datenbank. 
%       URL: https://www.elektromodellflug.de/oldpage/datenbank.htm
%       [last downloaded 20.02.2019]. 2012.
%
% Syntax:  [ SoC_full, SoC_nom, SoC_exp, V_full, V_exp, V_nom, C_rate ] = 
%   batteryAverageParams( DATA )
%
% Inputs:
%   DATA        matrix containing the data and four measuring points of the
%               discharge curve of every battery
%
% Outputs:
%   SoC_full    State of Charge with no current applied (scalar), -
%   SoC_nom     nominal State of Charge (scalar), -
%   SoC_exp     exponential State of Charge (scalar), -
%   V_full      cell voltage with no current applied (scalar), in V
%   V_exp       exponential cell voltage (scalar), in V
%   V_nom       nominal cell voltage (scalar), in V
%   C_rate      discharge rate (scalar), in 1/h
%
% See also: batteryDischargeParams, batteryVoltage
%
%   Copyright 2019 TU Braunschweig
% *************************************************************************


% Initialisations 
sum_1 = 0;
sum_2 = 0;
sum_3 = 0;
sum_4 = 0;
sum_5 = 0;
sum_6 = 0;
sum_7 = 0;
len = length(DATA);
SoC_full = zeros(len,1);
SoC_nom = zeros(len,1);
SoC_exp = zeros(len,1);

for n = 1:length(DATA)
    
    % Standardising the battery capacity to the State of Charge 
    % The capacity values are given in As, /3.6 transforms them into
    % mAh, finally the figures are standardised with the overall battery
    % capacity in mAh
    SoC_full(n) = DATA{n,3}(1) / 3.6 / DATA{n,5};
    SoC_nom(n) = DATA{n,3}(2) / 3.6 / DATA{n,5};
    SoC_exp(n) = DATA{n,3}(3) / 3.6 / DATA{n,5};

    % Forming the sum for all values to create the arithmetic mean value
    sum_1 = sum_1 + SoC_full(n);
    sum_2 = sum_2 + SoC_nom(n);
    sum_3 = sum_3 + SoC_exp(n);
    sum_4 = sum_4 + DATA{n,3}(4);
    sum_5 = sum_5 + DATA{n,3}(5);
    sum_6 = sum_6 + DATA{n,3}(6);
    sum_7 = sum_7 + DATA{n,3}(7);  

end

% Computing the arithmetic mean for the State of Charge and cell voltage
% points
SoC_full = sum_1 / length(DATA);
SoC_nom = sum_2 / length(DATA);
SoC_exp = sum_3 / length(DATA);
V_full = sum_4 / length(DATA);
V_exp = sum_5 / length(DATA);
V_nom = sum_6 / length(DATA);

% All points of the discharge curve have been determined at a discharge
% rate (C-Rate) of 10 1/h
C_rate = 10;


end

