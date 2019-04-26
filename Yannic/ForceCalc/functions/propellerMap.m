function [RPM, V, T, P, Tau] = propellerMap(DATA_APC,prop_name)
%propellerMap   computes regressions of given original APC propeller maps.
%
%   The function Propeller_map transforms the given propeller map with
%   non-uniformly velocity breakpoint vectors to a propeller map where the 
%   breakpoint vector is uniformly. This is done by polynomial regression.
%   The inputs of the function are a cell of given propeller maps and the
%   name of the propeller. Also a linear extrapolation of the APC propeller
%   concerning the airspeed is performed based on the derivative at the
%   borders of the map.
%   The outputs of the function are the RPM vector, the velocity vector V,
%   the thrust map (matrix) T and the power map (matrix) P.
%
% Syntax:  [RPM, V, T, P] = Propeller_map(DATA_APC,prop_name)
%
% Inputs:
% 	 DATA_APC               A matrix of cell array which contains propeller
%                           maps from APC for several propeller types.
%                           (matrix)
%    prop_name              The name of one specific propeller type within
%                           the first column of DATA_APC. (string)
%
% Outputs:
%    RPM                    A vector of the angular velocity of the
%                           propeller in revolutions per minute (RPM). This
%                           vector is extracted from the input DATA_APC and
%                           is not modified. (column vector), in RPM
%    V                      A vector of the velocity (m/s) of the propeller
%                           relative to the (not influenced) air from the 
%                           mininum velocity (first element) to the maximum
%                           velocity (last element). 
%                           The minimum  and maximum values in this vector 
%                           of this vector are defined in the parameters
%                           of this function. The step size between the
%                           minimum and maximum element is also defined in
%                           the parameters of this function. (row vector),
%                           in m/s
%    T                      A matrix of thrust (N) with the dimension of the
%                           length of RPM times the lengths of V. This data
%                           is approximated from the original data with a
%                           polynomial regression in the direction of V.
%                           The order of the polynomial regression is
%                           defined in the parameters of this function.
%                           (matrix), in N
%    P                      A matrix of power (W) with the dimension of the
%                           length of RPM times the lengths of V. This data
%                           is approximated from the original data with a
%                           polynomial regression in the direction of V.
%                           The order of the polynomial regression is
%                           defined in the parameters of this function.
%                           (matrix), in W
%
% 
%   Copyright 2019 TU Braunschweig
% *************************************************************************
 
%% extract the data of the specified propeller 

DATAp = DATA_APC;
ind = find(strcmp(DATAp(:,1),prop_name));
DATAp = DATAp(ind,:);

%% parameters

% order of the polynomial regression
order = 7;
% definition of the velocity vector V in m/s
x_beg = -10;                                % minimum velocity
x_end = round(max(DATAp{end,3})) + 10;	% maximum velocity
delta_x = 1;                                % distance of the breakpoints

%% computation

% properties of the given propeller map
rows = length(DATAp(:,2));
columns = (x_end - x_beg)/delta_x + 1;

% definition of the breakpoint vectors
RPM = [0;cell2mat(DATAp(2:end,2))];
V = x_beg:delta_x:x_end;
% initialization of the thrust and power matrix
T = zeros(rows,columns);
P = zeros(rows,columns);

% the original propeller map data is now interpolated for every given RPM
for k = 1:rows-1
    x = DATAp{k+1,3};   % one element of the V vector (scalar)
    y = DATAp{k+1,4};   % one row of the thrust matrix (vector)
    z = DATAp{k+1,5};   % one row of the power matrix (vector)
    
    % definition of special V breakpoints
    % second element (round)
    x_12 = interp1(V,V,x(2),'nearest') + 1;     
    % second last element (round)
    x_23 = interp1(V,V,x(end-1),'nearest') - 1; 
    if isnan(x_23)
        x_23 = V(end-1);
    end

    % interpolation of the thrust and power vector with a 5th degree
    % (order) polynomial
    % calculate the parameters of the thrust interpolation
    coeff = polyfit(x , y, order);     
    % calculate the parameters of the power interpolation
    coeffz = polyfit(x , z, order);        
    % calculation of the thrust polynomials
    yp_2 = polyval(coeff,x_12:delta_x:x_23);   
    % calculation of the power polynomials
    zp_2 = polyval(coeffz,x_12:delta_x:x_23);   

    % calculate the parameters of the derivative of the interpolation
    % polynomials (the derivative is needed for the extrapolation)
    coeff2 = coeff(1:end-1).*[length(coeff)-1:-1:1];
    coeffz2 = coeffz(1:end-1).*[length(coeffz)-1:-1:1];
    
    % calculation of the linear extrapolation functions in both
    % directions (the slope of the linear extrapolation function is
    % calculated with the derivative of the interpolation polynomial)
    yp_1 = sum(x(2).^(0:length(coeff2)-1).*coeff2(end:-1:1)) * ...
        ([x_beg:delta_x:x_12-delta_x] - x_12) + yp_2(1);
    yp_3 = sum(x(end-1).^(0:length(coeff2)-1).*coeff2(end:-1:1)) * ...
        ([x_23+delta_x:delta_x:x_end] - x_23) + yp_2(end);
    zp_1 = sum(x(2).^(0:length(coeffz2)-1).*coeffz2(end:-1:1)) * ...
        ([x_beg:delta_x:x_12-delta_x] - x_12) + zp_2(1);
    zp_3 = sum(x(end-1).^(0:length(coeffz2)-1).*coeffz2(end:-1:1)) * ...
        ([x_23+delta_x:delta_x:x_end] - x_23) + zp_2(end);
    
    % calculation of the complete thrust and power vector which consists of
    % the interpolation polynomial and the linear extrapolations in both
    % directions
    T(k+1,:) = [yp_1, yp_2, yp_3];
    P(k+1,:) = [zp_1, zp_2, zp_3];
end

% as the creation of the thrust and power map was achieved by stepwise 2D
% interpolations in V direction, the obtained maps are not necessarily
% smooth in RPM direction. The smooth function is used to smooth the maps
% in RPM direction.
for j = 1:length(T(1,:))
    T1 = smooth(T(:,j));
    P1 = smooth(P(:,j));
%     T1 = T(:,j);
%     P1 = P(:,j);
    T(:,j) = [T(1,j);T1(2:end-2);T(end-1:end,j)];
    P(:,j) = [P(1,j);P1(2:end-2);P(end-1:end,j)];
end

%calculation of the torque vector by stepwise calculating P/(RPM*2*pi/60)
for m = 1:length(V)                                                          
    Tau(:,m) = P(:,m)./(RPM*2*pi/60);
end