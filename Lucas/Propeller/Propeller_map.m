% The function Propeller_map transforms the given propeller map with
% non-uniformly velocity breakpoint vectors to a propeller map where the 
% breakpoint vector is uniformly. This is done by interpolation.
% The inputs of the function are a cell of given propeller maps and the
% name of the propeller.
% The outputs of the function are the RPM vector, the velocity vector V,
% the thrust map (matrix) T and the power map (matrix) P

% Creation of a propeller map by APC Propeller Test Data
function [RPM, V, T, P, Tau] = Propeller_map(DATA_APC,prop_name)

% extract the desired propeller map of the cell of multiple propeller maps
% with aid of the entered propeller name
DATAp = evalin('base','DATA_APC');
ind = find(strcmp(DATAp(:,1),prop_name));
DATAp = DATAp(ind,:);

% definition of the velocity vector V in m/s
x_beg = -10;        % minimum velocity
x_end = 35;         % maximum velocity
delta_x = 1;        % distance of the breakpoints

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
    x_12 = interp1(V,V,x(2),'nearest');     % second element (round)
    x_23 = interp1(V,V,x(end-1),'nearest'); % second last element (round)
    if isnan(x_23)
        x_23 = V(end-1);
    end

    % interpolation of the thrust and power vector with a 10th degree
    % polynomial
    coeff = polyfit(x , y, 10);     % calculate the parameters of the thrust interpolation
    coeffz = polyfit(x , z, 10);    % calculate the parameters of the power interpolation    
    yp_2 = polyval(coeff,x_12:delta_x:x_23);    % calculation of the thrust polynomials
    zp_2 = polyval(coeffz,x_12:delta_x:x_23);   % calculation of the power polynomials

    % calculate the parameters of the derivative of the interpolation
    % polynomials (the derivative is needed for the extrapolation)
    coeff2 = coeff(1:end-1).*[length(coeff)-1:-1:1];
    coeffz2 = coeffz(1:end-1).*[length(coeffz)-1:-1:1];
    
    % calculation of the linear extrapolation functions in both
    % directions (the slope of the linear extrapolation function is
    % calculated with the derivative of the interpolation polynomial)
    yp_1 = sum(x_12.^(0:length(coeff2)-1).*coeff2(end:-1:1)) * ([x_beg:delta_x:x_12-delta_x] - x_12) + yp_2(1);
    yp_3 = sum(x_23.^(0:length(coeff2)-1).*coeff2(end:-1:1)) * ([x_23+delta_x:delta_x:x_end] - x_23) + yp_2(end);
    zp_1 = sum(x_12.^(0:length(coeffz2)-1).*coeffz2(end:-1:1)) * ([x_beg:delta_x:x_12-delta_x] - x_12) + zp_2(1);
    zp_3 = sum(x_23.^(0:length(coeffz2)-1).*coeffz2(end:-1:1)) * ([x_23+delta_x:delta_x:x_end] - x_23) + zp_2(end);
    
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
    T(:,j) = [T(1,j);T1(2:end-2);T(end-1:end,j)];
    P(:,j) = [P(1,j);P1(2:end-2);P(end-1:end,j)];
end

%calculation of the torque vector by stepwise calculating P/(RPM*2*pi/60)
for m = 1:length(V)                                                           
    Tau(:,m) = P(:,m)./(RPM*2*pi/60);
end
