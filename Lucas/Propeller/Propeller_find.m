function [ output_args ] = Propeller_find( input_args )
% PROPELLER_FIND sucht in der Propellerdatenbank alle Propeller mit dem
% vorgeschriebenen Durchmesser 

%   ... geht die Propellerdatenbank von oben nach unten durch und 
%   vergleicht dabei die Propellerdurchmesser mit dem angegebenen 
%   Druchmesser. Dabei werden alle Propeller mit anderen Durchmessern aus
%   der Datenbank entfernt.

ind = 1;
len = length(DATA);
while ind < len
    prop_name = DATA{ind,1};
    ind = find(strcmp(DATA(:,1),prop_name));
    ind = max(ind)+1;
    
    [RPM, V, T, P, Tau] = Propeller_map(DATA,prop_name);
    
    counter = counter + 1;
end

end

