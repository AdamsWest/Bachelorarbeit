clear
close all
clc

load('Elektromodellflug');
% Herausnahme der Kapazität

Elektromodellflug(30,:) = [];

capacity = zeros(length(Elektromodellflug),1);
resistance = zeros(length(Elektromodellflug),1);

for i = 1:length(Elektromodellflug)
    
    capacity(i) = Elektromodellflug{i,5}/1000;
    resistance(i) = Elektromodellflug{i,3}(end);
    
end

plot(capacity,resistance,'rx')
xlabel('Kapazität in 1/1000Ah')
ylabel('Widerstand in mOhm')
hold on


fun = @(x,capacity) x(1)./capacity.^1; % Hier Funktion
x0 = [1];             % Startwert festlegen
x = lsqcurvefit(fun,x0,capacity,resistance);

Q = 0:0.1:max(capacity);
func = x./Q;
plot(Q,func)



