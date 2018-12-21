clear
% close all
clc

load('Elektromodellflug');
% Herausnahme der Kapazität
DATA = Elektromodellflug;
DATA(30,:) = [];

capacity = zeros(length(DATA),1);
resistance = zeros(length(DATA),1);
crate = zeros(length(DATA),1);

for i = 1:length(DATA)
    
    capacity(i) = DATA{i,5}/1000;
    resistance(i) = DATA{i,3}(end);
    crate(i) = DATA{i,6};
    
end

plot3(capacity,crate,resistance,'x')
xlabel('Kapazität in 1/1000Ah')
ylabel('C-Rate')
zlabel('Widerstand')
xlim([0 8]);
ylim([0 50]);

% surf(capacity,crate,resistance)

plot(capacity,resistance,'rx')
xlabel('Kapazität in 1/1000Ah')
ylabel('Widerstand in mOhm')
hold on


fun = @(x,capacity) x(1)./capacity.^1; % Hier Funktion
x0 = 1;             % Startwert festlegen
x = lsqcurvefit(fun,x0,capacity,resistance);

Q = 0:0.1:max(capacity);
func = x./Q;
plot(Q,func)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regression

% Punkte gleicher Kapazität
punkt = x/1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kapazität mit 2200mAh

for i = length(DATA):-1:1
       
    if DATA{i,5} ~= 2200
        DATA(i,:) = [];
    end
end
cratemax_2200 = zeros(length(DATA),1);
resist_2200 = zeros(length(DATA),1);

for i = 1:length(DATA)
       
    cratemax_2200(i) = DATA{i,6};
    resist_2200(i) = DATA{i,3}(end);%/(x/(DATA{i,5}/1000));
    
end

figure
plot(cratemax_2200, resist_2200,'bx')
hold on

fun = @(x_2200,cratemax_2200) x_2200(1)./cratemax_2200.^(1/2); % Hier Funktion
x0 = 1;             % Startwert festlegen
x_2200 = lsqcurvefit(fun,x0,cratemax_2200,resist_2200);

Q = 0:0.1:max(cratemax_2200);
func_2200 = x_2200./Q.^(1/2);
plot(Q,func_2200)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kapazität von 3300
DATA = Elektromodellflug;


for i = length(DATA):-1:1
       
    if DATA{i,5} ~= 3300
        DATA(i,:) = [];
    end
end
dim = size(DATA);
cratemax_3300 = zeros(dim(1),1);
resist_3300 = zeros(dim(1),1);

for i = 1:dim(1)
       
    cratemax_3300(i) = DATA{i,6};
    resist_3300(i) = DATA{i,3}(end);%/(x/(DATA{i,5}/1000));
    
end

figure
plot(cratemax_3300, resist_3300,'gx')
hold on

fun = @(x_3300,cratemax_3300) x_3300(1)./cratemax_3300.^(1/2); % Hier Funktion
x0 = 1;             % Startwert festlegen
x_3300 = lsqcurvefit(fun,x0,cratemax_3300,resist_3300);

Q = 0:0.1:max(cratemax_3300);
func_3300 = x_3300./Q.^(1/2);
plot(Q,func_3300)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kapazität von 3700
DATA = Elektromodellflug;

for i = length(DATA):-1:1
       
    if DATA{i,5} ~= 3700
        DATA(i,:) = [];
    end
end

dim = size(DATA);
cratemax_3700 = zeros(dim(1),1);
resist_3700 = zeros(dim(1),1);

for i = 1:dim(1)
       
    cratemax_3700(i) = DATA{i,6};
    resist_3700(i) = DATA{i,3}(end);%/(x/(DATA{i,5}/1000));
    
end

figure
plot(cratemax_3700, resist_3700,'kx')
hold on

fun = @(x_3700,cratemax_3300) x_3700(1)./cratemax_3700.^(1/2); % Hier Funktion
x0 = 1;             % Startwert festlegen
x_3700 = lsqcurvefit(fun,x0,cratemax_3700,resist_3700);

Q = 0:0.1:max(cratemax_3700);
func_3700 = x_3700./Q.^(1/2);
plot(Q,func_3700)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kapazität von 5000
DATA = Elektromodellflug;

for i = length(DATA):-1:1
       
    if DATA{i,5} ~= 5000
        DATA(i,:) = [];
    end
end

dim = size(DATA);
cratemax_5000 = zeros(dim(1),1);
resist_5000 = zeros(dim(1),1);

for i = 1:dim(1)
       
    cratemax_5000(i) = DATA{i,6};
    resist_5000(i) = DATA{i,3}(end);%/(x/(DATA{i,5}/1000));
    
end

figure
plot(cratemax_5000, resist_5000,'rx')
hold on

fun = @(x_5000,cratemax_5000) x_5000(1)./cratemax_5000.^(1/2); % Hier Funktion
x0 = 1;             % Startwert festlegen
x_5000 = lsqcurvefit(fun,x0,cratemax_5000,resist_5000);

Q = 0:0.1:max(cratemax_5000);
func_5000 = x_5000./Q.^(1/2);
plot(Q,func_5000)
