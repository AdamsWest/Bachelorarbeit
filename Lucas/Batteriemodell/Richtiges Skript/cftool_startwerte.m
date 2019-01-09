% Für cftool

clear
close all
clc

load('Elektromodellflug');
% Herausnahme der Kapazität
DATA = Elektromodellflug;
% DATA(63,:) = [];       % id_bat = 63
% DATA(40,:) = [];       % id_bat = 40
DATA(30,:) = [];       % id_bat = 30
% DATA(14,:) = [];       % id_bat = 14
% DATA(38,:) = [];       % id_bat = 38


capacity = zeros(length(DATA),1);
resistance = zeros(length(DATA),1);
crate = zeros(length(DATA),1);

for i = 1:length(DATA)
    
    capacity(i) = DATA{i,5}/1000;
    resistance(i) = DATA{i,3}(end);
    crate(i) = DATA{i,6};
    
end