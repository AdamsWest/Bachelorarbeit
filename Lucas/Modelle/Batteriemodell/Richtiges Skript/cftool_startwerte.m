% Für cftool

clear
close all
clc

figure_resis = figure;

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

figure(figure_resis)
subplot(121), plot(capacity,resistance,'rx','LineWidth',2), grid on
xlabel('Kapazität [Ah]'), ylabel('Widerstand [Ohm]')
% legend('Batterie')

subplot(122), plot(crate,resistance,'rx','LineWidth',2), grid on
xlabel('C-Rate [1/h]'), ylabel('Widerstand [Ohm]')
% xlim([0 50])
% legend('Batterie')

fig = gcf;
fig.PaperPositionMode = 'auto';
set(fig,'PaperUnits','centimeters', 'PaperPosition', [0 0 21 29.7]); 
fig_pos = fig.PaperPosition;
fig.PaperSize = [21 29.7];
print(fig,'Widerstand','-dpdf')
