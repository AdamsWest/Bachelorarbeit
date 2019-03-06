%% Motorregler

clear
close all
clc

%%
figure_PWM = figure;
% Initialisierung
U_Bat = 14.8;
Delta_U_Bat = 0.1;

lengthi = floor(abs(U_Bat-0)/Delta_U_Bat+1);
PWM = zeros(lengthi,1);
eta_PWM = zeros(lengthi,1);

i = 1;
for U_mot = 0:Delta_U_Bat:U_Bat
    
    [PWM(i),eta_PWM(i)] = ESC(U_mot,U_Bat);
    
end

figure(figure_PWM)
plot(PWM*100,eta_PWM*100)
grid on
xlabel('PWM [%]');
ylabel('\eta_{PWM} [%]');

saveas(gcf,'Motorreglermodell', 'pdf');