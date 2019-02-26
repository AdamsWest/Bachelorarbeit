% Die größten Abweichungen
clear
close all
clc

figure_durchschn_abweichung = figure;
run('Gesamt')


%% Durchschnittliche Abweichung aller Zellen 

Durchschnitt = zeros(C_Rate_max,1);

for i = 2:1:C_Rate_max+1
    
    Durchschnitt(i-1) = nanmean(tolerance_crate(1:end,i));
    
end

figure(figure_durchschn_abweichung)
subplot(121)
plot(1:C_Rate_max,Durchschnitt,'LineWidth',2)
grid on
xlabel('Entladerate [1/h]');
ylabel('Ø Batteriezellenabweichung [%]');
% flaeche = trapz(abs(Durchschnitt));
% figure
% plot(1:C_Rate_max,abs(Durchschnitt))
% xlabel('C-Rate');
% ylabel('durchschnt. Abweichung aller Zellen in %');
% disp(num2str(flaeche));


%% Verlauf der Standardabweichung über der C-Rate

% Standardabweichung = zeros(C_Rate_max,1);
% 
% for i = 2:1:C_Rate_max+1
%     
%     Standardabweichung(i-1) = std(tolerance_crate(1:end,i),'omitnan');
%     
% end
% 
% figure
% plot(1:C_Rate_max,Standardabweichung)
% xlabel('C-Rate');
% ylabel('Standardabweichung');



%% Copy and Paste Beispiel für eine C-Rate und alle Batterieabweichung

% quick and dirty copy and paste Bespiel Abweichung von C-Rate bei 20
subplot(122)
plot(1:length(DATA),tolerance_crate(1:end,21),'bx','LineWidth',1.5)
grid on
hold on 
bar = zeros(length(DATA),1);
for i = 1:length(DATA)
    bar(i) = nanmean(tolerance_crate(1:end,21));
end
plot(1:length(DATA),bar,'r', 'LineWidth',1.5) 
xlabel('Batterienummer','LineWidth',2)
ylabel('Abweichung der Zellenspannung [%]')
lgd = legend('Abweichung einer Batteriezelle','Durchschnittliche Abweichung'); title(lgd,'Spannungsabweichung')

% ImagesizeX
% ImagesizeY
% print('-bestfit','Abweichungen','-dpdf')

ImageSizeX = 14;
ImageSizeY = 20;
figure(figure_durchschn_abweichung)
set(gcf,'PaperUnits','centimeters', 'PaperPosition', [0 0 ImageSizeX ImageSizeY]); 
set(gcf,'Units','centimeters', 'PaperSize', [ImageSizeX ImageSizeY]); 
saveas(gcf,'Abweichungen', 'pdf');  