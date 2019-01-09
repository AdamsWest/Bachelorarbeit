% Die größten Abweichungen
% close all
run('Gesamt')


%% Durchschnittliche Abweichung aller Zellen 

Durchschnitt = zeros(C_Rate_max,1);

for i = 2:1:C_Rate_max+1
    
    Durchschnitt(i-1) = nanmean(tolerance_crate(1:end,i));
    
end

figure
plot(1:C_Rate_max,Durchschnitt)
xlabel('C-Rate');
ylabel('durchschnt. Abweichung aller Zellen in %');



%% Verlauf der Standardabweichung über der C-Rate

Standardabweichung = zeros(C_Rate_max,1);

for i = 2:1:C_Rate_max+1
    
    Standardabweichung(i-1) = std(tolerance_crate(1:end,i),'omitnan');
    
end

figure
plot(1:C_Rate_max,Standardabweichung)
xlabel('C-Rate');
ylabel('Standardabweichung');



%% Copy and Paste Beispiel für eine C-Rate und alle Batterieabweichung

% quick and dirty copy and paste Bespiel Abweichung von C-Rate bei 20
figure
plot(1:length(DATA),tolerance_crate(1:end,21),'rx')
hold on 
bar = zeros(length(DATA),1);
for i = 1:length(DATA)
    bar(i) = nanmean(tolerance_crate(1:end,21));
end
plot(1:length(DATA),bar) 
xlabel('Batterienummer (id\_bat)')
ylabel('Abweichung in %')