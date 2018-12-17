% Die größten Abweichungen

run('Gesamt')
% Initialisierung
Flaechen = zeros(length(DATA));


for i = 1:length(DATA)
    Flaechen(i) = trapz(tolerance_crate(i,2:end));
end

% figure
% plot(1:length(DATA),Flaechen)
% xlabel('Batterienummer (Batterie ID)')
% ylabel('Fläche der Abweichung')

% Darstellung der Durchscnittsabweichungen

Durchschnitt = zeros(C_Rate_max,1);

for i = 2:1:C_Rate_max+1

    Durchschnitt(i-1) = mean(tolerance_crate(1:end,i));
end

figure
plot(1:C_Rate_max,Durchschnitt)
xlabel('C-Rate');
ylabel('durchschnt. Abweichung aller Zellen');


plot(1:length(DATA),tolerance_crate(1:end,21))
hold on 
bar = zeros(length(DATA),1);
for i = 1:length(DATA)
    bar(i) = mean(tolerance_crate(1:end,21));
end
plot(1:length(DATA),bar) 