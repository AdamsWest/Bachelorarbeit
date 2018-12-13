% Die größten Abweichungen

run('Fehleruntersuchung_all')
% Initialisierung
Flaechen = zeros(length(Elektromodellflug));


for i = 1:length(Elektromodellflug)
    Flaechen(i) = trapz(tolerance_crate(i,2:end));
end


plot(1:length(Elektromodellflug),Flaechen)
xlabel('Batterienummer (Batterie ID)')
ylabel('Fläche der Abweichung')