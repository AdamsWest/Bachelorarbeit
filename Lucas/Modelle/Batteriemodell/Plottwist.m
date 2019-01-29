
run('Regression')
close all

faktor1 = x;
faktor2 = (x_2200+x_3300+x_3700+x_5000)/4;
faktor2 = faktor2;
kapa = 8000/1000;
entlade = 25;

% Widerstandsberechnung

punkt1 = faktor1/kapa;
punkt2 = (punkt1 + faktor2/entlade)/2;



%%
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

%%

plot3(capacity,crate,resistance,'rx')
xlabel('Kapazität in 1/1000Ah')
ylabel('C-Rate')
zlabel('Widerstand')
xlim([0 8]);
ylim([0 50]);

hold on 
[X Y] = meshgrid(0:0.25:8,0:0.25:50);
Z = (x./X + x_3700./Y.^(1/2))/2;
surface(X,Y,Z)

% plot(kapa,entlade,punkt2,'rx')