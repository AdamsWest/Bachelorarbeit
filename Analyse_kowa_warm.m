%% clean and load

close all;
clear all;
clc;
tic;

farbe.dgruen = [9 121 52] ./ 255;
farbe.rot = [224 1 13] ./ 255;
farbe.dblau = [50 11 155] ./ 255;
farbe.gelb = [241 226 0] ./ 255;
farbe.hblau = [46 194 240] ./ 255;
farbe.orange = [255 150 4] ./ 255;
farbe.hgruen = [33 211 43] ./ 255;
farbe.lila = [193 38 233] ./ 255;

cd('Korrigierte Daten');
Daten = dir;
cd('..');

cd('Abschnitte');
Abschnitte = dir;
cd('..');

deselected_flight = [];

for f = 3:length(Daten)
  
  if ~any((deselected_flight +2)==f)
    
    cd('Korrigierte Daten');
    load(Daten(f).name);        % Daten aus Korrigierte Daten holen
%     cellfun(@(x) evalin('base',[x '=B.' x ]),fieldnames(B));
%     clear B;
    cd('..');
    cd('Abschnitte');
    load(Abschnitte(f).name); % Abschnitte laden
    cd('..');
    
    for i=1:7                   % Fenster für Flugabschnitte generieren
    evalc([['x' int2str(i)] '=[Abschnitt(i).index(1);Abschnitt(i).index(1);Abschnitt(i).index(end);Abschnitt(i).index(end)]'] ); 
    end
    
    mb =(100: length(B.Bus_TOW)); % ersten 100 Werte abgeschnitten
    B.Bus_TOW = B.Bus_TOW(mb);
    B.h_a= B.h_a(mb);
    B.p = B.p(mb);
    bereich = [mb];


 r.tsy = B.tsy(bereich);
 r.fw0 = B.fw0(bereich) - mean(B.fw0(bereich)-B.tsy(bereich));
 r.fw1 = B.fw1(bereich) - mean(B.fw1(bereich)-B.tsy(bereich));
 r.pt0 = B.pt0(bereich) - mean(B.pt0(bereich)-B.tsy(bereich));
 r.pt1 = B.pt1(bereich) - mean(B.pt1(bereich)-B.tsy(bereich));
 r.humi = B.humi(bereich) - mean(B.humi(bereich)-B.tsy(bereich));
 r.ir = B.ir(bereich) - mean(B.ir(bereich)-B.tsy(bereich));
%  r.humiroh = B.humiroh(bereich) - mean(B.humiroh(bereich)-B.tsy(bereich));
        %% Eingangsmessdaten   Temperatur über die Zeit
  y=[-60;max(B.ir)+2;max(B.ir)+2;-60];
  
    %Tiefpass
 fs = 100;
 N    = 4;
 f_TP = (2);
 [Bt, Wt] = butter(N, f_TP /(fs / 2),'low');
 
  B.h_a =  filtfilt(Bt, Wt, B.h_a);
  
  
   q=figure('visible','on','Name',strcat('Eingangsdaten'));
   set(q, 'Units', 'normalized', 'Position',[0.1, 0.1, 0.8, 0.8] );
   AxesH = axes('Units', 'normalized', 'Position', [0.082,0.135, 0.81, 0.74]);
   p(1)=plot(B.Bus_TOW,r.fw0, 'color' ,farbe.dblau,'LineWidth', 2);
   hold on;
   p(2)=plot(B.Bus_TOW,r.fw1, 'color' ,farbe.dgruen,'LineWidth', 2);
   hold on;
   p(3)=plot(B.Bus_TOW,r.pt0, 'color' ,farbe.rot,'LineWidth', 2);
   hold on;
   p(4)=plot(B.Bus_TOW,r.pt1, 'color' ,farbe.gelb,'LineWidth', 2);
   hold on;
   p(5)=plot(B.Bus_TOW,r.tsy, 'color' ,farbe.lila,'LineWidth', 2);
   hold on;
   p(6)=plot(B.Bus_TOW,r.humi, 'color' ,farbe.orange,'LineWidth', 2);
   hold on;
   [hAx,hLine1,hLine2] = plotyy(B.Bus_TOW,r.ir,B.Bus_TOW,B.h_a );
   hold on;
   hLine1.Color = farbe.hblau;
   hLine2.Color = 'k';
   hLine1.LineWidth = 2;
   hLine2.LineWidth = 2;
   set(hAx,'YColor','k', 'FontSize', 30);
   h = legend( 'Finewire-0','Finewire-1','Pt-100-Sensor-0','Pt-100-Sensor-1','TSYS01','Humicap HMP-110','Infrarotsensor','Höhe','location','southwest');
   set(h,'FontSize',25)
   xlabel('Sekunden des Tages [s]',  'FontSize', 30);
   ylabel('Temperatur [°C]', 'FontSize', 40,'color' ,'k');
   ylabel((hAx(2)),'Höhe [m]','color' ,'k', 'FontSize', 30);
   set(gca, 'FontSize', 30)
   xx=fill(B.Bus_TOW(x3),y,'r',B.Bus_TOW(x4),y,'r');
   set(xx,'facealpha',.08);
   ylim([6 max(B.ir)+2]);
   ylim((hAx(2)),[-60 max(B.h_a)*3]);
  % title(strcat({['Eingangsmessdaten Temperatur Flug ' int2str(f-2)] ; ['ALICE Zarnekow']})); 
   grid on; zoom on;
   
     fig = gcf;
     fig.PaperUnits = 'inches';
     fig.PaperPosition = [0 0 22 13];
     print(fig,'-dpdf','-r0',['F:\ALICE_Zarnekow_2018-05-23\Plots\Eingangsmessdaten\Eingangsmessdaten_T(t)_V_F' int2str(f-2) '.pdf'])
%  
   cd('Plots');
%    cd('Eingangsmessdaten');
%    saveas(gcf,['Eingangsmessdaten_T(t)_F' int2str(f-2)],'png');
%    cd('../..'); 
    


        %% Eingangsmessdaten   Höhe über die Temperatur

        
        
%    q=figure('visible','on','Name',strcat('Eingangsdaten'));
%    set(q, 'Units', 'normalized', 'Position',[0.1, 0.1, 0.8, 0.8] );
%    AxesH = axes('Units', 'normalized', 'Position', [0.1,0.135, 0.81, 0.74]);
%    plot(r.fw0,B.h_a, 'color' ,farbe.dblau,'LineWidth', 2);
%    hold on;
%    plot(r.fw1,B.h_a, 'color' ,farbe.dgruen,'LineWidth', 2);
%    hold on;
%    plot(r.pt0,B.h_a, 'color' ,farbe.rot,'LineWidth', 2);
%    hold on;
%    plot(r.pt1,B.h_a, 'color' ,farbe.gelb,'LineWidth', 2);
%    hold on;
%    plot(r.tsy,B.h_a, 'color' ,farbe.lila,'LineWidth', 2);
%    hold on;
%    plot(r.humi,B.h_a, 'color' ,farbe.orange,'LineWidth', 2);
%    hold on;
%    plot(r.ir,B.h_a, 'color' ,farbe.hblau,'LineWidth', 2);
%    set(gca, 'FontSize', 30)
%    title(strcat({['Temperaturgradient Flug ' int2str(f-2)] ; ['ALICE Zarnekow']})); 
%    h = legend( 'Finewire-0','Finewire-1','Pt-100-Sensor-0','Pt-100-Sensor-1','TSYS01','Humicap HMP110','Infrarotsensor');
%    set(h,'FontSize',25)
%    xlabel('Temperatur [°C]', 'FontSize', 30);
%    ylabel('Höhe [m]', 'FontSize', 30);
%    grid on; zoom on;
%    
%      fig = gcf;
%      fig.PaperUnits = 'inches';
%      fig.PaperPosition = [0 0 22 13];
%      print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Eingangsmessdaten\Eingangsmessdaten_T(h)_F' int2str(f-2) '.png'])

    
 %% Filterung
 
 
 
 %Tiefpass
 fs = 100;
 N    = 4;
 f_TP = (2);
 [Bt, Wt] = butter(N, f_TP /(fs / 2),'low');
 
 tp.fw0 =  filtfilt(Bt, Wt, r.fw0); 
 tp.fw1 =  filtfilt(Bt, Wt, r.fw1); 
 tp.pt0 =  filtfilt(Bt, Wt, r.pt0); 
 tp.pt1 =  filtfilt(Bt, Wt, r.pt1); 
 tp.tsy =  filtfilt(Bt, Wt, r.tsy); 
 tp.humi =  filtfilt(Bt, Wt, r.humi); 
 tp.ir =  filtfilt(Bt, Wt, r.ir); 

%    q=figure('visible','on','Name',strcat('Eingangsdaten'));
%    set(q, 'Units', 'normalized', 'Position',[0.1, 0.1, 0.8, 0.8] );
%    AxesH = axes('Units', 'normalized', 'Position', [0.1,0.135, 0.81, 0.74]);
%    plot(tp.fw0,B.h_a, 'color' ,farbe.dblau,'LineWidth', 2);
%    hold on;
%    plot(tp.fw1,B.h_a, 'color' ,farbe.dgruen,'LineWidth', 2);
%    hold on;
%    plot(tp.pt0,B.h_a, 'color' ,farbe.rot,'LineWidth', 2);
%    hold on;
%    plot(tp.pt1,B.h_a, 'color' ,farbe.gelb,'LineWidth', 2);
%    hold on;
%    plot(tp.tsy,B.h_a, 'color' ,farbe.lila,'LineWidth', 2);
%    hold on;
%    plot(tp.humi,B.h_a, 'color' ,farbe.orange,'LineWidth', 2);
%    hold on;
%    plot(tp.ir,B.h_a, 'color' ,farbe.hblau,'LineWidth', 2);
%    title(strcat({['Temperaturgradient   -   2 Hz  tiefpassgefiltert'];[  'Flug '  int2str(f-2) ', ALICE Zarnekow']})); 
%    h = legend( 'Finewire-0','Finewire-1','Pt-100-Sensor-0','Pt-100-Sensor-1','TSYS01','Humicap HMP110','Infrarotsensor');
%    set(h,'FontSize',25)
%    xlabel('Temperatur [°C]', 'FontSize', 30);
%    ylabel('Höhe [m]', 'FontSize', 30);
%    set(gca, 'FontSize', 30)
%    grid on; zoom on;  
%    
%    fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_gefiltert\TP_Daten_T(h)_F' int2str(f-2) '.png'])


%% Spektrum Eingangsdaten (Temperatur)  Steigflug(Flug_7) = 1121 Punkte


%  %Hochpass
%  fs = 100;
%  N    = 4;
%  f_HP = (1/11);
%  [Bh, Wh] = butter(N, f_HP /(fs / 2),'high');
%  
%  r.fw0 =  filtfilt(Bh, Wh, r.fw0); 
%  r.fw1 =  filtfilt(Bh, Wh, r.fw1); 
%  r.pt0 =  filtfilt(Bh, Wh, r.pt0); 
%  r.pt1 =  filtfilt(Bh, Wh, r.pt1); 
%  r.tsy =  filtfilt(Bh, Wh, r.tsy); 
%  r.humi =  filtfilt(Bh, Wh, r.humi); 
%  r.ir =  filtfilt(Bh, Wh, r.ir); 

%  r.fw0 =  detrend(r.fw0, 'linear');
%  r.fw1 =  detrend(r.fw1, 'linear');
%  r.pt0 =  detrend(r.pt0, 'linear'); 
%  r.pt1 =  detrend(r.pt1, 'linear');
%  r.tsy =  detrend(r.tsy, 'linear'); 
%  r.humi =  detrend(r.humi, 'linear'); 
%  r.ir =  detrend(r.ir, 'linear'); 
%  B.humirohd =  detrend(B.humiroh, 'linear'); 

for i = 3:4
    
    
    nx = max(size(Abschnitt(i).index));
    na = 2;                        %%%%%%%%%% Glättungsfaktor 1..20    
    ww = hamming(2048);
    ww2 = hamming(64);
   
    if mod((length(Abschnitt(i).index)),(length(ww)))>= length(ww)/2
    fenster = floor(length(Abschnitt(i).index)/(length(ww)))+1;
    else
    fenster = floor(length(Abschnitt(i).index)/(length(ww)));
    end
    
    
    Abschnitth = [Abschnitt(i).index(1)+99:100:(Abschnitt(i).index(end))];
    Abschnitth = floor(Abschnitth./100);
    
    
     r.fw0d(Abschnitt(i).index) =  detrend(r.fw0(Abschnitt(i).index), 'linear');
     r.fw1d(Abschnitt(i).index) =  detrend(r.fw1(Abschnitt(i).index), 'linear');
     r.pt0d(Abschnitt(i).index) =  detrend(r.pt0(Abschnitt(i).index), 'linear'); 
     r.pt1d(Abschnitt(i).index) =  detrend(r.pt1(Abschnitt(i).index), 'linear');
     r.tsyd(Abschnitt(i).index) =  detrend(r.tsy(Abschnitt(i).index), 'linear'); 
     r.humid(Abschnitt(i).index) =  detrend(r.humi(Abschnitt(i).index), 'linear'); 
     r.ird(Abschnitt(i).index) =  detrend(r.ir(Abschnitt(i).index), 'linear'); 
     B.humirohd(Abschnitth) =  detrend(B.humiroh(Abschnitth), 'linear'); 
    
  
   
   
   
   [Tspektrum.P1, F1] = pwelch(r.fw0d(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P2, F2] = pwelch(r.fw1d(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P3, F3] = pwelch(r.pt0d(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P4, F4] = pwelch(r.pt1d(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P5, F5] = pwelch(r.tsyd(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P6, F6] = pwelch(B.humirohd(Abschnitth), ww2, [], [], 1);
   [Tspektrum.P7, F7] = pwelch(r.ird(Abschnitt(i).index), ww, [], [], 100);
   
   
   freq = 0.01:0.01:50;
   f2 = 0.01:0.01:0.1;
   fl  = (freq) .^ (-5/3) .* 1e-3;  
   fk  = (freq) .^ (-5/3) .* 1e-5;  
   fm  = (f2) .^ (-4)   .* 1e-7;  

% 
%    q=figure('visible','off','Name',strcat('Spektren Mischungsverhältnis'));
%    set(q, 'Units', 'normalized', 'Position',[0.1, 0.1, 0.8, 0.8] );
%    loglog( F1 , Tspektrum.P1, 'color' ,farbe.dblau,'LineWidth',2);
%    hold on;
%    loglog( F2 , Tspektrum.P2, 'color' ,farbe.dgruen,'LineWidth',2);
%    hold on;
%    loglog( F3 , Tspektrum.P3, 'color' ,farbe.rot,'LineWidth',2);
%    hold on;
%    loglog( F4 , Tspektrum.P4, 'color' ,farbe.gelb,'LineWidth',2);
%    hold on;
%    loglog( F5 , Tspektrum.P5, 'color' ,farbe.lila,'LineWidth',2);
%    hold on;
%    loglog( F6 , Tspektrum.P6, 'color' ,farbe.orange,'LineWidth',2);
%    hold on;
%    loglog( F7 , Tspektrum.P7, 'color' ,farbe.hblau,'LineWidth',2);
%    hold on;
%    loglog( freq , fl, 'k--');
%    hold on;
%    loglog( freq , fk, 'k--');
%    hold on;
%    title(strcat({' Leistungsdichtespektrum Temperatur 2048 DFT-Punkte (64 Humicap)' ,[Abschnitt(i).name ' Flug-' int2str(f-2) ' ALICE Zarnekow, ' int2str(fenster) ' Fenster']}), 'FontSize', 40); 
%    h=legend( 'Finewire-0','Finewire-1','Pt-100-Sensor-0','Pt-100-Sensor-1','TSYS01','Humicap HMP110','Infrarotsensor');
%    set(h,'location','northeast','FontSize',25)
%    xlabel('Frequenz [Hz]', 'FontSize', 35);
%    ylabel('Temperatur [°C]^2s', 'FontSize', 35);   
%    set(gca, 'FontSize', 30)
%    grid on; zoom on;
%    xlim([1e-3 60]);
%    ylim([1e-10 1e4]); 
% 
%    fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flugphase\Phase_' int2str(i) '\Spektrum_Eingang\Spektrum_Hamming2048_roh_F' int2str(f-2) '.png'])
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flug\Flug_' int2str(f-2) '\Spektrum_Eingang\Spektrum_Hamming2048_roh_Phase' int2str(i) '.png'])

   
%  close all;

end
%% Kohärenz- und Phasenspektrum (Temperatur)

for i = 3:4
   nx = max(size(Abschnitt(i).index));
   na = 2;                        %%%%%%%%%% Glättungsfaktor 1..20    
   ww = hamming(2048);
   
    if mod((length(Abschnitt(i).index)),(length(ww)))>= length(ww)/2
    fenster = floor(length(Abschnitt(i).index)/(length(ww)))+1;
    else
    fenster = floor(length(Abschnitt(i).index)/(length(ww)));
    end
    
   [P1coh, F1] = mscohere(r.fw1d(Abschnitt(i).index),r.fw0(Abschnitt(i).index),  ww, 0, [], 100); 
   [P2coh, F2] = mscohere(r.pt0d(Abschnitt(i).index),r.fw0(Abschnitt(i).index),  ww, 0, [], 100);
   [P3coh, F3] = mscohere(r.pt1d(Abschnitt(i).index),r.fw0(Abschnitt(i).index),  ww, 0, [], 100);
   [P4coh, F4] = mscohere(r.tsyd(Abschnitt(i).index),r.fw0(Abschnitt(i).index),  ww, 0, [], 100);
   [P5coh, F5] = mscohere(r.humid(Abschnitt(i).index),r.fw0(Abschnitt(i).index), ww, 0, [], 100);
   [P6coh, F6] = mscohere(r.ird(Abschnitt(i).index),r.fw0(Abschnitt(i).index),   ww, 0, [], 100);
   
   nb = 3;
   
   P1cohs = sgolayfilt(P1coh,1,nb);
   P2cohs = sgolayfilt(P2coh,1,nb);
   P3cohs = sgolayfilt(P3coh,1,nb);
   P4cohs = sgolayfilt(P4coh,1,nb);
   P5cohs = sgolayfilt(P5coh,1,nb);
   P6cohs = sgolayfilt(P6coh,1,nb);
   
   freq1 = find(P1coh >= 0.5); 
   freq2 = find(P2coh >= 0.5);
   freq3 = find(P3coh >= 0.5);
   freq4 = find(P4coh >= 0.5);
   freq5 = find(P5coh >= 0.5);
   freq6 = find(P6coh >= 0.5);
   
   f50 = ones(length(F1), 1) .* 0.5;
   
   
   
%    q=figure('visible','off','Name',strcat('Spektren Kohärenz'));
%    set(q, 'Units', 'normalized', 'Position',[0.1, 0.1, 0.8, 0.8] );
%    AxesH = axes('Units', 'normalized', 'Position', [0.1,0.135, 0.87, 0.74]);
%    semilogx(F1 , P1cohs,'color', farbe.dgruen,'LineWidth',2);hold on;
%    semilogx(F2 , P2cohs,'color', farbe.rot,'LineWidth',2);hold on;
%    semilogx(F3 , P3cohs,'color', farbe.gelb,'LineWidth',2);hold on;
%    semilogx(F4 , P4cohs,'color', farbe.lila,'LineWidth',2);hold on;
%    semilogx(F5 , P5cohs,'color', farbe.orange,'LineWidth',2);hold on;
%    semilogx(F6 , P6cohs,'color', farbe.hblau,'LineWidth',2);hold on;
%    semilogx(F1 , f50,'color', farbe.rot);
%    title(strcat( {[' Kohärenzspektrum 2048 DFT-Punkte - Auf FW0 '];[ Abschnitt(i).name, ' Flug ', int2str(f-2) ' ALICE Zarnekow, ' int2str(fenster) ' Fenster']}), 'FontSize', 40); 
%    h=legend('Finewire-1','Pt-100-Sensor-0','Pt-100-Sensor-1','TSYS01','Humicap HMP110','Infrarotsensor');
%    set(h,'location','southwest','FontSize',25)
%    xlabel('f [Hz]', 'FontSize', 35); 
%    ylabel('Kohärenz', 'FontSize', 35); 
%    set(gca, 'FontSize', 30)
%    grid on; zoom on;
%    ylim([0 1]);
%    xlim([1e-2 60]);
% 
% 
%   
%    fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flugphase\Phase_' int2str(i) '\Spektrum_Kohärenz\Spektrum_Hamming2048_SG_' int2str(nb) '_F' int2str(f-2) '.png'])
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flug\Flug_' int2str(f-2) '\Spektrum_Kohärenz\Spektrum_Hamming2048_SG_' int2str(nb) '_Phase' int2str(i) '.png'])



  clear nb;
   [Ph1, F1] = cpsd(r.fw1d(Abschnitt(i).index),r.fw0(Abschnitt(i).index), ww, 0, [], 100); 
   [Ph2, F2] = cpsd(r.pt0d(Abschnitt(i).index),r.fw0(Abschnitt(i).index), ww, 0, [], 100);
   [Ph3, F3] = cpsd(r.pt1d(Abschnitt(i).index),r.fw0(Abschnitt(i).index), ww, 0, [], 100);
   [Ph4, F4] = cpsd(r.tsyd(Abschnitt(i).index),r.fw0(Abschnitt(i).index), ww, 0, [], 100);
   [Ph5, F5] = cpsd(r.humid(Abschnitt(i).index),r.fw0(Abschnitt(i).index), ww, 0, [], 100);
   [Ph6, F6] = cpsd(r.ird(Abschnitt(i).index),r.fw0(Abschnitt(i).index), ww, 0, [], 100);

   nb = 3;
   
%    Ph1 = sgolayfilt(Ph1,1,nb);
%    Ph2 = sgolayfilt(Ph2,1,nb);
%    Ph3 = sgolayfilt(Ph3,1,nb);
%    Ph4 = sgolayfilt(Ph4,1,nb);
%    Ph5 = sgolayfilt(Ph5,1,nb);
%    Ph6 = sgolayfilt(Ph6,1,nb);
   
   
   
%    q=figure('visible','off','Name',strcat('Spektren Phase'));
%    set(q, 'Units', 'normalized', 'Position',[0.1, 0.1, 0.8, 0.8] );
%    AxesH = axes('Units', 'normalized', 'Position', [0.1,0.155, 0.87, 0.70]);
%     scatter( F1 , angle(Ph1).*(180/pi),70.*P1coh+1,farbe.dgruen,'filled', 'MarkerFaceColor', farbe.dgruen);hold on;
%     scatter( F2 , angle(Ph2).*(180/pi),70.*P2coh+1,farbe.rot,'filled', 'MarkerFaceColor', farbe.rot);hold on;
%     scatter( F3 , angle(Ph3).*(180/pi),70.*P3coh+1,farbe.gelb,'filled', 'MarkerFaceColor', farbe.gelb);hold on;
%     scatter( F4 , angle(Ph4).*(180/pi),70.*P4coh+1,farbe.lila,'filled', 'MarkerFaceColor', farbe.lila);hold on;
%     scatter( F5 , angle(Ph5).*(180/pi),70.*P5coh+1,farbe.orange,'filled', 'MarkerFaceColor', farbe.orange);
%     scatter( F6 , angle(Ph6).*(180/pi),70.*P6coh+1,farbe.hblau,'filled', 'MarkerFaceColor', farbe.hblau);
%    title(strcat({[' Phasenspektrum 2048 DFT-Punkte - Auf FW0 '];[ Abschnitt(i).name, ' Flug ', int2str(f-2) ' ALICE Zarnekow, ' int2str(fenster) ' Fenster']}), 'FontSize', 40); 
%    h=legend('Finewire-1','Pt-100-Sensor-0','Pt-100-Sensor-1','TSYS01','Humicap HMP110','Infrarotsensor');
%    set(h,'location','northwest','FontSize',25)
%    xlabel('f [Hz]', 'FontSize', 30);
%    set(gca,'xscale','log', 'FontSize', 35)
%    ylabel('Phasenverschiebung[°]', 'FontSize', 35);
%    Ticks = -90:45:90;
%    set(gca, 'YTickMode', 'manual', 'YTick', Ticks, 'ylim', [-200,200]);
%    grid on; zoom on;
%     xlim([5e-2 60]);
%    ylim([-180 180]);
% 
% 
%    fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flugphase\Phase_' int2str(i) '\Spektrum_Phase\Spektrum_Hamming2048_SG_' int2str(nb) '_F' int2str(f-2) '.png'])
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flug\Flug_' int2str(f-2) '\Spektrum_Phase\Spektrum_Hamming2048_SG_' int2str(nb) '_Phase' int2str(i) '.png'])

%    close all;

end 


%% Spektrum Eingangsdaten (Druck)

for i = 1:7
 
    
    B.p_diff1 = B.p - B.MHP_p1(mb);
    B.p_diff2 = B.p - B.MHP_p2(mb);
    
    
   nx = max(size(Abschnitt(i).index));
   na = 2;                        %%%%%%%%%% Glättungsfaktor 1..20    
   ww = hamming(1024);

   
   [Pspektrum.P1, F1] = pwelch(B.p(Abschnitt(i).index), ww, 0, [], 100);
   [Pspektrum.P2, F2] = pwelch(B.p_diff1(Abschnitt(i).index), ww, 0, [], 100);
   [Pspektrum.P3, F3] = pwelch(B.p_diff2(Abschnitt(i).index), ww, 0, [], 100);
   
   
   freq = 0.01:0.01:50;
   f2 = 0.01:0.01:0.1;
   fl  = (freq) .^ (-5/3) .* 1e-3;  
   fk  = (freq) .^ (-5/3) .* 1e-5;  
   fm  = (f2) .^ (-4)   .* 1e-7;  


%    figure('visible','off','Name',strcat('Spektren Druck'));
%    loglog( F1 , Pspektrum.P1, 'color' ,farbe.dblau,'LineWidth',1);
%    hold on;
%    loglog( F2 , Pspektrum.P2, 'color' ,farbe.dgruen,'LineWidth',1);
%    hold on;
%    loglog( F3 , Pspektrum.P3, 'color' ,farbe.rot,'LineWidth',1);
%    hold on;
%    loglog( freq , fl, 'k--');
%    hold on;
%    loglog( freq , fk, 'k--');
%    hold on;
%    title(strcat({' Spektrum Druck Abschnitt-' Abschnitt(i).name 'Flug-' int2str(f-2)}), 'FontSize', 40); 
%    h=legend('P_0-Referenz-Messstelle','P_1-Differenz','P_2-Differenz');
%    set(h,'location','northeast','FontSize',25)
%    xlabel('Frequenz [Hz]', 'FontSize', 35);
%    ylabel('Druck [hPa]^2s', 'FontSize', 35);   
%    set(gca, 'FontSize', 30)
%    grid on; zoom on;
%    xlim([1e-3 60]);
%    ylim([1e-10 1e4]); 
% 
%    
%    
%    fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flugphase\Phase_' int2str(i) '\Spektrum_Druck\Spektrum_Hamming1024_F' int2str(f-2) '.png'])
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flug\Flug_' int2str(f-2) '\Spektrum_Druck\Spektrum_Hamming1024_Phase' int2str(i) '.png'])

% %    
%   close all;

end

%% Kohärenz- und Phasenspektrum (Druck)

for i = 3:4
   nx = max(size(Abschnitt(i).index));
   na = 2;                        %%%%%%%%%% Glättungsfaktor 1..20    
   ww = hamming(1024);
   
    if mod((length(Abschnitt(i).index)),(length(ww)))>= length(ww)/2
    fenster = floor(length(Abschnitt(i).index)/(length(ww)))+1;
    else
    fenster = floor(length(Abschnitt(i).index)/(length(ww)));
    end
    
   [P1coh, F1] = mscohere(B.p_diff1(Abschnitt(i).index),B.p(Abschnitt(i).index),  ww, 0, [], 100); 
   [P2coh, F2] = mscohere(B.p_diff2(Abschnitt(i).index),B.p(Abschnitt(i).index),  ww, 0, [], 100);

   nb = 3;
   
   P1cohs = sgolayfilt(P1coh,1,nb);
   P2cohs = sgolayfilt(P2coh,1,nb);
   
   freq1 = find(P1coh >= 0.5); 
   freq2 = find(P2coh >= 0.5);
   
   f50 = ones(length(F1), 1) .* 0.5;
   
%     figure('visible','off','Name',strcat('Spektren Kohärenz'));
%    semilogx(F1 , P1cohs,'color', farbe.dgruen,'LineWidth',1);hold on;
%    semilogx(F2 , P2cohs,'color', farbe.rot,'LineWidth',1);hold on;
%    semilogx(F1 , f50,'color', farbe.rot);
%    title(strcat( {[' Kohärenzspektrum 2048 DFT-Punkte - Auf P_0 '];[ Abschnitt(i).name, ' Flug ', int2str(f-2) ' ALICE Zarnekow, ' int2str(fenster) ' Fenster']}), 'FontSize', 40); 
%    h=legend('P_1-Differenz','P_2-Differenz');
%    set(h,'location','southwest','FontSize',25)
%    xlabel('f [Hz]', 'FontSize', 35); 
%    ylabel('Kohärenz', 'FontSize', 35); 
%    set(gca, 'FontSize', 30)
%    grid on; zoom on;
%    ylim([0 1]);
%    xlim([1e-2 60]);


  
%    fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flugphase\Phase_' int2str(i) '\Spektrum_Druck\KSpektrum_Hamming2048_SG_' int2str(nb) '_F' int2str(f-2) '.png'])
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flug\Flug_' int2str(f-2) '\Spektrum_Druck\KSpektrum_Hamming2048_SG_' int2str(nb) '_Phase' int2str(i) '.png'])



  clear nb;
   [Ph1, F1] = cpsd(B.p_diff1(Abschnitt(i).index),B.p(Abschnitt(i).index), ww, [], [], 100); 
   [Ph2, F2] = cpsd(B.p_diff2(Abschnitt(i).index),B.p(Abschnitt(i).index), ww, [], [], 100);


   nb = 5;
   
%    Ph1 = sgolayfilt(Ph1,1,nb);
%    Ph2 = sgolayfilt(Ph2,1,nb);

   
%    figure('visible','off','Name',strcat('Spektren Phase '));
%     scatter( F1 , angle(Ph1).*(180/pi),70.*P1coh+1,farbe.dgruen,'filled', 'MarkerFaceColor', farbe.dgruen);hold on;
%     scatter( F2 , angle(Ph2).*(180/pi),70.*P2coh+1,farbe.rot,'filled', 'MarkerFaceColor', farbe.rot);hold on;
%    title(strcat({[' Phasenspektrum - auf P_0 '];[ Abschnitt(i).name, ' Flug ', int2str(f-2)]}), 'FontSize', 40); 
%    h=legend('P_1-Differenz','P_2-Differenz');
%    set(h,'location','northwest','FontSize',25)
%    xlabel('f [Hz]', 'FontSize', 30);
%    set(gca,'xscale','log', 'FontSize', 35)
%    ylabel('Phasenverschiebung[°]', 'FontSize', 35);
%    Ticks = -90:45:90;
%    set(gca, 'YTickMode', 'manual', 'YTick', Ticks, 'ylim', [-200,200]);
%    grid on; zoom on;
%     xlim([5e-2 60]);
%    ylim([-180 180]);
% 
% 
%    fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flugphase\Phase_' int2str(i) '\Spektrum_Druck\PSpektrum_Hamming1024_SG_' int2str(nb) '_F' int2str(f-2) '.png'])
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flug\Flug_' int2str(f-2) '\Spektrum_Druck\PSpektrum_Hamming1024_SG_' int2str(nb) '_Phase' int2str(i) '.png'])

%    close all;

end 


%% Spektrum (Feuchte )

 % Umrechnung von RH in MV
  B.humi_E = (B.humi_rh(mb)./100) .* (6.107 + 30 * 0.000001 .* B.p) .*10.^(((8.082 -(B.humi(mb)./556)).*B.humi(mb))./(256.1+B.humi(mb)));
  B.humi_mv = 621.98 .* (B.humi_E./(B.p-B.humi_E));
%   
%   B.humi_Eroh = (B.humi_rhroh(2:end-1)./100) .* (6.107 + 30 * 0.000001 .* [B.p(1:100:end)]) .*10.^(((8.082 -(B.humiroh(2:end-1)./556)).*B.humiroh(2:end-1))./(256.1+B.humiroh(2:end-1)));
%   B.humi_mvroh = 621.98 .* (B.humi_Eroh./([B.p(1:100:end)]-B.humi_Eroh));
%   
  B.p14_E = (B.p14rh(mb)./100) .* (6.107 + 30 * 0.000001 .* B.p) .*10.^(((8.082 -(B.tsy(mb)./556)).*B.tsy(mb))./(256.1+B.tsy(mb)));
  B.p14_mv = 621.98 .* (B.p14_E./(B.p-B.p14_E));
  
  r.p14_mv = B.p14_mv - mean(B.p14_mv-B.humi_mv);
  r.humi_mv = B.humi_mv;
  
 for i = 3:4
     
   ww = hamming(2048);
   
    if mod((length(Abschnitt(i).index)),(length(ww)))>= length(ww)/2
    fenster = floor(length(Abschnitt(i).index)/(length(ww)))+1;
    else
    fenster = floor(length(Abschnitt(i).index)/(length(ww)));
    end
    
%     Abschnitth = [Abschnitt(i).index(1)+99:100:(Abschnitt(i).index(end))];
%     Abschnitth = floor(Abschnitth./100);
%    
%     r.humi_mvd(Abschnitt(i).index) =  detrend(r.humi_mv(Abschnitt(i).index), 'linear');
%     r.humi_mvrohd(Abschnitth) =  detrend(B.humi_mvroh(Abschnitth), 'linear');
%     r.p14_mvd(Abschnitt(i).index) =  detrend(r.p14_mv(Abschnitt(i).index), 'linear');
%     
%    
%    [MVspektrum.P1, F1] = pwelch(r.humi_mvrohd(Abschnitth), ww2, [], [], 1);
%    [MVspektrum.P2, F2] = pwelch(r.p14_mvd(Abschnitt(i).index), ww, [], [], 100);
% 
%    
%    freq = 0.01:0.01:50;
%    f2 = 0.01:0.01:0.1;
%    fl  = (freq) .^ (-5/3) .* 1e-3;  
%    fk  = (freq) .^ (-5/3) .* 1e-5;  
%    fm  = (f2) .^ (-4)   .* 1e-7;  
% 
% 
%    figure('visible','on','Name',strcat('Spektren Mischungsverhältnis'));
%    loglog( F1 , MVspektrum.P1, 'color' ,farbe.dblau,'LineWidth',2);
%    hold on;
%    loglog( F2 , MVspektrum.P2, 'color' ,farbe.rot,'LineWidth',2);
%    hold on;
%    loglog( freq , fl, 'k--');
%    hold on;
%    loglog( freq , fk, 'k--');
%    hold on; 
%    title(strcat({[ 'Spektrum Feuchte 2048 DFT-Punkte (64 Humicap)'];[ Abschnitt(i).name, ' Flug ', int2str(f-2) ' ALICE Zarnekow, ' int2str(fenster) ' Fenster']}), 'FontSize', 40); 
%    h=legend('Humicap','Rapid-P14');
%    set(h,'location','northeast','FontSize',25)
%    xlabel('Frequenz [Hz]', 'FontSize', 35);
%    ylabel('Mischungsverhältnis [g/kg]^2s', 'FontSize', 35);   
%    set(gca, 'FontSize', 30)
%    grid on; zoom on;
%    xlim([1e-3 60]);
%    ylim([1e-10 1e4]); 
%    
%    fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flugphase\Phase_' int2str(i) '\Spektrum_Feuchte\Spektrum_Hamming2048_F' int2str(f-2) '.png'])
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flug\Flug_' int2str(f-2) '\Spektrum_Feuchte\Spektrum_Hamming2048_Phase' int2str(i) '.png'])




     
  
  
 end
 
 
 
         %% Eingangsmessdaten   Feuchte MV über die Zeit

%    q=figure('visible','on','Name',strcat('Eingangsdaten'));
%    set(q, 'Units', 'normalized', 'Position',[0.1, 0.1, 0.8, 0.8] );
%    AxesH = axes('Units', 'normalized', 'Position', [0.082,0.135, 0.81, 0.74]);
%    plot(B.Bus_TOW,r.humi_mv, 'color' ,farbe.dgruen,'LineWidth',2);
%    hold on;
%    [hAx,hLine1,hLine2] = plotyy(B.Bus_TOW,r.p14_mv,B.Bus_TOW,B.h_a );
%    hold on;
%    hLine1.Color = farbe.dblau;
%    hLine2.Color = 'k';
%    hLine1.LineWidth = 2;
%    hLine2.LineWidth = 2;
%    set(hAx,'YColor','k', 'FontSize', 30);
%    h = legend( 'Humicap','Rapid P14','Höhe');
%    set(h,'FontSize',25)
%    set(gca, 'FontSize', 30)
%    xlabel('Sekunden des Tages [s]',  'FontSize', 30);
%    ylabel('Feuchte [g/kg]', 'FontSize', 40,'color' ,'k');
%    ylabel((hAx(2)),'Höhe [m]','color' ,'k', 'FontSize', 30);
%    set(hAx(2), 'FontSize', 30)
%    xx=fill(B.Bus_TOW(x3),y,'r',B.Bus_TOW(x4),y,'r');
%    set(xx,'facealpha',.08);
%    ylim([-2 max(r.p14_mv)+3]);
%    ylim((hAx(2)),[-60 max(B.h_a)*3]);
%    title(strcat({['Eingangsmessdaten Feuchte Flug ' int2str((f-2))] ; ['ALICE Zarnekow']})); 
%    grid on; zoom on;
%    
%      fig = gcf;
%      fig.PaperUnits = 'inches';
%      fig.PaperPosition = [0 0 22 13];
%      print(fig,'-dpng','-r0',['F:\ALICE_Zarnekow_2018-05-23\Plots\Eingangsmessdaten\Eingangsmessdaten_MV(t)_V_F' int2str(f-2) '.png'])

        %% Eingangsmessdaten   Höhe über die Feuchte MV
    
%    q=figure('visible','on','Name',strcat('Eingangsdaten'));
%    set(q, 'Units', 'normalized', 'Position',[0.1, 0.1, 0.8, 0.8] );
%    AxesH = axes('Units', 'normalized', 'Position', [0.1,0.135, 0.81, 0.74]);
%    plot(r.humi_mv,B.h_a, 'color' ,farbe.dgruen,'LineWidth',3);
%    hold on;
%    plot(r.p14_mv,B.h_a, 'color' ,farbe.dblau,'LineWidth',3);
%    set(gca, 'FontSize', 30)
%    title(strcat({['Feuchtegradient Flug ' int2str((f-2))] ; ['ALICE Zarnekow']})); 
%    h = legend( 'Humicap','Rapid P14');
%    set(h,'FontSize',25)
%    xlabel('Feuchte [g/kg]', 'FontSize', 30);
%    ylabel('Höhe [m]', 'FontSize', 30);
%    grid on; zoom on;
%    
%    
%      fig = gcf;
%      fig.PaperUnits = 'inches';
%      fig.PaperPosition = [0 0 22 13];
%      print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Eingangsmessdaten\Eingangsmessdaten_MV(h)_F' int2str(f-2) '.png'])


%    q=figure('visible','on','Name',strcat('Eingangsdaten'));
%    set(q, 'Units', 'normalized', 'Position',[0.1, 0.1, 0.8, 0.8] );
%    AxesH = axes('Units', 'normalized', 'Position', [0.082,0.135, 0.81, 0.74]);
%    subplot(1,2,1);
%    plot(r.fw0,B.h_a, 'color' ,farbe.dblau,'LineWidth', 2);
%    hold on;
%    plot(r.fw1,B.h_a, 'color' ,farbe.dgruen,'LineWidth', 2);
%    hold on;
%    plot(r.pt0,B.h_a, 'color' ,farbe.rot,'LineWidth', 2);
%    hold on;
%    plot(r.pt1,B.h_a, 'color' ,farbe.gelb,'LineWidth', 2);
%    hold on;
%    plot(r.tsy,B.h_a, 'color' ,farbe.lila,'LineWidth', 2);
%    hold on;
%    plot(r.humi,B.h_a, 'color' ,farbe.orange,'LineWidth', 2);
%    hold on;
%    plot(r.ir,B.h_a, 'color' ,farbe.hblau,'LineWidth', 2);
%    h = legend({ 'Finewire-0','Finewire-1',  [['Pt-100- ']; ['Sensor-0']]     ,  [['Pt-100- ']; ['Sensor-1']]  ,  'TSYS01'  ,  [['Humicap']; ['HMP-110']]  ,  'IR-Sensor'});
%    set(gca, 'FontSize', 30)
%    set(h,'FontSize',20)
%    xlabel('Temperatur [°C]', 'FontSize', 30);
%    ylabel('Höhe [m]', 'FontSize', 30);
% %    title(strcat({['Temperaturgradient'] ; ['Flug ' int2str(f-2)  ', ALICE Zarnekow']})); 
%    grid on; zoom on;
%    subplot(1,2,2);
%    plot(r.humi_mv,B.h_a, 'color' ,farbe.dgruen,'LineWidth',2);
%    hold on;
%    plot(r.p14_mv,B.h_a, 'color' ,farbe.dblau,'LineWidth',2);
%    set(gca, 'FontSize', 30)
% %    title(strcat({['Feuchtegradient'] ; ['Flug ' int2str(f-2)  ', ALICE Zarnekow']})); 
%    h = legend( 'Humicap','Rapid P14');
%    set(h,'FontSize',25)
%    xlabel('Feuchte [g/kg]', 'FontSize', 30);
%    ylabel('Höhe [m]', 'FontSize', 30);
%    grid on; zoom on;
%    
%         fig = gcf;
%      fig.PaperUnits = 'inches';
%      fig.PaperPosition = [0 0 22 13];
%      print(fig,'-dpng','-r0',['F:\ALICE_Zarnekow_2018-05-23\Plots\Eingangsmessdaten\Eingang_T(h)_MV(h)_oT_F' int2str(f-2) '.png'])




%% Kohärenz- und Phasenspektrum (Feuchte)

for i = 7:7
 
   ww = hamming(2048);
    
   [P1coh, F1] = mscohere(B.humi_mv(Abschnitt(i).index),B.p14_mv(Abschnitt(i).index),  ww, 0, [], 100); 
   
   nb = 3;
   
   P1cohs = sgolayfilt(P1coh,1,nb);

   
   freq1 = find(P1coh >= 0.5); 

   
   f50 = ones(length(F1), 1) .* 0.5;
   
%     figure('visible','on','Name',strcat('Spektren Kohärenz'));
%    semilogx(F1 , P1cohs,'color', farbe.dblau,'LineWidth',1);hold on;
%    semilogx(F1 , f50,'color', farbe.rot);
%    title(strcat( {[' Kohärenzspektrum - Auf P14 '];[ Abschnitt(i).name, ' Flug ', int2str(f-2)]}), 'FontSize', 40); 
%    h=legend('Humicap');
%    set(h,'location','southwest','FontSize',25)
%    xlabel('f [Hz]', 'FontSize', 35); 
%    ylabel('Kohärenz', 'FontSize', 35); 
%    set(gca, 'FontSize', 30)
%    grid on; zoom on;
%    ylim([0 1]);
%    xlim([1e-2 60]);
% 
% 
%   
%    fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flugphase\Phase_' int2str(i) '\Spektrum_Feuchte\KSpektrum_Hamming1024_SG_' int2str(nb) '_F' int2str(f-2) '.png'])
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flug\Flug_' int2str(f-2) '\Spektrum_Feuchte\KSpektrum_Hamming1024_SG_' int2str(nb) '_Phase' int2str(i) '.png'])



  clear nb;
   [Ph1, F1] = cpsd(B.humi_mv(Abschnitt(i).index),B.p14_mv(Abschnitt(i).index), ww, 0, [], 100); 

   nb = 3;
   
   Ph1 = sgolayfilt(Ph1,1,nb);


   
%    figure('visible','off','Name',strcat('Spektren Phase '));
%     scatter( F1 , angle(Ph1).*(180/pi),70.*P1coh+1,farbe.dblau,'filled', 'MarkerFaceColor', farbe.dgruen);hold on;
%    title(strcat({[' Phasenspektrum - auf P14 '];[ Abschnitt(i).name, ' Flug ', int2str(f-2)]}), 'FontSize', 40); 
%    h=legend('Humicap');
%    set(h,'location','northwest','FontSize',25)
%    xlabel('f [Hz]', 'FontSize', 30);
%    set(gca,'xscale','log', 'FontSize', 35)
%    ylabel('Phasenverschiebung[°]', 'FontSize', 35);
%    Ticks = -90:45:90;
%    set(gca, 'YTickMode', 'manual', 'YTick', Ticks, 'ylim', [-200,200]);
%    grid on; zoom on;
%     xlim([5e-2 60]);
%    ylim([-180 180]);
% 
% 
%    fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flugphase\Phase_' int2str(i) '\Spektrum_Feuchte\PSpektrum_Hamming1024_SG_' int2str(nb) '_F' int2str(f-2) '.png'])
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flug\Flug_' int2str(f-2) '\Spektrum_Feuchte\PSpektrum_Hamming1024_SG_' int2str(nb) '_Phase' int2str(i) '.png'])


end 

%% Korrelation
for i=1:7
    
    
for k=1:100

%  %Tiefpass
%  fs = 100;
%  N    = 1;
%  f_TP = (1/(90)+(k/100)-0.01);
%  [Bt, Wt] = butter(N, f_TP /(fs / 2),'low');
%  
%  tpk.fw0 =  filtfilt(Bt, Wt, r.fw0);
%  tpk.fw1 =  filtfilt(Bt, Wt, r.fw1);
%  tpk.pt0 =  filtfilt(Bt, Wt, r.pt0);
%  tpk.pt1 =  filtfilt(Bt, Wt, r.pt1);
%  tpk.tsy =  filtfilt(Bt, Wt, r.tsy); 
%  tpk.humi =  filtfilt(Bt, Wt, r.humi);
%  
%   %Hochpass
%  fs = 100;
%  N    = 1;
%  f_HP = (1/(100)+(k/100)-0.01);
%  [Bh, Wh] = butter(N, f_HP /(fs / 2),'high');
%  
%  hpk.fw0 =  filtfilt(Bh, Wh, r.fw0); 
%  hpk.fw1 =  filtfilt(Bh, Wh, r.fw1);
%  hpk.pt0 =  filtfilt(Bh, Wh, r.pt0);
%  hpk.pt1 =  filtfilt(Bh, Wh, r.pt1);
%  hpk.tsy =  filtfilt(Bh, Wh, r.tsy); 
%  hpk.humi =  filtfilt(Bh, Wh, r.humi);
%  
%  
%    %Bandpass
%  
%  bpk.fw0 =  filtfilt(Bt, Wt, hpk.fw0); 
%  bpk.fw1 =  filtfilt(Bt, Wt, hpk.fw1);
%  bpk.pt0 =  filtfilt(Bt, Wt, hpk.pt0);
%  bpk.pt1 =  filtfilt(Bt, Wt, hpk.pt1);
%  bpk.tsy =  filtfilt(Bt, Wt, hpk.tsy); 
%  bpk.humi =  filtfilt(Bt, Wt, hpk.humi);
%  
% 
% % [cor.acor1,cor.lag1]   = xcorr(bpk.tsy(Abschnitt(4).index),bpk.fw0(Abschnitt(4).index), 'coeff'); 
% % [cor.acor1a,cor.lag1a] = xcorr(bpk.humi(Abschnitt(4).index),bpk.fw0(Abschnitt(4).index), 'coeff');
% % [cor.acor1b,cor.lag1b] = xcorr(bpk.fw1(Abschnitt(4).index),bpk.fw0(Abschnitt(4).index), 'coeff');
% % [cor.acor1c,cor.lag1c] = xcorr(bpk.pt0(Abschnitt(4).index),bpk.fw0(Abschnitt(4).index), 'coeff');
% % [cor.acor1d,cor.lag1d] = xcorr(bpk.pt1(Abschnitt(4).index),bpk.fw0(Abschnitt(4).index), 'coeff');
% % 
% % 
% % 
% % cor.max(k) = max(cor.acor1);
% % cor.maxa(k) = max(cor.acor1a);
% % cor.maxb(k) = max(cor.acor1b);
% % cor.maxc(k) = max(cor.acor1c);
% % cor.maxd(k) = max(cor.acor1d);
% % cor.freq(k)= ((1/(1000)+(k/1000)-0.001)+(1/(900)+(k/1000)-0.001))/2;

% a1 = corrcoef(bpk.tsy(Abschnitt(i).index),bpk.fw0(Abschnitt(i).index));
% a2 = corrcoef(bpk.fw1(Abschnitt(i).index),bpk.fw0(Abschnitt(i).index));
% a3 = corrcoef(bpk.pt1(Abschnitt(i).index),bpk.fw0(Abschnitt(i).index));
% a4 = corrcoef(bpk.pt0(Abschnitt(i).index),bpk.fw0(Abschnitt(i).index));
% a5 = corrcoef(bpk.humi(Abschnitt(i).index),bpk.fw0(Abschnitt(i).index));
% 
% 
% cor.coef1(k) = a1(2:2);
% cor.coef2(k) = a2(2:2);
% cor.coef3(k) = a3(2:2);
% cor.coef4(k) = a4(2:2);
% cor.coef5(k) = a5(2:2);
% 
% % cor.max(k) = max(cor.acor1);
% % cor.maxa(k) = max(cor.acor1a);
% cor.freq(k)= ((1/(1000)+(k/1000)-0.001)+(1/(900)+(k/1000)-0.001))/2;



%        figure('Name', 'Kreuzkorrelation Flug')
%         plot(cor.lag1./20,cor.acor1, 'color' ,farbe.dgruen,'LineWidth',3);
%         hold on;
%         plot(cor.lag1a./20,cor.acor1a, 'color' ,farbe.rot,'LineWidth',3);
%         hold on;
%         title(strcat({' Kreuzkorrelation'}), 'FontSize', 40); 
%         h=legend('TSY/FW0','Humi/FW0');
%         set(h,'FontSize',20)
%         set(gca,'FontSize',30)
%         xlabel('Zeit [s]', 'FontSize', 40);
%         ylabel('Korrelation', 'FontSize', 40);
%         grid on; zoom on;
%         ylim([-0.3 1]);
%         
% 
%     fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['G:\\Plots\Korrelation\Korr-' int2str(100-k*20) '-' int2str(195-k*20) 'Hz.png'])

end
 


%         figure('visible','off','Name', 'Korrelationsabfall')
%         semilogx(cor.freq,cor.coef1, 'color' ,farbe.dgruen,'LineWidth',3);
%         hold on;
%         semilogx(cor.freq,cor.coef2, 'color' ,farbe.dblau,'LineWidth',3);
%         hold on;
%         semilogx(cor.freq,cor.coef3, 'color' ,farbe.rot,'LineWidth',3);
%         hold on;
%         semilogx(cor.freq,cor.coef4, 'color' ,farbe.orange,'LineWidth',3);
%         hold on;
%         semilogx(cor.freq,cor.coef5, 'color' ,farbe.hblau,'LineWidth',3);
%         hold on;
%         title(strcat({['Maximale Kreuzkorrelation ' Abschnitt(i).name];[ 'Flug ' num2str(f-2)]}), 'FontSize', 40); 
%         h=legend('TSY/FW0','FW1/FW0','PT1/FW0','PT0/FW0','Humi/FW0');
%         set(h,'FontSize',20)
%         set(gca,'FontSize',30)
%         xlabel('Frequenz [Hz]', 'FontSize', 40);
%         ylabel('Korrelation', 'FontSize', 40);
%         grid on; zoom on;
%         ylim([-0.1 1]);
%         
%         
%            fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Plots\Daten_nach_Flugphase\Phase_' int2str(i) '\Korrelation\Korrelation_F' int2str(f-2) '.png'])
%      
        
end
clear f_TP f_HP




 
 
 
end 
end

clearvars -except farbe Daten r Abschnitt Abschnitte mb B f50 Abschnitth 

Anz.ges=0;Anz.Abschnitt=0;

for i=3:4
    
   Tspektrum.ges1 = 0; Tspektrum.ges2 = 0;Tspektrum.ges3 = 0;Tspektrum.ges4 = 0;Tspektrum.ges5 = 0;Tspektrum.ges6 = 0;Tspektrum.ges7 = 0;
   MVspektrum.ges1 = 0; MVspektrum.ges2 = 0;
   P1coh.ges1 = 0;P2coh.ges2 = 0;P3coh.ges3 = 0;P4coh.ges4 = 0;P5coh.ges5 = 0;P6coh.ges6 = 0;P1mvcoh.ges1 = 0;
   std.ges1 = 0;std.ges2 = 0;std.ges3 = 0;std.ges4 = 0;std.ges5 = 0;std.ges6 = 0;std.ges7 = 0;
   stdmv.ges1 = 0;stdmv.ges2 = 0;
   std.gescoh1 = 0;std.gescoh2 = 0;std.gescoh3 = 0;std.gescoh4 = 0;std.gescoh5 = 0;std.gescoh6 = 0;std.gescoh7 = 0;stdmv.gescoh1 = 0;
   
   
   
   
   
   for f=3:length(Daten) % Mittelwert für Spektren
       
    cd('Korrigierte Daten');
    load(Daten(f).name); 
    cd('..');
    cd('Abschnitte');
    load(Abschnitte(f).name); % Abschnitte laden
    cd('..');
    
    mb =(100: length(B.Bus_TOW)); % ersten 100 Werte abgeschnitten
    B.Bus_TOW = B.Bus_TOW(mb);
    B.h_a= B.h_a(mb);
    B.p = B.p(mb);
    bereich = [mb];
   

 r.tsy = B.tsy(bereich);
 r.fw0 = B.fw0(bereich) - mean(B.fw0(bereich)-B.tsy(bereich));
 r.fw1 = B.fw1(bereich) - mean(B.fw1(bereich)-B.tsy(bereich));
 r.pt0 = B.pt0(bereich) - mean(B.pt0(bereich)-B.tsy(bereich));
 r.pt1 = B.pt1(bereich) - mean(B.pt1(bereich)-B.tsy(bereich));
 r.humi = B.humi(bereich) - mean(B.humi(bereich)-B.tsy(bereich));
 r.ir = B.ir(bereich) - mean(B.ir(bereich)-B.tsy(bereich));
 
 B.humi_E = (B.humi_rh(mb)./100) .* (6.107 + 30 * 0.000001 .* B.p) .*10.^(((8.082 -(B.humi(mb)./556)).*B.humi(mb))./(256.1+B.humi(mb)));
  B.humi_mv = 621.98 .* (B.humi_E./(B.p-B.humi_E));
    
  u = length(B.humiroh)-length([B.p(1:100:end)]);
  
  B.humi_Eroh = (B.humi_rhroh(1:end-u)./100) .* (6.107 + 30 * 0.000001 .* [B.p(1:100:end)]) .*10.^(((8.082 -(B.humiroh(1:end-u)./556)).*B.humiroh(1:end-u))./(256.1+B.humiroh(1:end-u)));
  B.humi_mvroh = 621.98 .* (B.humi_Eroh./([B.p(1:100:end)]-B.humi_Eroh));
  B.p14_E = (B.p14rh(mb)./100) .* (6.107 + 30 * 0.000001 .* B.p) .*10.^(((8.082 -(B.tsy(mb)./556)).*B.tsy(mb))./(256.1+B.tsy(mb)));
  B.p14_mv = 621.98 .* (B.p14_E./(B.p-B.p14_E));
  
  r.p14_mv = B.p14_mv - mean(B.p14_mv-B.humi_mv);
  r.humi_mv = B.humi_mv;
  r.humi_mvroh = B.humi_mvroh;

   Abschnitth = [Abschnitt(i).index(1)+99:100:Abschnitt(i).index(end)];
   Abschnitth = floor(Abschnitth./100);
   
 r.fw0(Abschnitt(i).index) =  detrend(r.fw0(Abschnitt(i).index), 'linear');
 r.fw1(Abschnitt(i).index) =  detrend(r.fw1(Abschnitt(i).index), 'linear');
 r.pt0(Abschnitt(i).index) =  detrend(r.pt0(Abschnitt(i).index), 'linear'); 
 r.pt1(Abschnitt(i).index) =  detrend(r.pt1(Abschnitt(i).index), 'linear');
 r.tsy(Abschnitt(i).index) =  detrend(r.tsy(Abschnitt(i).index), 'linear'); 
 r.humi(Abschnitt(i).index) =  detrend(r.humi(Abschnitt(i).index), 'linear'); 
 r.ir(Abschnitt(i).index) =  detrend(r.ir(Abschnitt(i).index), 'linear');
 r.humiroh(Abschnitth) =  detrend(B.humiroh(Abschnitth), 'linear'); 
 
 r.humi_mvrohd(Abschnitth) =  detrend(r.humi_mvroh(Abschnitth), 'linear');
 r.p14_mvd(Abschnitt(i).index) =  detrend(r.p14_mv(Abschnitt(i).index), 'linear'); 

    ww = hamming(2048);
  
    
    ww2 = hamming(64);
   
    if mod((length(Abschnitt(i).index)),(length(ww)))>= length(ww)/2
    Anz.Abschnitt = floor(length(Abschnitt(i).index)/(length(ww)))+1;
    else
    Anz.Abschnitt = floor(length(Abschnitt(i).index)/(length(ww)));
    end
 
   Anz.ges = Anz.ges + Anz.Abschnitt ;
   
   
 
   
   [Tspektrum.P1, F1] = pwelch(r.fw0(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P2, F2] = pwelch(r.fw1(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P3, F3] = pwelch(r.pt0(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P4, F4] = pwelch(r.pt1(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P5, F5] = pwelch(r.tsy(Abschnitt(i).index), ww, [], [], 100);
   %[Tspektrum.P6, F6] = pwelch(r.humi(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P6, F6] = pwelch(r.humiroh(Abschnitth), ww2, [], [], 1);
   [Tspektrum.P7, F7] = pwelch(r.ir(Abschnitt(i).index), ww, [], [], 100);
   
   
   Tspektrum.ges1 = (Tspektrum.ges1 .* (f-3) + Tspektrum.P1) ./ (f-2) ;
   Tspektrum.ges2 = (Tspektrum.ges2 .* (f-3) + Tspektrum.P2) ./ (f-2) ;
   Tspektrum.ges3 = (Tspektrum.ges3 .* (f-3) + Tspektrum.P3) ./ (f-2) ;
   Tspektrum.ges4 = (Tspektrum.ges4 .* (f-3) + Tspektrum.P4) ./ (f-2) ;
   Tspektrum.ges5 = (Tspektrum.ges5 .* (f-3) + Tspektrum.P5) ./ (f-2) ;
   Tspektrum.ges6 = (Tspektrum.ges6 .* (f-3) + Tspektrum.P6) ./ (f-2) ;
   Tspektrum.ges7 = (Tspektrum.ges7 .* (f-3) + Tspektrum.P7) ./ (f-2) ;
   
   [P1coh.a, F1] = mscohere(r.fw1(Abschnitt(i).index),r.fw0(Abschnitt(i).index),  ww, 0, [], 100); 
   [P2coh.a, F2] = mscohere(r.pt0(Abschnitt(i).index),r.fw0(Abschnitt(i).index),  ww, 0, [], 100);
   [P3coh.a, F3] = mscohere(r.pt1(Abschnitt(i).index),r.fw0(Abschnitt(i).index),  ww, 0, [], 100);
   [P4coh.a, F4] = mscohere(r.tsy(Abschnitt(i).index),r.fw0(Abschnitt(i).index),  ww, 0, [], 100);
   [P5coh.a, F5] = mscohere(r.humi(Abschnitt(i).index),r.fw0(Abschnitt(i).index), ww, 0, [], 100);
   [P6coh.a, F6] = mscohere(r.ir(Abschnitt(i).index),r.fw0(Abschnitt(i).index),   ww, 0, [], 100);
   
   
   P1coh.ges1 = (P1coh.ges1 .* (f-3) + P1coh.a) ./ (f-2) ;
   P2coh.ges2 = (P2coh.ges2 .* (f-3) + P2coh.a) ./ (f-2) ;
   P3coh.ges3 = (P3coh.ges3 .* (f-3) + P3coh.a) ./ (f-2) ;
   P4coh.ges4 = (P4coh.ges4 .* (f-3) + P4coh.a) ./ (f-2) ;
   P5coh.ges5 = (P5coh.ges5 .* (f-3) + P5coh.a) ./ (f-2) ;
   P6coh.ges6 = (P6coh.ges6 .* (f-3) + P6coh.a) ./ (f-2) ;
   
   [MVspektrum.P1, F1] = pwelch(r.humi_mvrohd(Abschnitth), ww2, [], [], 1);
   [MVspektrum.P2, F2] = pwelch(r.p14_mvd(Abschnitt(i).index), ww, [], [], 100);

   MVspektrum.ges1 = (MVspektrum.ges1 .* (f-3) + MVspektrum.P1) ./ (f-2) ;
   MVspektrum.ges2 = (MVspektrum.ges2 .* (f-3) + MVspektrum.P2) ./ (f-2) ;
   
   [P1mvcoh.a, F1] = mscohere(r.humi_mv(Abschnitt(i).index),r.p14_mv(Abschnitt(i).index),  ww, 0, [], 100); 
   
   P1mvcoh.ges1 = (P1mvcoh.ges1 .* (f-3) + P1mvcoh.a) ./ (f-2) ;
    
   end
   
   
   for f=3:length(Daten) % Standardabweichung für Spektren
      
           cd('Korrigierte Daten');
    load(Daten(f).name); 
    cd('..');
    cd('Abschnitte');
    load(Abschnitte(f).name); % Abschnitte laden
    cd('..');
    
    mb =(100: length(B.Bus_TOW)); % ersten 100 Werte abgeschnitten
    B.Bus_TOW = B.Bus_TOW(mb);
    B.h_a= B.h_a(mb);
    B.p = B.p(mb);
    bereich = [mb];


 r.tsy = B.tsy(bereich);
 r.fw0 = B.fw0(bereich) - mean(B.fw0(bereich)-B.tsy(bereich));
 r.fw1 = B.fw1(bereich) - mean(B.fw1(bereich)-B.tsy(bereich));
 r.pt0 = B.pt0(bereich) - mean(B.pt0(bereich)-B.tsy(bereich));
 r.pt1 = B.pt1(bereich) - mean(B.pt1(bereich)-B.tsy(bereich));
 r.humi = B.humi(bereich) - mean(B.humi(bereich)-B.tsy(bereich));
 r.ir = B.ir(bereich) - mean(B.ir(bereich)-B.tsy(bereich));
 
  B.humi_E = (B.humi_rh(mb)./100) .* (6.107 + 30 * 0.000001 .* B.p) .*10.^(((8.082 -(B.humi(mb)./556)).*B.humi(mb))./(256.1+B.humi(mb)));
  B.humi_mv = 621.98 .* (B.humi_E./(B.p-B.humi_E));
    
  u = length(B.humiroh)-length([B.p(1:100:end)]);
  
  B.humi_Eroh = (B.humi_rhroh(1:end-u)./100) .* (6.107 + 30 * 0.000001 .* [B.p(1:100:end)]) .*10.^(((8.082 -(B.humiroh(1:end-u)./556)).*B.humiroh(1:end-u))./(256.1+B.humiroh(1:end-u)));
  B.humi_mvroh = 621.98 .* (B.humi_Eroh./([B.p(1:100:end)]-B.humi_Eroh));
  B.p14_E = (B.p14rh(mb)./100) .* (6.107 + 30 * 0.000001 .* B.p) .*10.^(((8.082 -(B.tsy(mb)./556)).*B.tsy(mb))./(256.1+B.tsy(mb)));
  B.p14_mv = 621.98 .* (B.p14_E./(B.p-B.p14_E));
  
  r.p14_mv = B.p14_mv - mean(B.p14_mv-B.humi_mv);
  r.humi_mv = B.humi_mv;
  r.humi_mvroh = B.humi_mvroh;
 
 r.fw0(Abschnitt(i).index) =  detrend(r.fw0(Abschnitt(i).index), 'linear');
 r.fw1(Abschnitt(i).index) =  detrend(r.fw1(Abschnitt(i).index), 'linear');
 r.pt0(Abschnitt(i).index) =  detrend(r.pt0(Abschnitt(i).index), 'linear'); 
 r.pt1(Abschnitt(i).index) =  detrend(r.pt1(Abschnitt(i).index), 'linear');
 r.tsy(Abschnitt(i).index) =  detrend(r.tsy(Abschnitt(i).index), 'linear'); 
 r.humi(Abschnitt(i).index) =  detrend(r.humi(Abschnitt(i).index), 'linear'); 
 r.ir(Abschnitt(i).index) =  detrend(r.ir(Abschnitt(i).index), 'linear');
 r.humiroh(Abschnitth) =  detrend(B.humiroh(Abschnitth), 'linear'); 
 
 r.humi_mvrohd(Abschnitth) =  detrend(r.humi_mvroh(Abschnitth), 'linear');
 r.p14_mvd(Abschnitt(i).index) =  detrend(r.p14_mv(Abschnitt(i).index), 'linear'); 
     
       
   [Tspektrum.P1, F1] = pwelch(r.fw0(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P2, F2] = pwelch(r.fw1(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P3, F3] = pwelch(r.pt0(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P4, F4] = pwelch(r.pt1(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P5, F5] = pwelch(r.tsy(Abschnitt(i).index), ww, [], [], 100);
   %[Tspektrum.P6, F6h] = pwelch(r.humi(Abschnitt(i).index), ww, [], [], 100);
   [Tspektrum.P6, F6h] = pwelch(r.humiroh(Abschnitth), ww2, [], [], 1);
   [Tspektrum.P7, F7] = pwelch(r.ir(Abschnitt(i).index), ww, [], [], 100);
   
   [P1coh.a, F1] = mscohere(r.fw1(Abschnitt(i).index),r.fw0(Abschnitt(i).index),  ww, 0, [], 100); 
   [P2coh.a, F2] = mscohere(r.pt0(Abschnitt(i).index),r.fw0(Abschnitt(i).index),  ww, 0, [], 100);
   [P3coh.a, F3] = mscohere(r.pt1(Abschnitt(i).index),r.fw0(Abschnitt(i).index),  ww, 0, [], 100);
   [P4coh.a, F4] = mscohere(r.tsy(Abschnitt(i).index),r.fw0(Abschnitt(i).index),  ww, 0, [], 100);
   [P5coh.a, F5] = mscohere(r.humi(Abschnitt(i).index),r.fw0(Abschnitt(i).index), ww, 0, [], 100);
   [P6coh.a, F6] = mscohere(r.ir(Abschnitt(i).index),r.fw0(Abschnitt(i).index),   ww, 0, [], 100);
   
   [MVspektrum.P1, F1] = pwelch(r.humi_mvrohd(Abschnitth), ww2, [], [], 1);
   [MVspektrum.P2, F2] = pwelch(r.p14_mvd(Abschnitt(i).index), ww, [], [], 100);
   
   [P1mvcoh.a, F1] = mscohere(r.humi_mv(Abschnitt(i).index),r.p14_mv(Abschnitt(i).index),  ww, 0, [], 100); 
       
    std.ges1 = std.ges1 + (Tspektrum.P1 - Tspektrum.ges1).^2 ;
    std.ges2 = std.ges2 + (Tspektrum.P2 - Tspektrum.ges2).^2 ;
    std.ges3 = std.ges3 + (Tspektrum.P3 - Tspektrum.ges3).^2 ;
    std.ges4 = std.ges4 + (Tspektrum.P4 - Tspektrum.ges4).^2 ;
    std.ges5 = std.ges5 + (Tspektrum.P5 - Tspektrum.ges5).^2 ;
    std.ges6 = std.ges6 + (Tspektrum.P6 - Tspektrum.ges6).^2 ;
    std.ges7 = std.ges7 + (Tspektrum.P7 - Tspektrum.ges7).^2 ;
    
    std.gescoh1 = std.gescoh1 + (P1coh.a - P1coh.ges1).^2 ;
    std.gescoh2 = std.gescoh2 + (P2coh.a - P2coh.ges2).^2 ;
    std.gescoh3 = std.gescoh3 + (P3coh.a - P3coh.ges3).^2 ;
    std.gescoh4 = std.gescoh4 + (P4coh.a - P4coh.ges4).^2 ;
    std.gescoh5 = std.gescoh5 + (P5coh.a - P5coh.ges5).^2 ;
    std.gescoh6 = std.gescoh6 + (P6coh.a - P6coh.ges6).^2 ;
    
    stdmv.ges1 = stdmv.ges1 + (MVspektrum.P1 - MVspektrum.ges1).^2 ;
    stdmv.ges2 = stdmv.ges2 + (MVspektrum.P2 - MVspektrum.ges2).^2 ;
    
    stdmv.gescoh1 = stdmv.gescoh1 + (P1mvcoh.a - P1mvcoh.ges1).^2 ;

 
   
    
   end
   
   std.ges1 = (1/(length(Daten)-3) .* std.ges1).^0.5;
   std.ges2 = (1/(length(Daten)-3) .* std.ges2).^0.5;
   std.ges3 = (1/(length(Daten)-3) .* std.ges3).^0.5;
   std.ges4 = (1/(length(Daten)-3) .* std.ges4).^0.5;
   std.ges5 = (1/(length(Daten)-3) .* std.ges5).^0.5;
   std.ges6 = (1/(length(Daten)-3) .* std.ges6).^0.5;
   std.ges7 = (1/(length(Daten)-3) .* std.ges7).^0.5;
   
   std.gescoh1 = (1/(length(Daten)-3) .* std.gescoh1).^0.5;
   std.gescoh2 = (1/(length(Daten)-3) .* std.gescoh2).^0.5;
   std.gescoh3 = (1/(length(Daten)-3) .* std.gescoh3).^0.5;
   std.gescoh4 = (1/(length(Daten)-3) .* std.gescoh4).^0.5;
   std.gescoh5 = (1/(length(Daten)-3) .* std.gescoh5).^0.5;
   std.gecohs6 = (1/(length(Daten)-3) .* std.gescoh6).^0.5;

   stdmv.ges1 = (1/(length(Daten)-3) .* stdmv.ges1).^0.5;
   stdmv.ges2 = (1/(length(Daten)-3) .* stdmv.ges2).^0.5;
   
   stdmv.gescoh1 = (1/(length(Daten)-3) .* stdmv.gescoh1).^0.5;
   

   nb = 3;
   
   P1cohs = sgolayfilt(P1coh.ges1,1,nb);
   P2cohs = sgolayfilt(P2coh.ges2,1,nb);
   P3cohs = sgolayfilt(P3coh.ges3,1,nb);
   P4cohs = sgolayfilt(P4coh.ges4,1,nb);
   P5cohs = sgolayfilt(P5coh.ges5,1,nb);
   P6cohs = sgolayfilt(P6coh.ges6,1,nb);
   
   P1mvcohs = sgolayfilt(P1mvcoh.ges1,1,nb);
   
   %% plot
for k=1
   X.T1=[F1(6:end);flipud(F1(6:end))];               
   Y.T1=[Tspektrum.ges1(6:end)+std.ges1(6:end);flipud(Tspektrum.ges1(6:end)-std.ges1(6:end))];     
   X.T2=[F2(6:end);flipud(F2(6:end))];               
   Y.T2=[Tspektrum.ges2(6:end)+std.ges2(6:end);flipud(Tspektrum.ges2(6:end)-std.ges2(6:end))]; 
   index = (Y.T2<0);
   Y.T2(index) = 10^-14;
   X.T3=[F3(6:end);flipud(F3(6:end))];               
   Y.T3=[Tspektrum.ges3(6:end)+std.ges3(6:end);flipud(Tspektrum.ges3(6:end)-std.ges3(6:end))]; 
   index = (Y.T3<0);
   Y.T3(index) = 10^-14;
   X.T4=[F4(6:end);flipud(F4(6:end))];               
   Y.T4=[Tspektrum.ges4(6:end)+std.ges4(6:end);flipud(Tspektrum.ges4(6:end)-std.ges4(6:end))];
   index = (Y.T4<0);
   Y.T4(index) = 10^-14;
   X.T5=[F5(6:end);flipud(F5(6:end))];               
   Y.T5=[Tspektrum.ges5(6:end)+std.ges5(6:end);flipud(Tspektrum.ges5(6:end)-std.ges5(6:end))]; 
   X.T6=[F6h(6:end);flipud(F6h(6:end))];               
   Y.T6=[Tspektrum.ges6(6:end)+std.ges6(6:end);flipud(Tspektrum.ges6(6:end)-std.ges6(6:end))]; 
   index = (Y.T6<0);
   Y.T6(index) = 10^-14;
   X.T7=[F7(6:end);flipud(F7(6:end))];               
   Y.T7=[Tspektrum.ges7(6:end)+std.ges7(6:end);flipud(Tspektrum.ges7(6:end)-std.ges7(6:end))]; 
   
   X.coh1=[F1(6:end);flipud(F1(6:end))];               
   Y.coh1=[P1coh.ges1(6:end)+std.gescoh1(6:end);flipud(P1coh.ges1(6:end)-std.gescoh1(6:end))];  
   index = (Y.coh1<0);
   Y.coh1(index) = 10^-14;
   X.coh2=[F2(6:end);flipud(F2(6:end))];               
   Y.coh2=[P2coh.ges2(6:end)+std.gescoh2(6:end);flipud(P2coh.ges2(6:end)-std.gescoh2(6:end))]; 
   index = (Y.coh2<0);
   Y.coh2(index) = 10^-14;
   X.coh3=[F3(6:end);flipud(F3(6:end))];               
   Y.coh3=[P3coh.ges3(6:end)+std.gescoh3(6:end);flipud(P3coh.ges3(6:end)-std.gescoh3(6:end))]; 
   index = (Y.coh3<0);
   Y.coh3(index) = 10^-14;
   X.coh4=[F4(6:end);flipud(F4(6:end))];               
   Y.coh4=[P4coh.ges4(6:end)+std.gescoh4(6:end);flipud(P4coh.ges4(6:end)-std.gescoh4(6:end))];
   index = (Y.coh4<0);
   Y.coh4(index) = 10^-14;
   X.coh5=[F5(6:end);flipud(F5(6:end))];               
   Y.coh5=[P5coh.ges5(6:end)+std.gescoh5(6:end);flipud(P5coh.ges5(6:end)-std.gescoh5(6:end))]; 
   index = (Y.coh5<0);
   Y.coh5(index) = 10^-14;
   X.coh6=[F6(6:end);flipud(F6(6:end))];               
   Y.coh6=[P6coh.ges6(6:end)+std.gescoh6(6:end);flipud(P6coh.ges6(6:end)-std.gescoh6(6:end))];
   index = (Y.coh6<0);
   Y.coh6(index) = 10^-14;
   
   Xmv.T1=[F6h(6:end);flipud(F6h(6:end))];               
   Ymv.T1=[MVspektrum.ges1(6:end)+stdmv.ges1(6:end);flipud(MVspektrum.ges1(6:end)-stdmv.ges1(6:end))];   
   index = (Ymv.T1<0);
   Ymv.T1(index) = 10^-14;
   Xmv.T2=[F2(6:end);flipud(F2(6:end))];               
   Ymv.T2=[MVspektrum.ges2(6:end)+stdmv.ges2(6:end);flipud(MVspektrum.ges2(6:end)-stdmv.ges2(6:end))]; 
   index = (Ymv.T2<0);
   Ymv.T2(index) = 10^-14;
   
   Xmv.coh1=[F1(6:end);flipud(F1(6:end))];               
   Ymv.coh1=[P1mvcoh.ges1(6:end)+stdmv.gescoh1(6:end);flipud(P1mvcoh.ges1(6:end)-stdmv.gescoh1(6:end))];  
   index = (Ymv.coh1<0);
   Ymv.coh1(index) = 10^-14;
end  %Fillgrenzen

   
   freq = 0.01:0.01:50;
   f2 = 0.01:0.01:0.1;
   fl  = (freq) .^ (-5/3) .* 1e-3;  
   fk  = (freq) .^ (-5/3) .* 1e-5;  
   fm  = (f2) .^ (-4)   .* 1e-7;  
   f50 = ones(length(F1), 1) .* 0.5;
   
   
   
   q=figure('visible','on','Name',strcat('Spektren Mischungsverhältnis'));
   set(q, 'Units', 'normalized', 'Position',[0.1, 0.1, 0.8, 0.8] );
   AxesH = axes('Units', 'normalized', 'Position', [0.12,0.145, 0.85, 0.74]);
   loglog( F1 , Tspektrum.ges1, 'color' ,farbe.dblau,'LineWidth',2);
   hold on;
   loglog( F2 , Tspektrum.ges2, 'color' ,farbe.dgruen,'LineWidth',2);
   hold on;
   loglog( F3 , Tspektrum.ges3, 'color' ,farbe.rot,'LineWidth',2);
   hold on;
   loglog( F4 , Tspektrum.ges4, 'color' ,farbe.gelb,'LineWidth',2);
   hold on;
   loglog( F5 , Tspektrum.ges5, 'color' ,farbe.lila,'LineWidth',2);
   hold on;
   loglog( F6h , Tspektrum.ges6, 'color' ,farbe.orange,'LineWidth',2);
   hold on;
   loglog( F7 , Tspektrum.ges7, 'color' ,farbe.hblau,'LineWidth',2);
   hold on;
   loglog( freq , fl, 'k--');
   hold on;
   loglog( freq , fk, 'k--');
   hold on;
   title(strcat({'ALICE Steigflüge'}), 'FontSize', 40); 
%    title(strcat({[' Leistungsdichtespektrum Temperatur ' ];[ 'Mittelung über ' int2str(Anz.ges) ' Fenster ' Abschnitt(i).name ' - ALICE Zarnekow']}), 'FontSize', 40); 
   h=legend( 'Finewire-0','Finewire-1','Pt-100-Sensor-0','Pt-100-Sensor-1','TSYS01','Humicap HMP110','Infrarotsensor');
   set(h,'location','northeast','FontSize',25)
   xlabel('Frequenz [Hz]', 'FontSize', 35);
   ylabel('Temperatur [°C]^2s', 'FontSize', 35);
%    xx=fill(X.T1,Y.T1,farbe.dblau,'LineStyle','none');  
%    set(xx,'facealpha',.2);
%    xx=fill(X.T2,Y.T2,farbe.dgruen,'LineStyle','none');  
%    set(xx,'facealpha',.2);
%    xx=fill(X.T3,Y.T3,farbe.rot,'LineStyle','none');  
%    set(xx,'facealpha',.2);
%    xx=fill(X.T4,Y.T4,farbe.gelb,'LineStyle','none');  
%    set(xx,'facealpha',.3);
%    xx=fill(X.T5,Y.T5,farbe.lila,'LineStyle','none');  
%    set(xx,'facealpha',.2);
%    xx=fill(X.T6,Y.T6,farbe.orange,'LineStyle','none');  
%    set(xx,'facealpha',.2);
%    xx=fill(X.T7,Y.T7,farbe.hblau,'LineStyle','none');  
%    set(xx,'facealpha',.3);
   set(gca, 'FontSize', 30)
   grid on; zoom on;
   xlim([1e-3 60]);
   ylim([1e-10 1e4]); 
   
%    fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['F:\ALICE_Zarnekow_2018-05-23\Mittelung\', Abschnitt(i).name, 'e_V' '.png'])

   
   
   q=figure('visible','on','Name',strcat('Spektren Kohärenz'));
   set(q, 'Units', 'normalized', 'Position',[0.1, 0.1, 0.8, 0.8] );
   AxesH = axes('Units', 'normalized', 'Position', [0.1,0.135, 0.87, 0.74]);
   semilogx(F1 , P1cohs,'color', farbe.dgruen,'LineWidth',2);hold on;
   semilogx(F2 , P2cohs,'color', farbe.rot,'LineWidth',2);hold on;
   semilogx(F3 , P3cohs,'color', farbe.gelb,'LineWidth',2);hold on;
   semilogx(F4 , P4cohs,'color', farbe.lila,'LineWidth',2);hold on;
   semilogx(F5 , P5cohs,'color', farbe.orange,'LineWidth',2);hold on;
   semilogx(F6 , P6cohs,'color', farbe.hblau,'LineWidth',2);hold on;
   semilogx(F1 , f50,'color', farbe.rot);
   title(strcat({'ALICE Steigflüge'}), 'FontSize', 40);
%    title(strcat( {[' Kohärenzspektrum Temperatur - Auf FW0 '];[' Mittelung über ' int2str(Anz.ges) ' Fenster ' Abschnitt(i).name ' - ALICE Zarnekow']}), 'FontSize', 40); 
   h=legend('Finewire-1','Pt-100-Sensor-0','Pt-100-Sensor-1','TSYS01','Humicap HMP110','Infrarotsensor');
   set(h,'location','southwest','FontSize',25)
   xlabel('Frequenz [Hz]', 'FontSize', 35); 
   ylabel('Kohärenz', 'FontSize', 35); 
%    xx=fill(X.coh1,Y.coh1,farbe.dgruen,'LineStyle','none');  
%    set(xx,'facealpha',.2);
%    xx=fill(X.coh2,Y.coh2,farbe.rot,'LineStyle','none');  
%    set(xx,'facealpha',.2);
%    xx=fill(X.coh3,Y.coh3,farbe.gelb,'LineStyle','none');  
%    set(xx,'facealpha',.2);
%    xx=fill(X.coh4,Y.coh4,farbe.lila,'LineStyle','none');  
%    set(xx,'facealpha',.3);
%    xx=fill(X.coh5,Y.coh5,farbe.orange,'LineStyle','none');  
%    set(xx,'facealpha',.2);
%    xx=fill(X.coh6,Y.coh6,farbe.hblau,'LineStyle','none');  
%    set(xx,'facealpha',.2);
   set(gca, 'FontSize', 30)
   grid on; zoom on;
   ylim([0 1]);
   xlim([1e-2 60]);
   
   
      fig = gcf;
   fig.PaperUnits = 'inches';
   fig.PaperPosition = [0 0 22 13];
   print(fig,'-dpng','-r0',['F:\ALICE_Zarnekow_2018-05-23\Mittelung\hamm2048-Koh-', Abschnitt(i).name, 'e_V' '.png']) 
 

%% Plots Feuchte gemittelt

%    q=figure('visible','on','Name',strcat('Spektren Mischungsverhältnis'));
%    set(q, 'Units', 'normalized', 'Position',[0.1, 0.1, 0.8, 0.8] );
%    AxesH = axes('Units', 'normalized', 'Position', [0.12,0.145, 0.85, 0.74]);
%    loglog( F6h , MVspektrum.ges1, 'color' ,farbe.dblau,'LineWidth',2);
%    hold on;
%    loglog( F1 , MVspektrum.ges2, 'color' ,farbe.rot,'LineWidth',2);
%    hold on;
%    loglog( freq , fl, 'k--');
%    hold on;
%    loglog( freq , fk, 'k--');
%    hold on;
%    title(strcat({'ALICE Steigflüge'}), 'FontSize', 40); 
% %    title(strcat({[' Leistungsdichtespektrum Mischungsverhältnis ' ];[ 'Mittelung über ' int2str(Anz.ges) ' Fenster ' Abschnitt(i).name ' - ALICE Zarnekow']}), 'FontSize', 40); 
%    h=legend( 'Humicap HMP110','Rapid P-14');
%    set(h,'location','northeast','FontSize',25)
%    xlabel('Frequenz [Hz]', 'FontSize', 35);
%    ylabel('Feuchte [g/kg]^2s', 'FontSize', 35);
% %    xx=fill(Xmv.T1,Ymv.T1,farbe.dblau,'LineStyle','none');  
% %    set(xx,'facealpha',.2);
% %    xx=fill(Xmv.T2,Ymv.T2,farbe.rot,'LineStyle','none');  
% %    set(xx,'facealpha',.2);
%    set(gca, 'FontSize', 30)
%    grid on; zoom on;
%    xlim([1e-3 60]);
%    ylim([1e-10 1e4]); 
%    
%    fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['F:\ALICE_Zarnekow_2018-05-23\Mittelung\', Abschnitt(i).name, 'e_MV_V' '.png'])

   
   
%    q=figure('visible','on','Name',strcat('Spektren Kohärenz'));
%    set(q, 'Units', 'normalized', 'Position',[0.1, 0.1, 0.8, 0.8] );
%    AxesH = axes('Units', 'normalized', 'Position', [0.1,0.125, 0.87, 0.74]);
%    semilogx(F1 , P1mvcohs,'color', farbe.dgruen,'LineWidth',2);hold on;
%    title(strcat( {[' Kohärenzspektrum Mischungsverhältnis - Auf P-14 '];[' Mittelung über ' int2str(Anz.ges) ' Fenster ' Abschnitt(i).name ' - ALICE Zarnekow']}), 'FontSize', 40); 
%    h=legend('Humicap HMP110');
%    set(h,'location','southwest','FontSize',25)
%    xlabel('Frequenz [Hz]', 'FontSize', 35); 
%    ylabel('Kohärenz', 'FontSize', 35); 
%    xx=fill(Xmv.coh1,Ymv.coh1,farbe.dgruen,'LineStyle','none');  
%    set(xx,'facealpha',.2);
%    grid on; zoom on;
%    set(gca, 'FontSize', 30)
%    ylim([0 1]);
%    xlim([1e-2 60]);
%    
%    
%       fig = gcf;
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 22 13];
%    print(fig,'-dpng','-r0',['G:\ALICE_Zarnekow_2018-05-23\Mittelung\hamm2048-Koh-', Abschnitt(i).name, 'e_MV' '.png']) 
 
end



