%% Untersuchung der Batterien %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear 
% close all 
% clc

%Matrix laden
load('Elektromodellflug.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NORMZELLE ERZEUGEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DATA = Elektromodellflug;

% Löschen der Ausreißer
DATA(63,:) = [];       % id_bat = 63
DATA(40,:) = [];       % id_bat = 40
DATA(30,:) = [];       % id_bat = 30
DATA(14,:) = [];       % id_bat = 14
DATA(38,:) = [];       % id_bat = 38





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FEHLERUNTERSUCHUNG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PWM = 0.80;
eta_PWM = 0.7;
I_mot = 8;
n_Prop = 4;
% C_Rate = 15;


%% Initialisierung für Fehlertoleranz
control = 0;
tolerance = zeros(length(DATA),1);
C_Rate_max = 50;                                        % Maximal zu untersuchende C_Rate festlegen
tolerance_crate = [1:length(DATA)]';       % Sammeln aller Abweichungen aller Batterien für jede C_Rate in 1er Schritten

for k = 1:1:C_Rate_max
    
    C_Rate = k;
    for l = 1:length(DATA)                 % für jede Batterie
        
        id_bat = l;
        
        
        
        %% set time step delta t in seconds
        step_dt = 1;
        

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %%  ORIGINALZELLE
        
        % DATA = evalin('base','Elektromodellflug');
        [Eo, A, K] = Batterie_parameter(cell2mat(DATA(id_bat,3)));
        
        % battery parameters
        Q = DATA{id_bat,3}(1);
        B = 3/DATA{id_bat,3}(3);
        R = DATA{id_bat,3}(8);
        N_el = DATA{id_bat,4};
        C = DATA{id_bat,5} / 1000;
        C_rate_max = DATA{id_bat,6};
        
        % initializing the variables
        i_int = 0;
        V_bat_1 = zeros(floor(3600*1.5/C_Rate),1);
        control = 0;
        bar = zeros(floor(3600*1.5/C_Rate),1);
        
        %     I_bat = PWM * I_mot / eta_PWM * n_Prop;
        I_bat = C_Rate*C;
        
        for n = 1:floor(3600*1.5/C_Rate)
            
            
            i_int = I_bat*step_dt + i_int;   % integral of the current
            
            
            % calculating the battery voltage (of one cell)
            V_bat_1(n) = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int + I_bat*0) + A * exp(-B*i_int);
            V_bat_1(n) = N_el * V_bat_1(n);   % the battery voltage of all cells
            
     
            
            if V_bat_1(n) < 3.1*N_el
                control = 1;
            end
            if control == 1
                V_bat_1(n) = NaN;
            end
            
            bar(n) = 3.1*N_el;
        end
             
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %% NORMZELLE
        
        [Q,Qnom,Qexp,Vfull,Vexp,Vnom,i,R] = Normcell(DATA);             % Normzelle generieren
        Cnom = DATA{id_bat,5}/1000;                                     % Kapazität
        B = 3/Qexp;
 
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% % % % %         R = 0.1824 /(Q*(0.06423*Q + 0.7362*C_Rate)^0.7725);             % MI MI MI hard coded
% % % %         R = 0.2853 / (Q*(0.1681*Q+ 1.426*C_Rate)^0.7549);
% % %       R = 0.2453 / (Q*(1.984* Q+ 2.253*C_Rate)^0.6245);                 % <-- bestes Ergebnis
% %        R = 0.09974 / (Q*(0.2548*Q+0.8755*C_Rate)^0.5424);                 % beide sehr gut
        R = 0.1077 / (Q*(0.1555*Q+0.9825*C_Rate)^0.5485);
%           R = 0.08974 / (Q*(0.1902*Q+0.724*C_Rate)^0.5419);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
        
        Batterie_data = [Q Qnom Qexp Vfull Vexp Vnom i R 0 0 0];        % Zwischenbelegung
        [Eo,A,K] = Batterie_parameter(Batterie_data);
        
        
        %% Vorgaben
        %     I_bat = PWM * I_mot / eta_PWM * n_Prop;                   % Batteriestrom
        I_bat = C*C_Rate;
        I_bat = I_bat/(Cnom);                                           % Normierung des Batteriestroms
        i_int = 0;
        V_bat_2 = zeros(floor(3600*1.5/C_Rate),1);
        control = 0;
        
        for n = 1:floor(3600*1.5/C_Rate)
            
            
            i_int = I_bat*step_dt/3600 + i_int;   % integral of the current
            
            
            % calculating the battery voltage (of one cell)
            V_bat_2(n) = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int + I_bat*0) + A * exp(-B*i_int);
            % the battery voltage of all cells
            V_bat_2(n) = N_el * V_bat_2(n)*1.0;  
            
            
            if V_bat_2(n) < 3.1*N_el
                control = 1;
            end
            if control == 1
                V_bat_2(n) = NaN;
            end
            
            
        end
        
        if k == 50 && l == 61
            aaa = 1;
        end
        % Berechnung der Toleranz nur bis zum ersten Mal, wenn V_bat unter
        % V_bat < 3.1 * N_el, ansonsten alle Ergebnisse entfernen
        flaeche1 = trapz(V_bat_1(~isnan(V_bat_1)));
        flaeche2 = trapz(V_bat_2(~isnan(V_bat_2)));
%         for n = floor(3600*1.5/C_Rate):-1:1
%             if isnan(V_bat_1(n)) == 1
%                 V_bat_1(n) = 1;
%             end
%             if isnan(V_bat_2(n)) == 1
%                 V_bat_2(n) = 1;
%             end
%         end

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        %% TOLERANZ 
        
        % Berechnung der Modellabweichung von der Originalzelle
        % Beachte dabei die max. C-Rate
        
        if DATA{id_bat,6} < C_Rate                          % Wenn C-Rate zu groß
            tolerance(l) = NaN;                               % Setze Abweichung auf 0
        else
            tolerance(l) = (flaeche1-flaeche2)/flaeche2*100;% ansonsten berechne Abweichung
        end
        
        
    end
    
    tolerance_crate = [tolerance_crate tolerance];
    % ^ Zeilen = Batterienummer(id_bat), Spalten = C_Rate in 1er Schritten,
    % Matrix = Abweichung Normzelle von Originalzelle
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    %% FEHLERTOLERANZ
    % Plotten der Abweichungen für jede Batterienummer
%     x = 1:length(DATA);
%     plot(x,tolerance)
%     grid on
%     xlabel('Batterienummer (Index in Matrix)');
%     ylabel('Durchschnittliche Abweichung in %');
%     average_tolerance = mean(tolerance);        % durchschnittliche Abweichung in % aller Batterien
end

