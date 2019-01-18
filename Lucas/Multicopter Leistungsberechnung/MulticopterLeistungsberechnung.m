% Initialisierung: Parameterberechnung 

U_Bat_nom = N_Bat_cell * U_Bat_cell;        % nominale Batteriespannung
U_Bat_min = N_Bat_cell * U_Bat_cell_min;    % minimale Batteriespannung
R = D * 0.0254 / 2;                         % Propellerradius in Meter
F = pi * R^2;                               % Fläche eines Propellers in Quadratmeter
Theta_75 = atan( 4*P_75 / (3*pi * D) );     % geometrischer Anstellwinkel des Propellers bei 75% des Radius

l = 0;

for V_variabel2 = V_Kg_2_min:V_Kg_2_Delta:V_Kg_2_max;             % Sinkgeschwindigkeit (positiv) in m/s

    % Initialisierungen
    lengthj = floor(abs(m_Bat_max - m_Bat_min) / m_Bat_Delta + 1);
    C_Rest_V_max = zeros(lengthj,1);
    V_1_opt = zeros(lengthj,1);
    h1 = zeros(lengthj,1);
    Theta_maxi = zeros(lengthj,1);

    j = 1;                                              % Initialisierung des Zaehlers j

    for m_Bat = m_Bat_max:-m_Bat_Delta:m_Bat_min         % Fuer alle Batteriemassen


        % Masse und induzierte Geschwindigkeit im Schwebeflug berechnen

        m = m_copter + m_Bat;                           % Gesamtmasse des Quadrocopters
        C_Bat = E_Dichte * m_Bat / U_Bat_nom;           % Kapazitaet der Batterie in As
        w_i0_f = sqrt(m*g / ( 2*rho*F*n_Prop ) );          % induzierte Geschwindigkeit im Schwebeflug v_i0


        % Initialisierungen
        lengthi = floor(abs(V_Kg_1_min - V_Kg_1_max) / V_Kg_1_Delta + 1);
        C_Rest_V = zeros(lengthi,1);
        U_mot_max = zeros(lengthi,1);
        U_mot_min = zeros(lengthi,1); 
        I_mot_max = zeros(lengthi,1); 
        C_Rate_m = zeros(lengthi,1);
        omega_max = zeros(lengthi,1); 
        Theta_max = zeros(lengthi,1); 
        alpha_max = zeros(lengthi,1);

        i = 1;

        for V_variabel1 = V_Kg_1_min:V_Kg_1_Delta:V_Kg_1_max          % Fuer alle Steiggeschwindigkeiten

            % Initialisierung
            Delta_C_Bat = 0;
            U_mot_0 = zeros(3,1);
            U_mot = zeros(3,1);
            I_mot_0 = zeros(3,1);
            I_mot = zeros(3,1);
            C_Rate = zeros(3,1);
            Thrust = zeros(3,1);
            Theta = zeros(3,1);
            omega_0 = zeros(3,1);
            omega = zeros(3,1);
            Throttle = zeros(3,1);
            alpha_75 = zeros(3,1);


            for k = 1:length(Abfrage_hovern)                                % Fuer alle Flugphasen

                % Erkennen der aktuellen Flugphase aus definierter Mission
                if Abfrage_hovern(k) == 0
                    if Abfrage_V_variabel1(k) == 0 && Abfrage_V_variabel2(k) == 0
                        V_Kg = Bahngeschwindigkeit(k);
                    elseif Abfrage_V_variabel1(k) ~= 0
                        V_Kg = V_variabel1;
                    else
                        V_Kg = V_variabel2;
                    end
                    t_Flug = Strecke(k) / V_Kg; 
                else
                    V_Kg = 0;
                    t_Flug = t_hover(k);
                end


                % Aerodynamik berechnen und
                % Schub berechnen 

                [Thrust(k),Theta(k),V_A,alpha] = Aerodynamik(u_Wg(k),V_Kg,gamma(k),c_W_copter_seitlich,c_W_copter_oben,c_A_copter_max,rho,A_copter,m,g);
                Thrust(k) = Thrust(k) / n_Prop;


                if Thrust(k) > max(Thrust_static)               % Schub zu gross (Ergebnis verwerfen)
                    omega_0(k) = NaN;
                    I_mot_0(k) = NaN;
                    C_Rate(k) = NaN;
                    C_Rest_V(i) = NaN;
                else


                % Motorzustand aus Kennfeld im Schwebeflug interpolieren

                    omega_0(k) = interp1(Thrust_static,omega_static,Thrust(k),'pchip');
                    I_bat_0 = n_Prop * interp1(Thrust_static,Amps_static,Thrust(k),'pchip');                   
                    Throttle(k) = interp1(Thrust_static,Throttle_static,Thrust(k),'pchip');
                    U_mot_0(k) = 0.01 * Throttle(k) * U_Bat_nom;
                    
                    PWM_0 = U_mot_0(k) / U_Bat_nom;
                    if PWM_0 > 0 && PWM_0 < 0.5 
                        eta_PWM = 0.7 * PWM_0 + 0.50;
                    elseif PWM_0 >= 0.5 && PWM_0 <= 1
                        eta_PWM = 0.2 * PWM_0 + 0.75;
                    else
                        eta_PWM = 1;
                    end
                    I_mot_0(k) = I_bat_0 / PWM_0 * eta_PWM / n_Prop;


                % Induzierte Geschwindigkeit berechnen,
                % Drehzahl neu bestimmen und
                % Motorzustand durch Steiggeschwindigkeit neu bestimmen

                [U_mot(k),I_mot(k),omega(k),alpha_75] = Antrieb(Thrust(k),omega_0(k),w_i0_f,V_A,alpha,R,Theta_75,rho,c_d0,a,I_mot_0(k),I_0,U_mot_0(k),K_V,R_i);
                %U_mot(k) = U_mot_0(k);
                %I_mot(k) = I_mot_0(k);
                

                % Zustand der Motorregler berechnen

                [PWM,eta_PWM] = ESC(U_mot(k),U_Bat_nom);


                % Batteriezustand berechnen

                [I_Bat,C_Rate(k),Delta_C_Bat,C_Rest_V(i)] = Batterie(PWM,eta_PWM,I_mot(k),n_Prop,C_Bat,P_Bat,Delta_C_Bat,t_Flug); 

                end

            end


            % Werden Grenzen ueberschritten?

            % maximal auftretende Groessen ermitteln
            U_mot_max(i) = max(U_mot); U_mot_min(i) = min(U_mot); I_mot_max(i) = max(I_mot); C_Rate_m(i) = max(C_Rate);
            omega_max(i) = max(omega); Theta_max(i) = max(Theta); alpha_max(i) = max(alpha_75);

            % Wenn Grenzen ueberschritten werden, Resultate entfernen
            if C_Rest_V(i) < 0.0 || U_mot_max(i) > U_Bat_nom || U_mot_min(i) <= 0 || C_Rate_m(i) > C_Rate_max || I_mot_max(i) > I_max || alpha_max(i) > alpha_stall
                C_Rest_V(i) = NaN;
                omega_max(i) = NaN;
                U_mot_max(i) = NaN;
                I_mot_max(i) = NaN;
            end

            i = i + 1;                                          % inkrementieren (Geschwindigkeit)

        end

        %% Plot

        V_1 = V_Kg_1_min:V_Kg_1_Delta:V_Kg_1_max;                     % Vektor x-Achse
        C_Rest_V_max(j) = max(C_Rest_V);                        % Maximum der Restladung ermitteln...
        if ~isnan(C_Rest_V_max(j))                              % ... wenn Fall existiert
            ind_opt = find(C_Rest_V == max(C_Rest_V));
            V_1_opt(j) = V_1(ind_opt);
        else
            V_1_opt(j) = NaN;                               % ... wenn Fall nicht existiert, entfernen
        end

        % Plot Restladung ueber Geschwindigkeit (Kurvenschar)
        figure(figure_C_Rest_V)
        subplot(3,1,V_variabel2/2)
        h1(j) = plot(V_1,C_Rest_V*100,'LineWidth',2);
        legendInfo{j} = ['m_{Bat} = ' num2str(m_Bat) 'kg'];
        grid on
        hold on

        % Plot max. RPM ueber Geschwindigkeit (Kurvenschar)
        figure(figure_omega)
        subplot(3,1,V_variabel2/2)
        plot(V_1,omega_max/(2*pi)*60,'LineWidth',2)
        grid on
        hold on  

        % Plot max. I_mot ueber Geschwindigkeit (Kurvenschar)
        figure(figure_I_mot)
        subplot(3,1,V_variabel2/2)
        plot(V_1,I_mot_max,'LineWidth',2)
        grid on
        hold on

        % Plot max. U_mot ueber Geschwindigkeit (Kurvenschar)
        figure(figure_U_mot)
        subplot(3,1,V_variabel2/2)
        plot(V_1,U_mot_max,'LineWidth',2)
        grid on
        hold on


        Theta_maxi(j) = max(Theta_max);
        
        display(['(',num2str(round((j + l*lengthj)/(lengthj * 3) * 100)),' %)'])

        j = j + 1;                                              % inkrementieren (Batteriemasse)
        
        

    end


    % Plot Maximum der Restladung
    figure(figure_C_Rest_V)
    subplot(3,1,V_variabel2/2)
    h2=plot(V_1_opt,C_Rest_V_max*100,'k-o','LineWidth',2);
    % maximales Theta anzeigen
    Theta_diag = max(Theta_maxi) * 180 / pi;
    text(16,25,['\Theta`_{min} = ' num2str(round(Theta_diag)) '°']);

    % Achsenbeschriftung Restladung
    figure(figure_C_Rest_V)
    title(['V_{Kg,2} = ' num2str(V_variabel2) 'm/s , V_{Wind,1} = ' num2str(u_Wg(1)) 'm/s'])
    legend([h1],legendInfo)
    xlim([0 V_Kg_1_max+5])
    ylim([20 80])
    xlabel('V_{Kg,1} [m/s]')
    ylabel('Restladung [%]')

    % Achsenbeschriftung RPM
    figure(figure_omega)
    xlabel('V_{Kg,1} [m/s]')
    ylabel('RPM')

    % Achsenbeschriftung RPM
    figure(figure_I_mot)
    xlabel('V_{Kg,1} [m/s]')
    ylabel('I_{mot} [A]')

    % Achsenbeschriftung RPM
    figure(figure_U_mot)
    xlabel('V_{Kg,1} [m/s]')
    ylabel('U_{mot} [V]')
    
    l = l + 1;

end


  %% Datei abspeichern
ImageSizeX = 14;
ImageSizeY = 24;
figure(figure_C_Rest_V)
set(gcf,'PaperUnits','centimeters', 'PaperPosition', [0 0 ImageSizeX ImageSizeY]); 
set(gcf,'Units','centimeters', 'PaperSize', [ImageSizeX ImageSizeY]); 
saveas(gcf,Dateiname, 'pdf');   