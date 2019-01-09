% Winddiskretisierung

v_W_min = 0;
v_W_Delta = 5;
v_W_max = 30;

% Massendiskretisierung

m_min = 1;
m_Delta = 1.5;
m_max = 10;

% Geschwindigkeitsdiskretisierung

V_Kg_min = 1;
V_Kg_Delta = 1;
V_Kg_max = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisierungen

% Batterie

U_Bat_nom = N_Bat_cell * U_Bat_cell;        % nominale Batteriespannung
U_Bat_min = N_Bat_cell * U_Bat_cell_min;    % minimale Batteriespannung
% C_Bat = E_Dichte * m_Bat / U_Bat_nom;       % Kapazitaet der Batterie in As
C_Bat = N_Bat_cell_p*C_Bat_cell*3600;
Delta_C_Bat = 0;                            % Initialisierung Batteriekapazität, die nach jedem delta_h gebraucht wird



% Propeller

%Entnahme des Durchmessers und des Pitches aus dem Propellernahmen
D = str2double(prop_name(1:strfind(prop_name,'x')-1));      % Durchmesser extrahieren
P_75 = prop_name(strfind(prop_name,'x')+1:end);             % Pitch extrahieren
while isnan(str2double(P_75)) == 1
        P_75(end) = [];
end
P_75 = str2double(P_75);                    % Pitch festlegen
R = D * 0.0254 / 2;                         % Propellerradius in Meter
F = pi * R^2;                               % Fläche eines Propellers in Quadratmeter
Theta_75 = atan( 4*P_75 / (3*pi * D) );     % geometrischer Anstellwinkel des Propellers bei 75% des Radius
[RPM_map, V_map, T_map, P_map, TAU_map] = Propeller_map(DATA_APC,prop_name);    % Aufbau des Kennfeldes





% Umgebungsparameter

T_11 = T_0 - 0.0065 *(11000-H_0);                   % T in 11000m Höhe                     
rho_11 = rho_0 * (1 - 0.0065*(11000/T_0))^4.256;    % Dichte in 11000m Höhe
p_11 = p_0 * (1 - 0.0065*(11000/T_0))^5.256;        % Druck in 11000m Höhe



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for v_W_variabel = v_W_min:v_W_Delta:v_W_max;

	v_W = v_W_variabel;	
	
	lenghtm = floor(abs(m_min-m_max) / m_Delta + 1);
	m_flugsystem = zeros(lengthm,1);
	h1 = zeros(lengthm,1);	


	p = 1;

	for m_variabel = m_min:m_Delta:m_max

		m_flugsystem(p) = m_variabel;
		

		


		lengthh = floor(abs(H_0 - H_max) / Delta_H + 1);
		V_Kg_opt = zeros(lengthh,1);
		C_Rest_V = zeros(lengthh,1);
		U_mot = zeros(lengthh,1);
		I_mot = zeros(lengthh,1);
		Omega = zeros(lengthh,1);
		C_Rate = zeros(lengthh,1);
		I_Bat = zeros(lengthh,1);
		PWM = zeros(lengthh,1);
		H = zeros(lengthh,1);			
		Thrust = zeros(lengthh,1);	
		tau = zeros(lengthh,1);
		alpha = zeros(lengthh,1);	
		Theta = zeros(lengthh,1);
		rho = zeros(lengthh,1);	
		P_Bat = zeros(lengthh,1);	
		M_tip = zeros(lengthh,1);
		eta_prop = zeros(lengthh,1);
		eta_ges = zeros(lengthh,1);
		
		i = 1;

		for h = H_0:Delta_H:H_max
		

			% ursprgl. Programm: Berechnung wichtiger Größen, Berechnung Schub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Umgebungsparameter als Funktion der Höhe
    			% Berechnung der Flughöhe für iterative Schritte und der mittleren
    			% Dichte zwischen den Diskretisierungspunkten
    
    			H_unten = h;
    			H_oben = h + Delta_H;
    			H_mitte = (H_oben + H_unten)/2;
    
    			% Berechnung der Dichte an den Intevallgrenzen nach Normatmosphärenbedingungen
    
    			if H_oben <= 11000
        
        			rho_unten = rho_0 *(1-0.0065*(H_unten/T_0))^4.256;
        			rho_oben = rho_0 *(1-0.0065*(H_oben/T_0))^4.256;
        
        			p_unten = p_0 *(1-0.0065*(H_unten/T_0))^5.256;
        			p_oben = p_0 *(1-0.0065*(H_oben/T_0))^5.256;
        
    			elseif H_unten <= 11000 && H_oben >11000
        
        			rho_unten = rho_0 *(1 - 0.0065*(H_unten/T_0))^4.256;
        			rho_oben = rho_11 * exp(-g/(287*T_11)*(H_oben-11000));
        
        			p_unten = p_0 *(1-0.0065*(H_unten/T_0))^5.256;
        			p_oben = p_11 * exp(-g/(287.1*T_11)*(H_oben-11000));
        
    			else
        
        			rho_unten = rho_11 * exp(-g/(287*T_11)*(H_unten-11000));
        			rho_oben = rho_11 * exp(-g/(287*T_11)*(H_oben-11000));
        
        			p_unten = p_11 * exp(-g/(287.1*T_11)*(H_unten-11000));
        			p_oben = p_11 * exp(-g/(287.1*T_11)*(H_oben-11000));
        
    			end
    
    			if H_mitte <= 11000                                             % mittlere Temperatur im Intervall
        			T = T_0 - 0.0065 * H_mitte;        
    			else        
        			T = T_11;
    			end
    
    			rho(i) = rho_unten + (rho_oben - rho_unten)/2;                  % Berechnung der mittleren Dichte im Intervall
    			p = (p_oben + p_unten)/2;                                       % mittlerer Druck im Intervall
    			a = sqrt(kappa*p/rho(i));                                       % Schallgeschwindigkeit in m/s
    

    			% Dichte an der oberen (_2) und unteren (_1) Intervallgrenze
    			if i == 1
        			rho_1 = rho_0;
        			rho_2 = rho(i);
    			else
        			rho_1 = rho(i-1);
        			rho_2 = rho(i);
    			end
    
    
    			T_map = T_map * rho(i)/rho_1;                                   % Anpassung des Schubkennfeldes an die sich ändernde Dichte
    			P_map = P_map * rho(i)/rho_1;                                   % Anpassung des Leistungskennfeldes an die sich ändernde Dichte
    			TAU_map = TAU_map * rho(i)/rho_1;                               % Anpassung des Drehmomentkennfeldes an die sich ändernde Dichte
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
			lengthvkg = floor(abs(V_Kg_min - V_Kg_max) / V_Kg_Delta + 1);
			V_Kg = zeros(lengthvkg,1);
			
			% Hier Vektor der Variablen zum Auswählen der opt. Steiggeschw. speichert

			z = 1;

			for V_Kg_variabel = V_Kg_min:V_Kg_Delta:V_Kg_max
			
				V_Kg(z) = V_Kg_variabel;

				t_Flug = Delta_H/V_Kg(z);
				
				% Berechne hier alle weiteren Größen des Flugsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Abfrage_Flugsystem == 1                                              % handelt es sich um Multicopter oder Flächenflugzeug
    
    % MULTICOPTER
    m = m_copter + m_Bat + m_Mot * n_Prop + m_nutz;                     % Gesamtmasse des Quadrocopters
    
    % Aerodynamik berechnen
    
    [Thrust(i),Theta(i),V_A,alpha(i)] = MulticopterAerodynamik(u_Wg,V_Kg,gamma,c_W_copter_seitlich,c_W_copter_oben,c_A_copter_max,rho(i),A_copter,m,g);
    
else
    
    % FLÄCHENFLUGZEUG
    m = m_flugzeug + m_Bat + m_Mot * n_Prop + m_nutz;                   % Gesamtmasse des Flächenflugzeugs
    
    % Aerodynamik
    
    [Thrust(i),V_A] = FlaechenflugzeugAerodynamik(m,g,epsilon,V_Kg);
    
    
end
    
    
    				Thrust(i) = Thrust(i) / n_Prop;                                         % Schub auf n Propeller verteilen
    
    
    
    				if Thrust(i) > max(max(T_map))                                          % wenn Schub zu gross (Ergebnis verwerfen)
        				Omega(i) = NaN;
        				I_mot(i) = NaN;
        				C_Rate(i) = NaN;
        				C_Rest_V(i) = NaN;
    				else
        
        
       	 				% Drehzahl und Drehmoment bestimmen
        
        				[Omega(i),tau(i)] = Propeller(V_A, alpha(i), Thrust(i), RPM_map, V_map, T_map, TAU_map);
        
        
        				% Wie groß ist die Blattspitzengeschwindigkeit?
        
        				M_tip(i) = (Omega(i) * R)/a;                                        % Blattspitzengeschwindigkeit in Ma
        
        
        				% Motorzustand berechnen
        
        				[U_mot(i),I_mot(i)] = Motor(tau(i),K_V,I_0,R_i,Omega(i));
        
        
        				% Zustand der Motorregler berechnen
        
        				[PWM(i),eta_PWM] = ESC(U_mot(i),U_Bat_nom);
        
        
        				% Batteriezustand berechnen
        
        				[I_Bat(i),C_Rate(i),Delta_C_Bat,C_Rest_V(i)] = Batterie(PWM(i),eta_PWM,I_mot(i),n_Prop,C_Bat,P_Bat_Peukert,Delta_C_Bat,t_Flug);
        
        
        				% Gesamtwirkungsgrad       
        
        				% Berechnung der induzierten Geschwindigkeiten nach van der Wall
        				% (Grundlagen der Hubschrauber-Aerodynamik) (2015) S. 153
        
        				vi0 = sqrt(m*g / ( 2*rho(i)*F*n_Prop ) );                          % induzierte Geschwindigkeit im Schwebeflug v_i0 
        				v = vi0;
        				mu_z = -V_A*sin(alpha(i));                                          % Geschwindigkeit durch die Rotorebene
        				mu = V_A*cos(alpha(i));                                             % Geschwindigkeit entlang Rotorebene
        				krit = 1;
        				while krit > 0.0005
            					f = v - mu_z - vi0^2 / sqrt(mu^2 + v^2);
            					fs = 1 + v * vi0^2 / (mu^2 + v^2)^(3/2);
            					v_i_neu = v - f/fs;
            					krit = abs(v_i_neu - v) / v_i_neu;
            					v = v_i_neu;
        				end
        				vi_vi0 = (v - mu_z) / vi0;
        				vi = vi0 * vi_vi0;                                                  % induzierte Geschwindigkeit im stationaeren Steigflug
      
        				% Figure of Merit des Rotors, Bezug auf van der Wall (Grundlagen der Hubschrauber-Aerodynamik) (2015) (S.122)
        				eta_prop(i) = (Thrust(i) * (V_A + vi))/(tau(i) .* Omega(i));  
        
        				eta_ges(i) = (n_Prop * Thrust(i) * (mu_z + vi))/(I_Bat(i) * U_Bat_nom);         % Leistung, die in Schub umgesetzt wird im Verhältnis zur aufgebrachten Leistung
        
        
    				end
    
    
    
    
    				% Werden Grenzen ueberschritten?
    
   				% Flugbereichsgrenzen für das Flächenflugzeug innerhalb der
   				% Flugenvellope
    
    				if Abfrage_Flugsystem == 0
        
        				% aerodynamische Grenze
       					% herausnahmen, Annahme Flug mit V* und nicht V_min, i.e. epsilon_opt
        				n_z = cos(atan(epsilon));
        				V_min = sqrt(2*n_z*m*g/(c_A_plane_max*S*rho(i)));
        
        				% Leistungsgrenze
        
        				% c_A = C_A0 + alpha * C_Aalpha;
        				% c_W = c_W0 + k * C_A^2
        				W = V_A^2 * rho(i)/2 * S * c_W_plane;
        				T_erf = W;
        
        				% Temperaturgrenze
        
        				T_zul = 273.15 + t_zul;
        				% T_max = 2 * T_zul / ((kappa - 1) * (V_A/a)^2 +2);
        				V_max_T = sqrt((T_zul-T)*2/(T*(kappa-1))) * a;
        
        				% Begrenzung durch Festigkeit
        
        				q_max = rho(i)/2 * V_A^2;
        				V_max_q = sqrt(q_zul * 2 /rho(i));
        
        				if  V_A > V_max_q || V_A >= V_max_T || V_min > V_A || T_erf > max(max(T_map)) % || T_max > T_zul
            					C_Rest_V(i) = NaN;
            					Omega(i) = NaN;
            					U_mot(i) = NaN;
            					I_mot(i) = NaN;
            					I_Bat(i) = NaN;
            					PWM(i) = NaN;
            					eta_ges(i) = NaN;
        				end
    				end
    
    
    				% Wenn Grenzen ueberschritten werden, Resultate entfernen
    	
    				if C_Rest_V(i) < 0.0 || U_mot(i) > U_Bat_nom || U_mot(i) <= 0 || C_Rate(i) > C_Rate_max || I_mot(i) > I_max || alpha(i) > alpha_stall || M_tip(i) >= 1
 	       				C_Rest_V(i) = NaN;
        				Omega(i) = NaN;
        				U_mot(i) = NaN;
        				I_mot(i) = NaN;
        				I_Bat(i) = NaN;
        				PWM(i) = NaN;
        				eta_ges(i) = NaN;
    				end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				z = z + 1;
				H(i) = H_oben;			% Speichern der Höhe im Vektor
    				i = i+1;


			end
		
			% Hier Kriterium für optimale Steiggeschwindigkeit 

			ind_opt = find();			% Index der opt. Steiggeschw. finden

			V_Kg_opt(h) = V_Kg(ind_opt);		% Bestimmung opt. Stgeschw. und speichern in Vektor für jeden Höhenabschnitt

		end
		
		

		% Plot Steiggeschwindigkeit über der Höhe (Kurvenschar)
		figure(figure_V_Kg)
		subplot(4,2,v_W_variabel/v_W_Delta)
		h1(p) = plot(H,V_Kg_opt,'LineWidth',2);
		legendInfo{p} = ['m_{Bat} = ' num2str(m_flugsystem) 'kg'];		
		grid on 
		hold on 

		% Plot Restladung über der Höhe (Kurvenschar)
		figure(figure_C_Rest_V)
		subplot(4,2,v_W_variabel/v_W_Delta)
		plot(H,C_Rest_V,'LineWidth',2)
		grid on 
		hold on

		% Plot Drehzahl über der Höhe (Kurvenschar)
		figure(figure_omega)
		subplot(4,2,v_W_variabel/v_W_Delta)
		plot(H,omega,'LineWidth',2)
		grid on 
		hold on

		% Plot PWM über der Höhe (Kurvenschar)
		figure(figure_PWM)
		subplot(4,2,v_W_variabel/v_W_Delta)
		plot(H,PWM,'LineWidth',2)
		grid on 
		hold on

		p = p + 1;

	end
	
	% Beschriftung der Diagramme
	
	% Achenbeschriftung Steiggeschwindigkeit
	figure(figure_V_Kg)
	title(['v_W = ' num2str(v_W_variabel) 'm/s']) % , V_{Wind} = ' num2str(m]);
	legend([h1], legendInfo)
	xlabel('Höhe [m]')
	ylabel('Steiggeschwindigkeit [m/s]')

	% Achenbeschriftung Restladung
	figure(figure_C_Rest_V)
	xlabel('Höhe [m]')
	ylabel('Restladung [%]')

	% Achenbeschriftung RPM
	figure(figure_omega)
	xlabel('Höhe [m]')
	ylabel('Drehzahl [U/min]')

	% Achenbeschriftung PWM
	figure(figure_PWM)
	xlabel('Höhe [m]')
	ylabel('Pulsweitenmodulation [%]')


end
