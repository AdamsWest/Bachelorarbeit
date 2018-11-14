function [Thrust,Theta_1,V_A,alpha] = MulticopterAerodynamik(u_Wg,V_Kg,gamma,c_W_copter_seitlich,c_W_copter_oben,c_A_max,rho,A_copter,m,g)

if u_Wg == 0
    u_Wg = 0.01;
end

w_Kg = V_Kg * sin(gamma);
u_Kg = V_Kg * cos(gamma);

V_A = sqrt((u_Kg + u_Wg)^2 + w_Kg^2);     % Absolute Flugwindgeschwindigkeit
gamma_a = atan(w_Kg / (u_Kg + u_Wg));
Theta_1 = 0;
Delta_Theta = 2;
i = 0;
    while Delta_Theta > 0.001*pi/180
        alpha = - gamma_a + Theta_1;
        c_W = - (c_W_copter_oben - c_W_copter_seitlich)/2 * cos(2*alpha) + (c_W_copter_oben + c_W_copter_seitlich)/2;
        c_A = c_A_max * sin(2*alpha);
        W = c_W * rho/2 * V_A^2 * A_copter;
        A = c_A * rho/2 * V_A^2 * A_copter;
        X_g = - W * cos(gamma_a) - A * sin(gamma_a);
        Z_g = W * sin(gamma_a) - A * cos(gamma_a) + m*g;
        Theta_2 = -atan(-X_g / Z_g);
        Delta_Theta = abs(Theta_2 - Theta_1);
        Theta_1 = Theta_2;  
        i = i + 1;
    end
    Thrust = sqrt(X_g^2 + Z_g^2);
    
end