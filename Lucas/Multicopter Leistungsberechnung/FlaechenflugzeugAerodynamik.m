function [Thrust,Theta_1,V_A,alpha] = FlaechenflugzeugAerodynamik(u_Wg,v_Wg,w_Wg,V_Kg,gamma,C_A0,C_Aalpha,epsilon,m,g)

if u_Wg == 0
    u_Wg = 0.01;
end 

if v_Wg == 0
    v_Wg = 0.1;
end

if w_Wg == 0
    w_Wg = 0.1;
end

u_Kg = sin(gamma) * V_Kg;
w_Kg = - cos(gamma) * V_Kg;

V_A = sqrt((u_Kg+u_Wg)^2 + (v_Wg)^2 + (w_Kg+w_Wg)^2);
gamma_a = atan(v_Wg / (u_Kg+u_Wg));
chi_a = atan((w_Kg+w_Wg) / sqrt((u_Kg+u_Wg)^2+ v_Wg^2));
Theta_1 = 0;
Delta_Theta = 2;
i=0;
    while Delta_Theta > 0.001*pi/180;
        alpha = Theta_1 - gamma_a;
        c_A = C_A0 + C_Aalpha * alpha;
        c_W = epsilon * c_A;
        W = c_W * rho/2 * V_A^2 * S;
        A = c_A * rho/2 * V_A^2 * S;
        X_g = - W * cos(gamma_a) * cos(chi_a) + A * sin(gamma_a);
        Y_g = W * sin(chi_a);
        Z_g = - W * sin(gamma_a) * cos(chi_a) - A * cos(gamma_a) + m*g;
        Theta_2 = -atan(-X_g / Z_g);
        Delta_Theta = abs(Theta_2 - Theta_1);
        Theta_1 = Theta_2;  
        i = i + 1;
        if i > 100
            break
        end
    end
    Thrust = sqrt(X_g^2 + Y_g^2 + Z_g^2);
    
end


