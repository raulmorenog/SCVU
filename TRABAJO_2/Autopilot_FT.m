%% Función que calcula las FT del sistema realimentado
%La función saca las funciones de transferencia del autopiloto. Parte de
%las funciones del SAS, para eso llama la función 'Aumented_FT', suponiendo
%un K_DL = 1, porque no aparece en el autopiloto. 

function [FT_AP_CL_1, FT_AP_OL_1] = Autopilot_FT(F_beta,F_r,G_act,G_gyro,G_vane,G_wash,p,FT_lat,K_P_p_deltaA)

    % Llamamos a la función del SAS
    [FT_SAS_CL, FT_SAS_OL] = Aumented_FT(F_beta,F_r,G_act,G_gyro,G_vane,G_wash,1,p,FT_lat);
    
    % Defino el sensor para p
    G_sp = G_gyro;
    
    % Asigno a nuevas variables las del SAS
    G_beta_deltaA   = FT_SAS_CL.beta_deltaS;
    G_phi_deltaA    = FT_SAS_CL.phi_deltaS;
    G_r_deltaA      = FT_SAS_CL.r_deltaS;
    G_p_deltaA      = FT_SAS_CL.p_deltaS;

    % Defino la open loop a mano
    FT_AP_OL_1 = G_sp*G_p_deltaA*K_P_p_deltaA;

    % Defino las close loop a mano del autopiloto
    DEN = 1 + FT_AP_OL_1;  % El denominador, que es común en todas
    FT_AP_CL_1.p = minreal(G_p_deltaA*K_P_p_deltaA/DEN,0.001);
    FT_AP_CL_1.beta = minreal((G_beta_deltaA*K_P_p_deltaA + G_beta_deltaA*...
        G_p_deltaA*K_P_p_deltaA^2*(G_sp - 1))/DEN,0.001);
    FT_AP_CL_1.r = minreal((G_r_deltaA*K_P_p_deltaA + G_r_deltaA*...
        G_p_deltaA*K_P_p_deltaA^2*(G_sp - 1))/DEN,0.001);
    FT_AP_CL_1.phi = minreal((G_phi_deltaA*K_P_p_deltaA + G_phi_deltaA*...
        G_p_deltaA*K_P_p_deltaA^2*(G_sp - 1))/DEN,0.001);
    
end