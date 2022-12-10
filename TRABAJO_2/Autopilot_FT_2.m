function [FT_AP_CL_2, FT_AP_OL_2] = Autopilot_FT_2(FT_AP_CL_1,G_gyro,K_P_phi_p)
    % Defino el sensor para phi
    G_sphi = G_gyro;
    
    %Asigno a nuevas variables las del inner loop
    G_beta_deltaA   = FT_AP_CL_1.beta ;
    G_phi_deltaA    = FT_AP_CL_1.phi ;
    G_r_deltaA      = FT_AP_CL_1.r ;
    G_p_deltaA      = FT_AP_CL_1.p ;

    %Defino la open loop a mano
    FT_AP_OL_2 = G_sphi*G_phi_deltaA*K_P_phi_p;

    %Defino las close loop a mano del autopiloto
    DEN = 1 + FT_AP_OL_2;  %El denominador, que es común en todas
    FT_AP_CL_2.phi = minreal(G_phi_deltaA*K_P_phi_p/DEN,0.001);
    FT_AP_CL_2.beta = minreal((G_beta_deltaA*K_P_phi_p + G_beta_deltaA*...
        G_phi_deltaA*K_P_phi_p^2*(G_sphi - 1))/DEN,0.001);
    FT_AP_CL_2.r = minreal((G_r_deltaA*K_P_phi_p + G_r_deltaA*...
        G_phi_deltaA*K_P_phi_p^2*(G_sphi - 1))/DEN,0.001);
    FT_AP_CL_2.p = minreal((G_p_deltaA*K_P_phi_p + G_p_deltaA*...
        G_phi_deltaA*K_P_phi_p^2*(G_sphi - 1))/DEN,0.001);
    
end