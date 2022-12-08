%% Función que calcula las FT del sistema realimentado
function [FT_CL, FT_OL] = Aumented_FT(F_beta,F_r,G_act,G_gyro,G_vane,G_wash,K_DL,p,FT_lat)
    % Calculamos todas las FT a mano en lazo cerrado
        % Funciones de transferencia de los elementos
    Ga_deltaA = G_act; Ga_deltaR = G_act;       % FT actuadores
    Gs_r = G_gyro; Gs_beta = G_vane;            % FT sensores
    Gf_r = G_wash;                              % FT filtro wash-out
    
    % Cálculo de las ganancias de realimentación en función de F
    K_deltaRbeta = -(F_beta - 1)*p.Cn_beta/p.Cn_deltaR; 
    K_deltaRr = -(F_r - 1)*(p.Cn_r/p.Cn_deltaR)*(0.5*p.b/p.Us);

    % Reasignamos a nuevas variables las FT de la planta libre
    G_betaDeltaA = FT_lat.fact.deltaA_beta; 
    G_betaDeltaR = FT_lat.fact.deltaR_beta; 
    G_rDeltaA = FT_lat.fact.deltaA_r; 
    G_rDeltaR = FT_lat.fact.deltaR_r; 
    G_phiDeltaA = FT_lat.fact.deltaA_phi; 
    G_phiDeltaR = FT_lat.fact.deltaR_phi; 
    G_pDeltaA = FT_lat.fact.deltaA_p; 
    G_pDeltaR = FT_lat.fact.deltaR_p;
    
        % Construcción de las FT en lazo cerrado 
    num_rDeltaA = K_DL*Ga_deltaA*G_rDeltaR+...
        K_DL*Ga_deltaA*Ga_deltaR*Gs_beta*K_deltaRbeta*(...
        G_rDeltaA*G_betaDeltaR-G_rDeltaR*G_betaDeltaA);
    den_rDeltaA = 1+Ga_deltaR*(Gs_beta*K_deltaRbeta*G_betaDeltaR+...
        Gf_r*Gs_r*K_deltaRr*G_rDeltaR); 
    FT_CL.r_deltaS =  num_rDeltaA/den_rDeltaA;
    
    num_betaDeltaA = K_DL*Ga_deltaA*G_betaDeltaA+...
        K_DL*Ga_deltaA*Ga_deltaR*Gf_r*Gs_r*K_deltaRr*(...
        G_betaDeltaA*G_rDeltaR-G_betaDeltaR*G_rDeltaA);
    den_betaDeltaA = 1+Ga_deltaR*(Gs_beta*K_deltaRbeta*G_betaDeltaR+...
        Gf_r*Gs_r*K_deltaRr*G_rDeltaR); 
    FT_CL.beta_deltaS =  num_betaDeltaA/den_betaDeltaA;
    
    FT_CL.p_deltaS = K_DL*Ga_deltaA*G_pDeltaA-...
        Ga_deltaR*Gf_r*Gs_r*K_deltaRr*G_pDeltaR*FT_CL.r_deltaS-...
        Ga_deltaR*Gs_beta*K_deltaRbeta*G_pDeltaR*FT_CL.beta_deltaS; 
    
    FT_CL.phi_deltaS = K_DL*Ga_deltaA*G_phiDeltaA-...
        Ga_deltaR*Gf_r*Gs_r*K_deltaRr*G_phiDeltaR*FT_CL.r_deltaS-...
        Ga_deltaR*Gs_beta*K_deltaRbeta*G_phiDeltaR*FT_CL.beta_deltaS;
    
    
    FT_CL.beta_deltaS =  minreal(FT_CL.beta_deltaS,0.001);
    FT_CL.r_deltaS =  minreal(FT_CL.r_deltaS,0.001);
    FT_CL.phi_deltaS =  minreal(FT_CL.phi_deltaS,0.001);
    FT_CL.p_deltaS =  minreal(FT_CL.p_deltaS,0.001);
    
        % FT en lazo abierto
    FT_OL = minreal(Ga_deltaR*(Gs_beta*K_deltaRbeta*G_betaDeltaR+...
        Gf_r*Gs_r*K_deltaRr*G_rDeltaR));
end