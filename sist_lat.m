%% Funciones de Transferencia Lat-Dir
function lat = sist_lat(p)
        % Adimensionalización lateral-direccional y coeficientes
    mu_latdir = p.m/(p.rhos*p.Sw*p.b/2);
    I_xxnd = p.I_xx/(p.rhos*p.Sw*(p.b/2)^3);
    I_zznd = p.I_zz/(p.rhos*p.Sw*(p.b/2)^3);
    I_xznd = p.I_xx/(p.rhos*p.Sw*(p.b/2)^3);
    
    Cxs = p.Cxs;
    Czs = p.Czs;

    Cy_beta = p.Cy_beta;    
    Cy_p = p.Cy_p;    
    Cy_r = p.Cy_r;    
    Cy_deltaA = p.Cy_deltaA;    
    Cy_deltaR = p.Cy_deltaR;   
    
    Cl_beta = p.Cl_beta;    
    Cl_p = p.Cl_p;    
    Cl_r = p.Cl_r;    
    Cl_deltaA = p.Cl_deltaA;
    Cl_deltaR = p.Cl_deltaR;
    
    Cn_beta = p.Cn_beta;
    Cn_p = p.Cn_p;    
    Cn_r = p.Cn_r;    
    Cn_deltaA = p.Cn_deltaA;
    Cn_deltaR = p.Cn_deltaR;
    
        % Construcción del Sistema Lat-Dir: [C]*{Y}=[D]*{X} 
    syms s z
    sympref('FloatingPointOutput',true)
    
    C=[2*mu_latdir*s-Cy_beta,Czs-Cy_p*s,2*mu_latdir-Cy_r;...
        -Cl_beta,I_xxnd*s^2-Cl_p*s,-Cl_r-I_xznd*s;...
        -Cn_beta,-Cn_p*s-I_xznd*s^2,I_zznd*s-Cn_r];
    D=[0,Cy_deltaR;Cl_deltaA,Cl_deltaR;Cn_deltaA,Cn_deltaR];
    
    %Cramer y factorización con s ADIMENSIONAL (s gorro) --> FT con s
    %adimensional
    DenLat = sym2poly(det(C));    %Cuártica de estabilidad Lat-Dir
    
    Num_deltaA_beta = sym2poly(det([D(:,1),C(:,2),C(:,3)]));
    TFdeltaA_beta = tf([Num_deltaA_beta],[DenLat]); 

    Num_deltaR_beta = sym2poly(det([D(:,2),C(:,2),C(:,3)])); 
    TFdeltaR_beta = tf([Num_deltaR_beta],[DenLat]);                                               
    
    Num_deltaA_phi = sym2poly(det([C(:,1),D(:,1),C(:,3)]));
    TFdeltaA_phi = tf([Num_deltaA_phi],[DenLat]); 

    Num_deltaR_phi = sym2poly(det([C(:,1),D(:,2),C(:,3)]));  
    TFdeltaR_phi = tf([Num_deltaR_phi],[DenLat]);
    
    Num_deltaA_r = sym2poly(det([C(:,1),C(:,2),D(:,1)]));  
    TFdeltaA_r_nd = tf([Num_deltaA_r],[DenLat]);  

    Num_deltaR_r = sym2poly(det([C(:,1),C(:,2),D(:,2)]));      
    TFdeltaR_r_nd = tf([Num_deltaR_r],[DenLat]);

    TFdeltaA_r = TFdeltaA_r_nd*tf(2*p.Us/p.b);         
    TFdeltaR_r = TFdeltaR_r_nd*tf(2*p.Us/p.b);
    
    % Desadimensionalización de s (s = z*b/2Us) --> FT con s dimensional
    DenLatDim = sym2poly(subs(det(C),s,z*p.b/(2*p.Us)));
    
    Num_deltaA_betaDim = sym2poly(subs(det([D(:,1),C(:,2),C(:,3)]),s,z*p.b/(2*p.Us))); 
    TFdeltaA_betaDim = tf([Num_deltaA_betaDim],[DenLatDim]); 

    Num_deltaR_betaDim = sym2poly(subs(det([D(:,2),C(:,2),C(:,3)]),s,z*p.b/(2*p.Us)));
    TFdeltaR_betaDim = tf([Num_deltaR_betaDim],[DenLatDim]);
    
    Num_deltaA_phiDim = sym2poly(subs(det([C(:,1),D(:,1),C(:,3)]),s,z*p.b/(2*p.Us)));      
    TFdeltaA_phiDim = tf([Num_deltaA_phiDim],[DenLatDim]);     

    Num_deltaR_phiDim = sym2poly(subs(det([C(:,1),D(:,2),C(:,3)]),s,z*p.b/(2*p.Us)));
    TFdeltaR_phiDim = tf([Num_deltaR_phiDim],[DenLatDim]);
        
    Num_deltaA_rDim = sym2poly(subs(det([C(:,1),C(:,2),D(:,1)]),s,z*p.b/(2*p.Us)));
    TFdeltaA_rDim = tf([Num_deltaA_rDim],[DenLatDim])*tf(2*p.Us/p.b);

    Num_deltaR_rDim = sym2poly(subs(det([C(:,1),C(:,2),D(:,2)]),s,z*p.b/(2*p.Us)));                
    TFdeltaR_rDim = tf([Num_deltaR_rDim],[DenLatDim])*tf(2*p.Us/p.b);


    % Cálculo de las propiedades de las FT
        % Zeros, polos y ganancias
    [zeros_deltaA_beta,polos_deltaA_beta,ganancia_deltaA_beta] = tf2zp(Num_deltaA_betaDim ,DenLatDim);
    [zeros_deltaR_betaDim,polos_deltaR_betaDim,ganancia_deltaR_betaDim] = tf2zp(Num_deltaR_betaDim ,DenLatDim);  
    [zeros_deltaA_phiDim,polos_deltaA_phiDim,ganancia_deltaA_phiDim] = tf2zp(Num_deltaA_phiDim ,DenLatDim); 
    [zeros_deltaR_phiDim,polos_deltaR_phiDim,ganancia_deltaR_phiDim] = tf2zp(Num_deltaR_phiDim ,DenLatDim); 
    [zeros_deltaA_rDim,polos_deltaA_rDim,ganancia_deltaA_rDim] = tf2zp(Num_deltaA_rDim*2*p.Us/p.b ,DenLatDim); 
    [zeros_deltaR_rDim,polos_deltaR_rDim,ganancia_deltaR_rDim] = tf2zp(Num_deltaR_rDim*2*p.Us/p.b ,DenLatDim); 

        % Amortiguamientos y frecuencias de los modos
    [w_natural,amort,polos] = damp(TFdeltaA_betaDim);

    omega_blnc = w_natural(4);  
    chi_blnc = amort(4);
    t12_blnc = -0.693/real(polos(4));
    omega_blnch = w_natural(2); 
    chi_blnch = amort(2); 
    t12_blnch = -0.693/real(polos(2));
    omega_esp = w_natural(1); 
    chi_esp = amort(1); 
    t12_esp = -0.693/real(polos(1));


    % Salida de la función 
        % Funciones de transferencia no factorizadas dimensionales
    lat.FT_nfact = struct( ...
    'FTdeltaA_beta',TFdeltaA_betaDim, ...
    'FTdeltaA_phi',TFdeltaA_phiDim, ...
    'FTdeltaA_r',TFdeltaA_rDim, ...
    'FTdeltaR_beta',TFdeltaR_betaDim, ...
    'FTdeltaR_phi',TFdeltaR_phiDim, ...
    'FTdeltaR_r',TFdeltaR_rDim);
        
        % Funciones de tranferencia factorizadas dimensionales
    lat.FT_fact = struct( ...
    'FTdeltaA_beta',zpk(zeros_deltaA_beta,polos_deltaA_beta,ganancia_deltaA_beta), ...
    'FTdeltaA_phi',zpk(zeros_deltaA_phiDim,polos_deltaA_phiDim,ganancia_deltaA_phiDim), ...
    'FTdeltaA_r',zpk(zeros_deltaA_rDim,polos_deltaA_rDim,ganancia_deltaA_rDim), ...
    'FTdeltaR_beta',zpk(zeros_deltaR_betaDim,polos_deltaR_betaDim,ganancia_deltaR_betaDim), ...
    'FTdeltaR_phi',zpk(zeros_deltaR_phiDim,polos_deltaR_phiDim,ganancia_deltaR_phiDim), ...
    'FTdeltaR_r',zpk(zeros_deltaR_rDim,polos_deltaR_rDim,ganancia_deltaR_rDim));

        % Propiedades de los modos 
    lat.Prop_mod = struct(...
        'Omega_Bal',omega_blnc,...
        'Amort_Bal',chi_blnc,...
        'T12_Bal',t12_blnc,...
        'Omega_BalH',omega_blnch,...
        'Amort_BalH',chi_blnch,...
        'T12_BalH',t12_blnch,...
        'Omega_Esp',omega_esp,...
        'Amort_Esp',chi_esp,...
        'T12_Esp',t12_esp);

end