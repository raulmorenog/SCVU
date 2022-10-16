%% Obtención de las FT para el Beech 99
clc
clear all 
close all 

%% Parámetros --> Los guardamos en una estrucutra p
    % Cambios de unidades
ft2m = 0.3048;      % ft a m
lb2kg = 0.4536;     % lb a kg
slg2kg = 14.5939;   % slg a kg 
deg2rad = pi/180;   % º a rad

    % Condición de vuelo
p.Us = 450*ft2m;
p.Ms = 0.434;
p.qs = 128.2*lb2kg/(ft2m^2)*9.81;
p.rhos = 2*p.qs/(p.Us^2);
p.aplhas = 0*deg2rad; 
p.epsilons = 0*deg2rad;
p.thetas = 0*deg2rad;       % Crucero
g = 9.81;

p.Cls = 0.306; 
p.Cds = 0.0298; 
p.CTxs = 0.0298; 
p.Cms = 0; 
p.CmTs = 0; 

    % Geometría
p.Sw = 280*ft2m^2;
p.b = 46*ft2m;
p.c = 6.5*ft2m;
p.cdg_c = 0.16;     % Localización del CDG fracción de la cuerda

    % Parámetros másicos
p.m = 11e3*lb2kg;
p.I_xx = 15.189e3*slg2kg*ft2m^2 ;
p.I_yy = 20.250e3*slg2kg*ft2m^2;
p.I_zz = 34.141e3*slg2kg*ft2m^2;
p.I_xz = 4.371e3*slg2kg*ft2m^2;

p.Czs = -p.m*g*cos(p.thetas)/(p.qs*p.Sw); 
p.Cxs = p.Czs*tan(p.thetas);

    % Derivadas de estabilidad longitudinales
p.Cd_0 = 0.0270; p.Cd_u = 0; p.Cd_alpha = 0.131; 
p.CT_xu = -0.0596;
p.Cl_0 = 0.201; p.Cl_u = 0.020; p.Cl_alpha = 5.48; p.Cl_alphap = 2.5; p.Cl_q = 8.1;
p.Cm_0 = 0.05; p.Cm_u = 0; p.Cm_alpha = -1.89; p.Cm_alphap = -9.1; p.Cm_q = -34;
p.Cm_Tu = 0; p.Cm_Taplha = 0;
p.Cd_deltae = 0; p.Cl_deltae = 0.60; p.Cm_deltae = -2; p.Cm_deltaep = 0;
     

    % Derivadas de estabilidad lateral-direccionales 
p.Cl_beta = -0.130; p.Cl_p = -0.50; p.Cl_r = 0.14; 
p.Cy_beta = -0.59; p.Cy_p = -0.19; p.Cy_r = 0.39;
p.Cn_beta = 0.08; p.Cn_p = -0.019; p.Cn_r = -0.197; 
p.Cn_Tbeta = 0; 
p.Cl_deltaA = 0.156; p.Cl_deltaR = -0.0106; 
p.Cy_deltaA = 0; p.Cy_deltaR = -0.144; 
p.Cn_deltaA = -0.0012 ; p.Cn_deltaR = 0.0758;

%% Adimensionalización longitudinal y coeficientes
mu_long = p.m/(p.rhos*p.Sw*p.c/2);
I_yynd = p.I_yy/(p.rhos*p.Sw*(p.c/2)^3);

Cx_u = p.CT_xu*cos(p.epsilons)-p.Cd_u;     Cx_alpha = p.Cls-p.Cd_alpha;     Cx_alphap = 0;                                  Cx_deltae = -p.Cd_deltae;                                   Cxs = p.Cxs;                             
Cz_u = -p.CT_xu*sin(p.epsilons)-p.Cl_u;    Cz_alpha = -p.Cl_alpha-p.Cds;    Cz_alphap = -p.Cl_alphap;    Cz_q = -p.Cl_q;    Cz_deltae = -p.Cl_deltae;                                   Czs = p.Czs;
Cm_u = p.Cm_u;                             Cm_alpha = p.Cm_alpha;           Cm_alphap = p.Cm_alphap;     Cm_q = p.Cm_q;     Cm_deltae = p.Cm_deltae;     Cm_deltaep = p.Cm_deltaep;        
 
%% Construcción del Sistema Longitudinal: [A]*{Y}=[b]*{X} ----------> Funciones de Transferencia Longitudinales
syms s z
sympref('FloatingPointOutput',true)

A=[2*mu_long*s-Cx_u,-Cx_alpha,-Czs;...
    -(Cz_u+2*Czs),(2*mu_long-Cz_alphap)*s-Cz_alpha,-(2*mu_long+Cz_q)*s;...
    -Cm_u,-(Cm_alphap*s+Cm_alpha),I_yynd*s^2-Cm_q*s];
b=[Cx_deltae;Cz_deltae;Cm_deltaep*s+Cm_deltae];

%Cramer y factorización con s ADIMENSIONAL (s gorro) --> FT con s
%adimensional
Den=sym2poly(det(A)); %Cuártica de estabilidad Longitudinal

Num_deltae_u=sym2poly(det([b,A(:,2),A(:,3)]));
TFdeltae_u_nd=tf([Num_deltae_u],[Den]); %Adimensional!!                                                   
TFdeltae_u= TFdeltae_u_nd*tf([p.Us]); %Dimensional en u, todavia no en s

Num_deltae_alpha=sym2poly(det([A(:,1),b,A(:,3)]));
TFdeltae_alpha=tf([Num_deltae_alpha],[Den]); 

Num_deltae_theta=sym2poly(det([A(:,1),A(:,2),b]));
TFdeltae_theta=tf([Num_deltae_theta],[Den]); 

%Desadimensionalización de s (s = z*c/2Us) --> FT con s dimensional
DenDim = sym2poly(subs(det(A),s,z*p.c/(2*p.Us)));

Num_deltae_uDim = sym2poly(subs(det([b,A(:,2),A(:,3)]),s,z*p.c/(2*p.Us)));
TFdeltae_uDim = tf([Num_deltae_uDim],[DenDim])*tf([p.Us]);

Num_deltae_alphaDim = sym2poly(subs(det([A(:,1),b,A(:,3)]),s,z*p.c/(2*p.Us)));
TFdeltae_alphaDim = tf([Num_deltae_alphaDim],[DenDim]);

Num_deltae_thetaDim = sym2poly(subs(det([A(:,1),A(:,2),b]),s,z*p.c/(2*p.Us)));
TFdeltae_thetaDim = tf([Num_deltae_thetaDim],[DenDim]);

%Paso a forma factorizada
% [zeros_u,polos_u,ganancia_u] = tf2zp(Num_deltae_uDim*p.Us,DenDim)
% [zeros_alpha,polos_alpha,ganancia_alpha] = tf2zp(Num_deltae_alphaDim,DenDim)
% [zeros_theta,polos_theta,ganancia_theta] = tf2zp(Num_deltae_thetaDim,DenDim)

%% Adimensionalización longitudinal y coeficientes
mu_latdir = p.m/(p.rhos*p.Sw*p.b/2);
I_xxnd = p.I_xx/(p.rhos*p.Sw*(p.b/2)^3);
I_zznd = p.I_zz/(p.rhos*p.Sw*(p.b/2)^3);
I_xznd = p.I_xx/(p.rhos*p.Sw*(p.b/2)^3);

Cy_beta = p.Cy_beta;    Cy_p = p.Cy_p;    Cy_r = p.Cy_r;    Cy_deltaA = p.Cy_deltaA;    Cy_deltaR = p.Cy_deltaR;    Czs = p.Czs;
Cl_beta = p.Cl_beta;    Cl_p = p.Cl_p;    Cl_r = p.Cl_r;    Cl_deltaA = p.Cl_deltaA;    Cl_deltaR = p.Cl_deltaR;    Cxs = p.Cxs;
Cn_beta = p.Cn_beta;    Cn_p = p.Cn_p;    Cn_r = p.Cn_r;    Cn_deltaA = p.Cn_deltaA;    Cn_deltaR = p.Cn_deltaR;
%% Construcción del Sistema Lat-Dir: [C]*{Y}=[D]*{X} ----------> Funciones de Transferencia Lat-Dir

syms s z
sympref('FloatingPointOutput',true)

C=[2*mu_latdir*s-Cy_beta,Czs-Cy_p*s,2*mu_latdir-Cy_r;...
    -Cl_beta,I_xxnd*s^2-Cl_p*s,-Cl_r-I_xznd*s;...
    -Cn_beta,-Cn_p*s-I_xznd*s^2,I_zznd*s-Cn_r];
D=[0,Cy_deltaR;Cl_deltaA,Cl_deltaR;Cn_deltaA,Cn_deltaR];

%Cramer y factorización con s ADIMENSIONAL (s gorro) --> FT con s
%adimensional
DenLat=sym2poly(det(C)); %Cuártica de estabilidad Lat-Dir

Num_deltaA_beta = sym2poly(det([D(:,1),C(:,2),C(:,3)]));    Num_deltaR_beta=sym2poly(det([D(:,2),C(:,2),C(:,3)]));
TFdeltaA_beta=tf([Num_deltaA_beta],[DenLat]);             TFdeltaR_beta=tf([Num_deltaR_beta],[DenLat]);                                               

Num_deltaA_phi = sym2poly(det([C(:,1),D(:,1),C(:,3)]));     Num_deltaR_phi=sym2poly(det([C(:,1),D(:,2),C(:,3)]));
TFdeltaA_phi = tf([Num_deltaA_phi],[DenLat]);               TFdeltaR_phi=tf([Num_deltaR_phi],[DenLat]);

Num_deltaA_r=sym2poly(det([C(:,1),C(:,2),D(:,1)]));       Num_deltaR_r=sym2poly(det([C(:,1),C(:,2),D(:,2)]));
TFdeltaA_r_nd=tf([Num_deltaA_r],[DenLat]);                TFdeltaR_r_nd=tf([Num_deltaR_r],[DenLat]);
TFdeltaA_r=TFdeltaA_r_nd*tf(2*p.Us/p.b);                  TFdeltaR_r=TFdeltaR_r_nd*tf(2*p.Us/p.b);

% Desadimensionalización de s (s = z*b/2Us) --> FT con s dimensional
DenLatDim = sym2poly(subs(det(C),s,z*p.b/(2*p.Us)));

Num_deltaA_betaDim = sym2poly(subs(det([D(:,1),C(:,2),C(:,3)]),s,z*p.b/(2*p.Us)));      Num_deltaR_betaDim = sym2poly(subs(det([D(:,2),C(:,2),C(:,3)]),s,z*p.b/(2*p.Us)));
TFdeltaA_betaDim = tf([Num_deltaA_betaDim],[DenLatDim]);                                TFdeltaR_betaDim = tf([Num_deltaR_betaDim],[DenLatDim]);

Num_deltaA_phiDim = sym2poly(subs(det([C(:,1),D(:,1),C(:,3)]),s,z*p.b/(2*p.Us)));       Num_deltaR_phiDim = sym2poly(subs(det([C(:,1),D(:,2),C(:,3)]),s,z*p.b/(2*p.Us)));
TFdeltaA_phiDim = tf([Num_deltaA_phiDim],[DenLatDim]);                                  TFdeltaR_phiDim = tf([Num_deltaR_phiDim],[DenLatDim]);
    
Num_deltaA_rDim = sym2poly(subs(det([C(:,1),C(:,2),D(:,1)]),s,z*p.b/(2*p.Us)));         Num_deltaR_rDim = sym2poly(subs(det([C(:,1),C(:,2),D(:,2)]),s,z*p.b/(2*p.Us)));
TFdeltaA_rDim = tf([Num_deltaA_rDim],[DenLatDim])*tf(2*p.Us/p.b);                       TFdeltaR_rDim = tf([Num_deltaR_rDim],[DenLatDim])*tf(2*p.Us/p.b);


% Agruparmos las FT en estructuras --> Conclusión de momento: Las TF_Long coinciden con las del libro y las de Roskam, las TF_Lat-Dirc coinciden con
% las del libro pero no con Roskam porque tiene signos cambiados en algunas derivadas (Cl_deltaR, Cn_deltaR, Cy_deltaR, Cn_p). Falta factorizar guay, diagramas y respuesta

FT_Longitudinales = struct( ...
    'FTdeltaE_u',TFdeltae_uDim, ...
    'FTdeltaE_alpha',TFdeltae_alphaDim, ...
    'FTdeltaE_theta',TFdeltae_thetaDim);

FT_Lateral_Direccionales = struct( ...
    'FTdeltaA_beta',TFdeltaA_betaDim, ...
    'FTdeltaA_phi',TFdeltaA_phiDim, ...
    'FTdeltaA_r',TFdeltaA_rDim, ...
    'FTdeltaR_beta',TFdeltaR_betaDim, ...
    'FTdeltaR_phi',TFdeltaR_phiDim, ...
    'FTdeltaR_r',TFdeltaR_rDim);


%% Obtención de las propiedades de las funciones de transferencia
    %Paso a forma factorizada
[zeros_u,polos_u,ganancia_u] = tf2zp(Num_deltae_uDim*p.Us,DenDim)
G = zpk(zeros_u,polos_u,ganancia_u)         
[w_natural,amort,polos] = damp(TFdeltae_uDim)

[zeros_u,polos_u,ganancia_u] = tf2zp(Num_deltaR_betaDim ,DenLatDim)
G = zpk(zeros_u,polos_u,ganancia_u)         
[w_natural,amort,polos] = damp(TFdeltaR_betaDim)

%% Prueba de las funciones

LO = sist_long(p)
LT = sist_lat(p)


