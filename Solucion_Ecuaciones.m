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
%p.Cd_ih = 0; p.Cl_ih = 0; p.Ci_deltah = 0;     

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

Cx_u = p.CT_xu*cos(p.epsilons)-p.Cd_u; 
Cz_u = -p.CT_xu*sin(p.epsilons)-p.Cl_u; 
Cm_u = p.Cm_u;
Cx_alpha = p.Cls-p.Cd_alpha;
Cz_alpha = -p.Cl_alpha-p.Cds; 
Cm_alpha = p.Cm_alpha;
Cx_alphap = 0;
Cz_alphap = -p.Cl_alphap;
Cm_alphap = p.Cm_alphap;
Cz_q = -p.Cl_q; 
Cm_q = p.Cm_q;
Cx_deltae = -p.Cd_deltae; 
Cz_deltae = -p.Cl_deltae; 
Cm_deltae = p.Cm_deltae;
Cm_deltaep = p.Cm_deltaep;

Czs = p.Czs;
Cxs = p.Cxs;
%% Construcción del Sistema Longitudinal: [A]*{Y}=[b]*{X} ----------> Funciones de Transferencia Longitudinales
syms s
sympref('FloatingPointOutput',true)

A=[2*mu_long*s-Cx_u,-Cx_alpha,-Czs;...
    -(Cz_u+2*Czs),(2*mu_long-Cz_alphap)*s-Cz_alpha,-(2*mu_long+Cz_q)*s;...
    -Cm_u,-(Cm_alphap*s+Cm_alpha),I_yynd*s^2-Cm_q*s];
b=[Cx_deltae;Cz_deltae;Cm_deltaep*s+Cm_deltae];

%Cramer y factorización
Den=sym2poly(det(A));

Num_deltae_u=sym2poly(det([b,A(:,2),A(:,3)]));
TFdeltae_u_nd=tf([Num_deltae_u],[Den]) %Adimensional!!
TFdeltae_u= TFdeltae_u_nd*tf([p.Us])
[zeros_u,polos_u,ganancia_u] = tf2zp(Num_deltae_u*p.Us,Den)


Num_deltae_alpha=sym2poly(det([A(:,1),b,A(:,3)]));
TFdeltae_alpha=tf([Num_deltae_alpha],[Den]) 
[zeros_alpha,polos_alpha,ganancia_alpha] = tf2zp(Num_deltae_alpha,Den)

Num_deltae_theta=sym2poly(det([A(:,1),A(:,2),b]));
TFdeltae_theta=tf([Num_deltae_theta],[Den]) 
[zeros_theta,polos_theta,ganancia_theta] = tf2zp(Num_deltae_theta,Den)

