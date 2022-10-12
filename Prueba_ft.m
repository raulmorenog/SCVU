%% Script para probar el funcinamiento de las funciones que obtienen FT
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
p.Us = 677*ft2m;
p.Ms = 0.7;
p.qs = 134.6*lb2kg/(ft2m^2);
p.rhos = 2*p.qs/(p.Us^2);
p.aplhas = 2.7*deg2rad; 
p.epsilons = 0*deg2rad;
p.thetas = 0*deg2rad;       
g = 9.81;

p.Cls = 0.41; 
p.Cds = 0.0335; 
p.CTxs = 0.0335; 
p.Cms = 0; 
p.CmTs = 0; 

    % Geometría
p.Sw = 230*ft2m^2;
p.b = 34*ft2m;
p.c = 7*ft2m;

    % Parámetros másicos
p.m = 13e3*lb2kg;
p.I_xx = 28e3*slg2kg*ft2m^2 ;
p.I_yy = 18.8e3*slg2kg*ft2m^2;
p.I_zz = 47e3*slg2kg*ft2m^2;
p.I_xz = 1.3e3*slg2kg*ft2m^2;

p.Czs = p.m*g*cos(p.thetas)/(p.qs*p.Sw);
p.Cxs = p.Czs*tan(p.thetas);

    % Derivadas de estabilidad longitudinales
p.Cd_0 = 0.0216; p.Cd_u = 0.104; p.Cd_alpha = 0.3; 
p.CT_xu = -0.07;
p.Cl_0 = 013; Cp.l_u = 0.4; p.Cl_alpha = 5.84; p.Cl_alphap = 2.2; p.Cl_q = 4.7;
p.Cm_0 = 0.05; p.Cm_u = 0.05; p.Cm_alpha = -0.64; p.Cm_alphap = -6.7; p.Cm_q = -15.5;
p.Cm_Tu = -0.003; p.Cm_Taplha = 0;
p.Cd_deltae = 0; p.Cl_deltae = 0.46; p.Cm_deltae = -1.24; p.Cm_deltaep = 0;
p.Cd_ih = 0; p.Cl_ih = 0.94; p.Ci_deltah = -2.5;

    % Derivadas de estabilidad lateral-direccionales 
p.Cl_beta = -0.11; p.Cl_p = -0.45; p.Cl_r = 0.16; 
p.Cy_beta = -0.73; p.Cy_p = 0; p.Cy_r = 0.4;
p.Cn_beta = 0.127; p.Cn_beta = -0.008; p.Cn_r = -0.2; 
p.Cn_Tbeta = 0; 
p.Cl_deltaA = 0.178; p.Cl_deltaR = -0.019; 
p.Cy_deltaA = 0; p.Cy_deltaR = -0.14; 
p.Cn_deltaA = -0.02 ; p.Cn_deltaR = 0.074;

%% Obtenemos las funciones de transferencia 
    % Movimiento Longitudinal
long = longitudinal(p);
long.coef_cuartica; 

[long.fugoide.w_fug long.fugoide.chi_fug];
[long.cortoperiodo.w_cp long.cortoperiodo.chi_cp];
    % Funciones de transferencia
long.Gu_deltae
long.Galpha_deltae
long.Gtheta_deltae

