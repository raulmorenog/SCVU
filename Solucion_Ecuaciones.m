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
p.Cn_beta = 0.08; p.Cn_p = 0.019; p.Cn_r = -0.197;  %Cn_p puesto positivo aunque en los datos era negativo.
p.Cn_Tbeta = 0; 
p.Cl_deltaA = 0.156; p.Cl_deltaR = -0.0106; 
p.Cy_deltaA = 0; p.Cy_deltaR = -0.144; 
p.Cn_deltaA = -0.0012 ; p.Cn_deltaR = 0.0758;

%% Funciones de Transferencia
LO = sist_long(p);
LT = sist_lat(p);
    % Las guardamos en un cell para trabajar más fácil con ellas (en bucles)
FT_long{1} = LO.FT_fact.FTdeltaE_u;
FT_long{2} = LO.FT_fact.FTdeltaE_alpha;
FT_long{3} = LO.FT_fact.FTdeltaE_theta;

FT_lat{1} = LT.FT_fact.FTdeltaA_beta; 
FT_lat{2} = LT.FT_fact.FTdeltaA_phi;
FT_lat{3} = LT.FT_fact.FTdeltaA_r;
FT_lat{4} = LT.FT_fact.FTdeltaR_beta; 
FT_lat{5} = LT.FT_fact.FTdeltaR_phi;
FT_lat{6} = LT.FT_fact.FTdeltaR_r;

    % Mostramos las FT
disp('Funciones de transferencia longitudinales, [u, alpha, theta])')
for i = 1:length(FT_long)
    FT_long{i}
end
disp('Funciones de transferencia lateral-direccional, [beta, phi, r][deltaA, deltaR])')
for i = 1:length(FT_lat)
    FT_lat{i}
end

%% Diagramas de Bode 
label_long = {'$G_{u \delta_{E}}$','$G_{\alpha \delta_{E}}$','$G_{\theta \delta_{E}}$'};
label_lat = {'$G_{\beta \delta_{A}}$','$G_{\phi \delta_{A}}$','$G_{r \delta_{A}}$',...
    '$G_{\beta \delta_{R}}$','$G_{\phi \delta_{R}}$','$G_{r \delta_{R}}$'};

for i = 1:length(FT_long)
    figure(1);
    bode(FT_long{i})
    grid on; hold all 
%     margin(FT_long{i})
%     hold all
    set(gcf,'Color',[1 1 1])
end
legend(label_long,'Interpreter','latex','Location','best')
for i = 1:length(FT_lat)
    figure(2);
    bode(FT_lat{i})
    grid on; hold all;
%     margin(FT_lat{i})
%     hold all
    set(gcf,'Color',[1 1 1])
end
legend(label_lat,'Interpreter','latex','Location','best')

%% Diagramas de Nichols

for i = 1:length(FT_long)
    figure(3);
    nichols(FT_long{i})
    grid on; hold all 
    set(gcf,'Color',[1 1 1])
end
legend(label_long,'Interpreter','latex','Location','best')
for i = 1:length(FT_lat)
    figure(4);
    nichols(FT_lat{i})
    grid on; hold all;
    set(gcf,'Color',[1 1 1])
end
legend(label_lat,'Interpreter','latex','Location','best')

%% Respuesta escalón 
for i = 1:length(FT_long)
    figure(5);
    step(FT_long{i})
    grid on; hold all 
    set(gcf,'Color',[1 1 1])
end
legend(label_long,'Interpreter','latex','Location','best')
for i = 1:length(FT_lat)
    figure(6);
    step(FT_lat{i})
    grid on; hold all;
    set(gcf,'Color',[1 1 1])
end
legend(label_lat,'Interpreter','latex','Location','best')

