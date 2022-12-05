%% SAS_LATERAL
clc
clear all 
close all 

%% Metemos todos los datos
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
p.alphas = 0*deg2rad; 
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
p.Cxs = -p.Czs*tan(p.thetas);

% Derivadas de estabilidad lateral-direccionales 
p.Cl_beta = -0.130; p.Cl_p = -0.50; p.Cl_r = 0.14; 
p.Cy_beta = -0.59; p.Cy_p = -0.19; p.Cy_r = 0.39;
p.Cn_beta = 0.08; p.Cn_p = 0.019; p.Cn_r = -0.197;  %Cn_p puesto positivo aunque en los datos era negativo.
p.Cn_Tbeta = 0; 
p.Cl_deltaA = 0.156; p.Cl_deltaR = -0.0106; 
p.Cy_deltaA = 0; p.Cy_deltaR = -0.144; 
p.Cn_deltaA = -0.0012 ; p.Cn_deltaR = 0.0758;

%% FT de la planta totalmente libre
FT_lat = FT_lat_function_elegante(p);   %Calcularmos las funciones de transferencia de la planta sin aumentar

%% Actuadores EMA (Electro-Mechanical Actuator)
Lc = 4.5e-3; Rc = 0.64; tauc = Lc/Rc;       % Propiedades eléctricas
Kmv = 0.0426; Jm= 3.36e-3; taum = Jm/Kmv;   % Propiedades mecánicas
Ge = tf([1/tauc],[1,1/tauc]);   % Delay eléctrico ~ 7ms
Gm = tf([1/taum],[1,1/taum]);   % Delay mecánico ~ 80ms
G_act = tf([1],[tauc+taum,1]);  %Delay eléctrico + mecánico // Equivalente a Ge*Gm

figure(1); bode(G_act); 
figure(2); step(G_act);
% Opción B: Delay puro de wn = 20rad/s (Literatura)
% H = tf([1],[1/20,1])

%% Sensores
    % Sensor de beta --> Veleta
Delay_dist = 0.65;      % m
V_wind_tunnel = 10;     % m/s
delay_vane = Delay_dist/V_wind_tunnel;  % s

[num_vane,den_vane] = pade(delay_vane,2); % Aproximacion de Pade de orden 2
G_vane = tf(num_vane,den_vane);
% G_vane = tf([-delay_vane/2,1],[delay_vane/2,1]); % Forma menos elegante

figure(3); bode(G_vane); 
figure(4); step(G_vane);

    % Sensor de r --> Giróscopo (IMU)
delay_gyro = 10e-3; %s
[num_gyro,den_gyro] = pade(delay_gyro,2); % Aproximacion de Pade de orden 2
G_gyro = tf(num_gyro,den_gyro);
% G_gyro = tf([-delay_gyro/2,1],[delay_gyro/2,1]); % Forma menos elegante

figure(5); bode(G_gyro); 
figure(6); step(G_gyro);

%% Análisis de sensibilidad
X_p = [real(FT_lat.Poles)]; % Parte real de los polos del la planta libre
Y_p = [imag(FT_lat.Poles)]; % Parte imaginaria de los polos de la planta libre
F = 0:0.5:5;    % Valores para los cuales se va a efectuar el barrido
% Variamos los coeficientes
Cn_beta = p.Cn_beta*F;
Cn_r = p.Cn_r*F;
POLES = 1e20;
X_dr = []; Y_dr = [];
X_s = []; Y_s = [];
X_r = []; Y_r = [];
for i = 1:length(F)
    for j = 1:length(F)
        s = p;
        s.Cn_beta = Cn_beta(i);
        s.Cn_r = Cn_r(j);
        FT_sensibilidad(i,j) = FT_lat_function_elegante(s);
        % Representamos los polos para el barrido de coeficientes
        X_dr = [X_dr;real(FT_sensibilidad(i,j).dutchroll.poles)];
        Y_dr = [Y_dr;imag(FT_sensibilidad(i,j).dutchroll.poles)];
        X_s = [X_s;real(FT_sensibilidad(i,j).espiral.poles)];
        Y_s = [Y_s;imag(FT_sensibilidad(i,j).espiral.poles)];
        X_r = [X_r;real(FT_sensibilidad(i,j).balance.poles)];
        Y_r = [Y_r;imag(FT_sensibilidad(i,j).balance.poles)];
        POLES = [POLES;FT_sensibilidad(i,j).Poles];
    end
end
Contorno_raices = figure(7);    % Entregable
p1 = plot(X_dr,Y_dr,'dm','markersize',4,'markerfacecolor','m'); hold on;
p2 = plot(X_s,Y_s,'dc','markersize',4,'markerfacecolor','c'); hold on;
p3 = plot(X_r,Y_r,'dg','markersize',4,'markerfacecolor','g'); hold on;
grid on;
axis([-4 0.5 -5.5 5.5]);
xa = [.87 .87]; ya = [.65 .9]; 
annotation('arrow',xa,ya,'color','k'); hold on;
xa = [.8 .6]; ya = [.9 .9]; 
annotation('arrow',xa,ya,'color','k'); hold on;
text1 = '$$C_{n\beta}$$'; text2 = '$$C_{nr}$$'; 
text(-0.06,4.8,text1,'fontsize',14,'interpreter','latex'); hold on; 
text(-0.75,4.87,text2,'fontsize',14,'interpreter','latex');

xlabel('$$Re$$ $$[\mathrm{s^{-1}}]$$','interpreter','latex','fontsize',14); 
ylabel('$$Im$$ $$[\mathrm{s^{-1}}]$$','interpreter','latex','fontsize',14);
legend([p1 p2 p3],{'Balanceo Holandes','Espiral','Convergencia balance'},...
    'location', 'northwest', 'orientation','vertical','interpreter','latex',...
    'fontsize',14);

% Representamos los límites de las normas para elegir los valores target
% deseados
A = -5:0.1:0.5;
wnDR_lim = 1;           % Límite de normas
chiDR_lim = 0.19;       % Límite de normas
figure(8)
plot(X_dr,Y_dr,'dm','markersize',4,'markerfacecolor','m'); hold on;
grid on; 
axis([-2 1 0 5]);
plot(A,-tan(acos(chiDR_lim)).*A,'k-','linewidth',1); hold on;
viscircles([0 0],wnDR_lim,'linewidth',1,'color','k'); hold on;
xline(-0.35,'k-','linewidth',1); hold on;

s.Cn_beta = Cn_beta(3);    % Valor deseado de Cn_beta
s.Cn_r = Cn_r(11);         % Valor deseado de Cn_r_target
FT_22 = FT_lat_function_elegante(s);
X_m = [real(FT_22.Poles)];
Y_m = [imag(FT_22.Poles)];
p1 = plot(X_m(2),Y_m(2),'pr','markersize',10,'markerfacecolor','r'); hold on;
p2 = plot(X_p(2),Y_p(2),'ob','markersize',8,'markerfacecolor','b'); hold on;
xlabel('$$Re$$ $$[\mathrm{s^{-1}}]$$','interpreter','latex','fontsize',14); 
ylabel('$$Im$$ $$[\mathrm{s^{-1}}]$$','interpreter','latex','fontsize',14);
legend([p1 p2],{'Punto objetivo','Planta libre'},'location', 'northeast',...
    'orientation','vertical','interpreter','latex','fontsize',14)

% Elección de las derivadas de estabilidad target y nuevas características
% del modo
Cn_beta_target = Cn_beta(3);
Cn_r_target = Cn_r(11);
wnDR_target = FT_22.dutchroll.wn;
chiDR_target = FT_22.dutchroll.amort;

%% Direct link
% Nos basamos en el modelo de 1gdl como indican las diapos, calculado
% analíticamente
FT_conv_balance_1gdl.nofact = tf([p.Cl_deltaA*2*p.Us/p.b],[p.I_xx/(p.rhos*p.Sw*(p.b/2)^3)*(p.b/2/p.Us), -p.Cl_p]);
% Dimensionalizamos
FT_conv_balance_1gdl.fact = zpk(FT_conv_balance_1gdl.nofact);
FT_conv_balance_1gdl.Kstatic = -FT_conv_balance_1gdl.fact.K/FT_conv_balance_1gdl.fact.P{1,1};   %SAco la ganancia estática
K_DL = 1/FT_conv_balance_1gdl.Kstatic;  % Valor de la ganancia del direct link

%% Respuesta temporal con el direct link y actuador (compilar este apartado cuando el Simulink esté correcto)
% Habrá que definir la duración del escalón, es de 20s que lo hemos sacado a
% ojo de las gráficas de respuesta a escalón.

% Para la respuesta teniendo en cuenta el sistema de 3gdl
RT_OL_133 = sim('KDL_act_planta_133',40);     % Llamo a modelo de simulink para apartado 1.3.3      

figure(9)
plot(RT_OL_133.tout,RT_OL_133.deltaA_stick); grid on;   % Deflexión elegida para el stick
xlabel('$$t \mathrm{[s]}$$','interpreter','latex','FontSize',14)
ylabel('$$\delta_{a,stick} \mathrm{[^o]}$$','interpreter','latex','FontSize',14)

figure(10)
subplot(2,2,1)
plot(RT_OL_133.tout,RT_OL_133.beta); grid on;  % Ángulo de resbalamiento
xticks(0:10:40);
yticks(-0.1:0.05:0.2);
axis([0 40 -0.05 0.2]);
xlabel('$$t \mathrm{[s]}$$','interpreter','latex','FontSize',14)
ylabel('$$\beta \mathrm{[^o]}$$','interpreter','latex','FontSize',14)

subplot(2,2,2)
plot(RT_OL_133.tout,RT_OL_133.phi); grid on;   % Ángulo de balance
xticks(0:10:40);
yticks(0:3:15);
axis([0 40 0 16]);
xlabel('$$t \mathrm{[s]}$$','interpreter','latex','FontSize',14)
ylabel('$$\phi \mathrm{[^o]}$$','interpreter','latex','FontSize',14)

subplot(2,2,3)
plot(RT_OL_133.tout,RT_OL_133.r); grid on;      % Velocidad de guiñada
xticks(0:10:40);
yticks(0:0.2:1.2);
axis([0 40 -0.1 1.2]);
xlabel('$$t \mathrm{[s]}$$','interpreter','latex','FontSize',14)
ylabel('$$r \mathrm{[^o/s]}$$','interpreter','latex','FontSize',14)

subplot(2,2,4)
plot(RT_OL_133.tout,RT_OL_133.p);  grid on;     % Velocidad de balance
xticks(0:10:40);
yticks(-0.5:0.5:1.5);
axis([0 40 -0.5 1.5]);
xlabel('$$t \mathrm{[s]}$$','interpreter','latex','FontSize',14)
ylabel('$$p \mathrm{[^o/s]}$$','interpreter','latex','FontSize',14)
sgtitle('Variables de estado para modelo 3 gdl','interpreter','latex',...
    'fontsize',14)

    % Respuesta con el modelo de 1gdl
RT_OL_133_simp = sim('modelo_convergencia_balance_1gdl',40);
figure(11)
plot(RT_OL_133_simp.tout,RT_OL_133_simp.deltaA_stick); grid on;   % Deflexión elegida para el stick
xlabel('$$t \mathrm{[s]}$$','interpreter','latex','FontSize',14)
ylabel('$$\delta_{a,stick} \mathrm{[^o]}$$','interpreter','latex','FontSize',14)

figure(12)
subplot(1,2,1)
plot(RT_OL_133_simp.tout,RT_OL_133_simp.phi); grid on;  % Ángulo de balance
xticks(0:5:40);
yticks(0:3:15);
axis([0 40 0 16]);
xlabel('$$t \mathrm{[s]}$$','interpreter','latex','FontSize',14)
ylabel('$$\phi \mathrm{[^o]}$$','interpreter','latex','FontSize',14)

subplot(1,2,2)
plot(RT_OL_133_simp.tout,RT_OL_133_simp.p); grid on;   % Velocidad de balance
xticks(0:5:40);
yticks(-0.5:0.25:1.5);
axis([0 40 -0.5 1.5]);
xlabel('$$t \mathrm{[s]}$$','interpreter','latex','FontSize',14)
ylabel('$$p \mathrm{[^o/s]}$$','interpreter','latex','FontSize',14)
sgtitle('Variables de estado para modelo 1 gld','interpreter','latex',...
    'fontsize',14)


%% SAS 
% Calculamos todas las FT a mano en lazo cerrado
k_deltar_beta = 1; k_deltar_r = 1; G_w = 1;     % Variables todavía no definidas
s_lp = tf([1 0],[1]);   % Defino la s como función de transferencia para la derivada en integral necesaria
FT_CL.beta_deltaS = (FT_lat.fact.deltaA_beta*G_act*K_DL)/(1 + ...
    G_act*FT_lat.fact.deltaR_beta*(k_deltar_beta*G_vane + s_lp*k_deltar_r*G_gyro*G_w));
FT_CL.r_deltaS = (FT_lat.fact.deltaA_r*G_act*K_DL)/(1 + ...
    G_act*FT_lat.fact.deltaR_r*(k_deltar_beta*G_vane/s_lp + k_deltar_r*G_gyro*G_w));
FT_CL.phi_deltaS = FT_lat.fact.deltaA_phi*G_act*K_DL - G_act*FT_lat.fact.deltaR_phi*...
    (k_deltar_beta*G_vane*FT_CL.beta_deltaS + k_deltar_r*G_gyro*G_w*FT_CL.r_deltaS);
FT_CL.p_deltaS = FT_lat.fact.deltaA_p*G_act*K_DL - G_act*FT_lat.fact.deltaR_p*...
    (k_deltar_beta*G_vane*FT_CL.beta_deltaS + k_deltar_r*G_gyro*G_w*FT_CL.r_deltaS);

FT_CL.beta_deltaS = minreal(FT_CL.beta_deltaS,0.1);
FT_CL.r_deltaS = minreal(FT_CL.r_deltaS,0.1);
FT_CL.phi_deltaS = minreal(FT_CL.phi_deltaS,0.1);
FT_CL.p_deltaS = minreal(FT_CL.p_deltaS,0.1);

% figure
% zplane(FT_CL.beta_deltaS.Z{1, 1},FT_CL.beta_deltaS.P{1, 1})   
% figure
% zplane(FT_CL.r_deltaS.Z{1, 1},FT_CL.r_deltaS.P{1, 1})  
% figure
% zplane(FT_CL.phi_deltaS.Z{1, 1},FT_CL.phi_deltaS.P{1, 1}) 
% figure
% zplane(FT_CL.p_deltaS.Z{1, 1},FT_CL.p_deltaS.P{1, 1})  


% FT en lazo cerrado Hugo y Raúl
Ga_deltaA = G_act; Ga_deltaR = G_act;       % FT actuadores
Gs_r = G_gyro; Gs_beta = G_vane;            % FT sensores
Gf_r = G_w;                                 % FT filtro wash-out (supuesto 1)
K_deltaRbeta = 1; K_deltaRr = 1;            % Ganancias de realimentación
    % Funciones de transferencia
G_betaDeltaA = FT_lat.fact.deltaA_beta; G_betaDeltaR = FT_lat.fact.deltaR_beta; 
G_rDeltaA = FT_lat.fact.deltaA_r; G_rDeltaR = FT_lat.fact.deltaR_r; 
G_phiDeltaA = FT_lat.fact.deltaA_phi; G_phiDeltaR = FT_lat.fact.deltaR_phi; 
G_pDeltaA = FT_lat.fact.deltaA_p; G_pDeltaR = FT_lat.fact.deltaR_p;

    % Construcción de las FT en lazo cerrado 
M = 1+Ga_deltaR*G_betaDeltaR*K_deltaRbeta*Gs_beta; 
num_rDeltaA = K_DL*Ga_deltaA*G_rDeltaA -...
    Ga_deltaR*Ga_deltaA*K_DL*G_rDeltaR*G_betaDeltaA*Gs_beta*K_deltaRbeta/M;
den_rDeltaA = 1+Ga_deltaR*Gf_r*Gs_r*K_deltaRr*G_rDeltaR*(1+...
    K_deltaRbeta*Gs_beta*Ga_deltaR*G_betaDeltaR/M); 
FT_cl.r_deltaS =  num_rDeltaA/den_rDeltaA;

num_betaDeltaA = K_DL*Ga_deltaA*G_betaDeltaA-...
    Ga_deltaR*G_betaDeltaR*K_deltaRr*Gf_r*Gs_r*FT_cl.r_deltaS;
den_betaDeltaA = M; 
FT_cl.beta_deltaS =  num_betaDeltaA/den_betaDeltaA;

FT_cl.p_deltaS = K_DL*Ga_deltaA*G_pDeltaA-...
    Ga_deltaR*K_deltaRr*Gf_r*Gs_r*G_pDeltaR*FT_cl.r_deltaS-...
    Ga_deltaR*K_deltaRbeta*Gs_beta*G_pDeltaR*FT_cl.beta_deltaS; 

FT_cl.phi_deltaS = K_DL*Ga_deltaA*G_phiDeltaA-...
    Ga_deltaR*K_deltaRr*Gf_r*Gs_r*G_phiDeltaR*FT_cl.r_deltaS-...
    Ga_deltaR*K_deltaRbeta*Gs_beta*G_phiDeltaR*FT_cl.beta_deltaS;


FT_cl.beta_deltaS =  minreal(FT_cl.beta_deltaS,0.1);
FT_cl.r_deltaS =  minreal(FT_cl.r_deltaS,0.1);
FT_cl.phi_deltaS =  minreal(FT_cl.phi_deltaS,0.1);
FT_cl.p_deltaS =  minreal(FT_cl.p_deltaS,0.1);

    % FT en lazo abierto
FT_ol = minreal(den_rDeltaA);

    % Mapas de polos
% figure
% zplane(FT_cl.beta_deltaS.Z{1, 1},FT_cl.beta_deltaS.P{1, 1})
% figure
% zplane(FT_cl.r_deltaS.Z{1, 1},FT_cl.r_deltaS.P{1, 1})
% figure
% zplane(FT_cl.phi_deltaS.Z{1, 1},FT_cl.phi_deltaS.P{1, 1})
% figure
% zplane(FT_cl.p_deltaS.Z{1, 1},FT_cl.p_deltaS.P{1, 1})