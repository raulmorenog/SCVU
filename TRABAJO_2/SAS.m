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
Lc = 4.5e-3; Rc = 0.64; tauc = Lc/Rc; %Propiedades eléctricas
Kmv = 0.0426; Jm= 3.36e-3; taum = Jm/Kmv; %Propiedades mecánicas
Ge = tf([1/tauc],[1,1/tauc]); %Delay eléctrico ~ 7ms
Gm = tf([1/taum],[1,1/taum]); %Delay mecánico ~ 80ms
G_act = tf([1],[tauc+taum,1]); %Delay eléctrico + mecánico // Equivalente a Ge*Gm

figure(1); bode(G_act); 
figure(2); step(G_act);
%Opción B: Delay puro de wn = 20rad/s (Literatura)
%H = tf([1],[1/20,1])

%% Sensores
%Sensor de beta --> Veleta
Delay_dist = 0.65; %m
V_wind_tunnel = 10; %m/s
delay_vane = Delay_dist/V_wind_tunnel; %s

[num_vane,den_vane] = pade(delay_vane,2); %Aproximacion de Pade de orden 2
G_vane = tf(num_vane,den_vane);
%G_vane = tf([-delay_vane/2,1],[delay_vane/2,1]); %Forma menos elegante

figure(3); bode(G_vane); 
figure(4); step(G_vane);

%Sensor de r --> Giróscopo (IMU)
delay_gyro = 10e-3; %s
[num_gyro,den_gyro] = pade(delay_gyro,2); %Aproximacion de Pade de orden 2
G_gyro = tf(num_gyro,den_gyro);
%G_gyro = tf([-delay_gyro/2,1],[delay_gyro/2,1]); %Forma menos elegante

figure(5); bode(G_gyro); 
figure(6); step(G_gyro);

%% Análisis de sensibilidad
X_p = [real(FT_lat.Poles)]; %Parte real de los polos del la planta libre
Y_p = [imag(FT_lat.Poles)]; %Parte imaginaria de los polos de la planta libre
F = 0:0.5:5;    %Valores para los cuales se va a efectuar el barrido
% Variamos los coeficientes
Cn_beta = p.Cn_beta*F;
Cn_r = p.Cn_r*F;
POLES = 1e20;
X = 1e20; Y = 1e20;
for i = 1:length(F)
    for j = 1:length(F)
        s = p;
        s.Cn_beta = Cn_beta(i);
        s.Cn_r = Cn_r(j);
        FT_sensibilidad(i,j) = FT_lat_function_elegante(s);
        %Representamos los polos para el barrido de coeficientes
        X = [X;real(FT_sensibilidad(i,j).dutchroll.poles)];
        Y = [Y;imag(FT_sensibilidad(i,j).dutchroll.poles)];
        POLES = [POLES;FT_sensibilidad(i,j).Poles];
    end
end
Contorno_raices = figure(7); %Entregable
scatter(X,Y,'x'); 
hold on;
scatter(X_p,Y_p,'*','k');grid on; axis([-2.5 0.2 -4.5 4.5]);
xlabel('$$Re [s^{-1}]$$','interpreter','latex'); ylabel('$$Im [s^{-1}]$$','interpreter','latex');
wnDR_lim = 1;   %Límite de normas
chiDR_lim = 0.19;    %Límite de normas
%Representamos los límites de las normas para elegir los valores target
%deseados
A = -5:0.1:0.5;
figure(8)
scatter(X,Y,'x'); grid on; axis([-2 1 0 5]);
hold on
plot(A,-tan(acos(chiDR_lim)).*A,'k-')
viscircles([0 0],wnDR_lim)

s.Cn_beta = Cn_beta(3);    %Varlor deseado de Cn_beta
s.Cn_r = Cn_r(11);         %Valor deseado de Cn_r_target
FT_22 = FT_lat_function_elegante(s);
X_m = [real(FT_22.Poles)];
Y_m = [imag(FT_22.Poles)];
scatter(X_m,Y_m,'o','r')
scatter(X_p,Y_p,'*','k')

%Elección de las derivadas de estabilidad target y nuevas características
%del modo
Cn_beta_target = Cn_beta(3);
Cn_r_target = Cn_r(11);
wnDR_target = FT_22.dutchroll.wn;
chiDR_target = FT_22.dutchroll.amort;

%% Direct link
% Nos basamos en el modelo de 1gdl como indican las diapos, calculado
% analíticamente
FT_conv_balance_1gdl.nofact = tf([p.Cl_deltaA*2*p.Us/p.b],[p.I_xx/(p.rhos*p.Sw*(p.b/2)^3)*(p.b/2/p.Us), -p.Cl_p]);
%Dimensionalizamos
FT_conv_balance_1gdl.fact = zpk(FT_conv_balance_1gdl.nofact);
FT_conv_balance_1gdl.Kstatic = -FT_conv_balance_1gdl.fact.K/FT_conv_balance_1gdl.fact.P{1,1};   %SAco la ganancia estática
K_DL = 1/FT_conv_balance_1gdl.Kstatic;  %Valor de la ganancia del direct link

%% Respuesta temporal con el direct link y actuador (compilar este apartado cuando el Simulink esté correcto)
%Habrá que definir la duración del escalón, es de 20s que lo hemos sacado a
%ojo de las gráficas de respuesta a escalón.
RT_OL_133 = sim('KDL_act_planta_133',40);     %Llamo a modelo de simulink para apartado 1.3.3      
figure(9)
plot(RT_OL_133.tout,RT_OL_133.deltaA_stick); grid on;   %Deflexión elegida para el stick
xlabel('$$t \mathrm{[s]}$$','interpreter','latex','FontSize',14)
ylabel('$$\delta_{a,stick} \mathrm{[^o]}$$','interpreter','latex','FontSize',14)
figure(10)
plot(RT_OL_133.tout,RT_OL_133.beta); grid on;      %Ángulo de resbalamiento
xlabel('$$t \mathrm{[s]}$$','interpreter','latex','FontSize',14)
ylabel('$$\beta \mathrm{[^o]}$$','interpreter','latex','FontSize',14)
figure(11)
plot(RT_OL_133.tout,RT_OL_133.phi); grid on;   %Ángulo de balance
xlabel('$$t \mathrm{[s]}$$','interpreter','latex','FontSize',14)
ylabel('$$\phi \mathrm{[^o]}$$','interpreter','latex','FontSize',14)
figure(12)
plot(RT_OL_133.tout,RT_OL_133.r); grid on;   %Velocidad de guiñada
xlabel('$$t \mathrm{[s]}$$','interpreter','latex','FontSize',14)
ylabel('$$r \mathrm{[^o]}$$','interpreter','latex','FontSize',14)
figure(13)
plot(RT_OL_133.tout,RT_OL_133.p);  grid on;   %Velocidad de balance
xlabel('$$t \mathrm{[s]}$$','interpreter','latex','FontSize',14)
ylabel('$$p \mathrm{[^o]}$$','interpreter','latex','FontSize',14)

%% SAS 
