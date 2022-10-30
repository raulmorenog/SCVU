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
%FT_long{4} = LO.FT_fact.FTdeltaE_q;

FT_lat{1} = LT.FT_fact.FTdeltaA_beta; 
FT_lat{2} = LT.FT_fact.FTdeltaA_phi;
FT_lat{3} = LT.FT_fact.FTdeltaA_p;
FT_lat{4} = LT.FT_fact.FTdeltaA_r;
FT_lat{5} = LT.FT_fact.FTdeltaR_beta; 
FT_lat{6} = LT.FT_fact.FTdeltaR_phi;
FT_lat{7} = LT.FT_fact.FTdeltaR_r;


    % Mostramos las FT
disp('Funciones de transferencia longitudinales, [u, alpha, theta])')
for i = 1:length(FT_long)
    FT_long{i}
end
disp('Funciones de transferencia lateral-direccional, [beta, phi, p, r][deltaA, deltaR])')
for i = 1:length(FT_lat)
    FT_lat{i}
end

%% Diagramas de Bode 
label_long = {'$$G\_{u\delta\_{e}}$$','$$G\_{\alpha\delta\_{e}}$$','$$G\_{\theta \delta\_{e}}$$'};
label_lat = {'$$G\_{\beta\delta\_{a}}$$','$$G\_{\phi \delta\_{a}}$$','$$G\_{p \delta\_{a}}$$',...
    '$$G\_{r \delta\_{a}}$$','$$G\_{\beta \delta\_{r}}$$','$$G\_{\phi \delta\_{r}}$$','$$G\_{r \delta\_{r}}$$'};

k = 1;
iFig = 1;
for i = 1:length(FT_long)
    figure(iFig);
    bode(FT_long{i},{10^-3,10^4},'b-')
    grid on; hold all 
    margin(FT_long{i})
    set(gcf,'Color',[1 1 1])
    %legend(label_long{i},'Interpreter','latex','Location','best')
    iFig = i+k; 
end

k = iFig;
for i = 1:length(FT_lat)
    figure(iFig);
    bode(FT_lat{i},'b-')
    grid on; hold all;
    margin(FT_lat{i},{10^-5,10^2})
    set(gcf,'Color',[1 1 1])
    %legend(label_lat{i},'Interpreter','latex','Location','best')
    iFig = i+k;
end


%% Diagramas de Nichols
k = iFig;
for i = 1:length(FT_long)
    [modulo_bode fase_bode] = bode(FT_long{i},{10^-3,10^4});
    modulodB_bode = squeeze(20*log10(modulo_bode));
    faseDeg_bode = squeeze(fase_bode);

    figure(iFig)
    plot(faseDeg_bode,modulodB_bode,'Color','b','Linewidth',2)
    grid on; hold on; 
    plot(faseDeg_bode(1),modulodB_bode(1),'ro','Linewidth',2)
    hold on; 
    plot([180 180],[0 100],'k-','Linewidth',2)
    plot([-180 -180],[0 100],'k-','Linewidth',2)
    xlabel('Open-Loop Phase [deg])'); 
    ylabel('Open-Loop Gain [dB]');
    set(gca,'XLim',[min(faseDeg_bode) max(faseDeg_bode)]);
    % Ticks de separacion
    hold on; set(gca,'XTick',[-720:45:720]); % Grados
    hold on; set(gca,'YTick',[-150:10:100]); % dB
    set(gcf,'Color',[1 1 1])
    %legend(label_long{i},'Interpreter','latex','Location','best')
    iFig = i+k;

end

k = iFig;
for i = 1:length(FT_lat)
    [modulo_bode fase_bode] = bode(FT_lat{i},{10^-5,10^2});
    modulodB_bode = squeeze(20*log10(modulo_bode));
    faseDeg_bode = squeeze(fase_bode);

    figure(iFig)
    plot(faseDeg_bode,modulodB_bode,'Color','b','Linewidth',2)
    grid on; hold on; 
    plot(faseDeg_bode(1),modulodB_bode(1),'ro','Linewidth',2)
    hold on; 
    plot([180 180],[0 100],'k-','Linewidth',2)
    plot([-180 -180],[0 100],'k-','Linewidth',2)
    xlabel('Open-Loop Phase [deg]'); 
    ylabel('Open-Loopn Gain [dB]');
    set(gca,'XLim',[min(faseDeg_bode) max(faseDeg_bode)]);
    % Ticks de separacion
    hold on; set(gca,'XTick',[-720:45:720]); % Grados
    hold on; set(gca,'YTick',[-150:10:100]); % dB
    set(gcf,'Color',[1 1 1])
    %legend(label_lat{i},'Interpreter','latex','Location','best')
    iFig = i+k;
end

%% Respuesta escalón 
    % Análisis de la respuesta estacionaria
k = iFig;
for i = 1:length(FT_long)
    figure(iFig);
    step(FT_long{i})
    grid on;  

    stepInfo_Long{i} = stepinfo(FT_long{i});
    %legend(label_long{i},'Interpreter','latex','Location','best')
    iFig = i+k;
end
k = iFig;
for i = 1:length(FT_lat)
    figure(iFig);
    step(FT_lat{i})
    grid on;

    stepInfo_Lat{i} = stepinfo(FT_lat{i});
    %legend(label_lat{i},'Interpreter','latex','Location','best')
    iFig = i+k;
end
    % Análsis de la respueta transitoria
k = iFig;
for i = 1:length(FT_long)
    figure(iFig);
    step(FT_long{i},20)
    grid on;  
    %legend(label_long{i},'Interpreter','latex','Location','best')
    iFig = i+k;
end
k = iFig;
for i = 1:length(FT_lat)
    figure(iFig);
    step(FT_lat{i},20)
    grid on;
    %legend(label_lat{i},'Interpreter','latex','Location','best')
    iFig = i+k;
end


%% Respuesta a rampa unitaria
    % Análisis de la respuesta estacionaria
t = (0:1:400)';
y = ones(length(t),1); 
for i = 1:11
    y(i) = (i-1)/10;
end

k = iFig;
for i = 1:length(FT_long)
    figure(iFig);
    lsim(FT_long{i},y,t)
    grid on;

    [Y,t_] = lsim(FT_long{i},y,t);
    rampaInfo_Long{i} = lsiminfo(Y,t_);
    %legend(label_long{i},'Interpreter','latex','Location','best')
    iFig = i+k;
end
k = iFig;
for i = 1:length(FT_lat)
    figure(iFig);
    lsim(FT_lat{i},y,t)
    grid on; 

    [Y,t_] = lsim(FT_lat{i},y,t);
    rampaInfo_Lat{i} = lsiminfo(Y,t_);
    %legend(label_lat{i},'Interpreter','latex','Location','best')
    iFig = i+k;
end
    % Análisis de la respuesta transitoria
t = (0:1:400)';
y = ones(length(t),1); 
for i = 1:11
    y(i) = (i-1)/10;
end


k = iFig;
for i = 1:length(FT_long)
    figure(iFig);
    lsim(FT_long{i},y,t)
    grid on; 
    axis([0 40 0 10])
    %legend(label_long{i},'Interpreter','latex','Location','best')
    iFig = i+k;
end
k = iFig;
for i = 1:length(FT_lat)
    figure(iFig);
    lsim(FT_lat{i},y,t)
    grid on; 
    axis([0 40 0 10])
    %legend(label_lat{i},'Interpreter','latex','Location','best')
    iFig = i+k;
end
