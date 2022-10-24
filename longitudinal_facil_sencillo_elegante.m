%% Obtención de las FT para el Beech 99
%Una vez finalizada la ejecución todos los datos se almacenan en la
%estructura FT_long
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


%% Derivadas de estabilidad y adimensionalización --> los guardamos en la estructura s
s.mu_long = p.m/(p.rhos*p.Sw*p.c/2);
s.I_yynd = p.I_yy/(p.rhos*p.Sw*(p.c/2)^3);

s.Cxs = p.Cxs;
s.Cx_u = p.CT_xu*cos(p.epsilons)-p.Cd_u;     
s.Cx_alpha = p.Cls-p.Cd_alpha;     
s.Cx_alphap = 0;                                  
s.Cx_deltae = -p.Cd_deltae;                                   

s.Czs = p.Czs;
s.Cz_u = -p.CT_xu*sin(p.epsilons)-p.Cl_u;    
s.Cz_alpha = -p.Cl_alpha-p.Cds;    
s.Cz_alphap = -p.Cl_alphap;    
s.Cz_q = -p.Cl_q;    
s.Cz_deltae = -p.Cl_deltae; 

s.Cm_u = p.Cm_u;
s.Cm_alpha = p.Cm_alpha;           
s.Cm_alphap = p.Cm_alphap;
s.Cm_q = p.Cm_q;     
s.Cm_deltae = p.Cm_deltae;
s.Cm_deltaep = p.Cm_deltaep;     
    
%% Definición del sistema --> Se definen las matrices de la ecuación de estado x_t = F*x + B*u
sys.F(1,1) = s.Cx_u/2/s.mu_long;
sys.F(1,2) = s.Cx_alpha / (2*s.mu_long);
sys.F(1,3) = 0;
sys.F(1,4) = s.Czs/2/s.mu_long;
sys.F(2,1) = (s.Cz_u + 2*s.Czs) /(2*s.mu_long - s.Cz_alphap);
sys.F(2,2) = s.Cz_alpha / (2*s.mu_long - s.Cz_alphap);
sys.F(2,3) = (2*s.mu_long + s.Cz_q) /(2*s.mu_long - s.Cz_alphap);
sys.F(2,4) = 0;
sys.F(3,1) = s.Cm_u / s.I_yynd + s.Cm_alphap / s.I_yynd * (s.Cz_u + 2*s.Czs) / (2*s.mu_long - s.Cz_alphap); 
sys.F(3,2) = s.Cm_alpha / s.I_yynd + s.Cm_alphap / s.I_yynd * s.Cz_alpha / (2*s.mu_long - s.Cz_alphap);
sys.F(3,3) = s.Cm_q / s.I_yynd + s.Cm_alphap / s.I_yynd * (2*s.mu_long + s.Cz_q) /(2*s.mu_long - s.Cz_alphap);
sys.F(3,4) = 0;
sys.F(4,1) = 0; 
sys.F(4,2) = 0;
sys.F(4,3) = 1;
sys.F(4,4) = 0;

sys.B(1,1) = s.Cx_deltae / 2 / s.mu_long;
sys.B(2,1) = s.Cz_deltae / (2*s.mu_long - s.Cz_alphap);
sys.B(3,1) = s.Cm_deltae / s.I_yynd + s.Cm_alphap / ...
    s.I_yynd * s.Cz_deltae / (2*s.mu_long - s.Cz_alphap);
sys.B(4,1) = 0;
    
%% Saco las funciones de transferencia adimensionalizadas
sys.sistema = ss(sys.F,sys.B,eye(4),zeros(4,1));
FT_nd = tf(sys.sistema);
% Guardo las tf en un cell para operar mejor
for i = 1:4
   FT_long_nd{i} = tf([FT_nd.Numerator{i}],[FT_nd.Denominator{i}]); 
end

%% Dimensionalizamos
% El denominador
for i = 1:length(FT_nd.Denominator{1})
    j = i - 1;
   FT.Denominator(length(FT_nd.Denominator{1}) - j) = ...
       FT_nd.Denominator{1}(length(FT_nd.Denominator{1}) - j)*(p.c/(2*p.Us))^j; 
end
% Los numeradores 
for i = 1:4
    for j = 1:length(FT_nd.Numerator{i})
        k = j - 1;
        FT.Numerator{i}(length(FT_nd.Numerator{i}) - k) = ...
            FT_nd.Numerator{i}(length(FT_nd.Numerator{i}) - k)*(p.c/(2*p.Us))^k; 
    end
    FT_long.nofact{i} = tf(FT.Numerator{i},FT.Denominator);
end
FT_long.nofact{1} = FT_long.nofact{1}*p.Us*pi/180;  %Damos dimensiones de velocidad entre grados

%% Factorizamos las funciones de transferencia
FT_long.fact = struct('deltae_u',zpk(FT_long.nofact{1}),...
                     'deltae_alpha',zpk(FT_long.nofact{2}),...
                     'deltae_q',zpk(FT_long.nofact{3}),...
                     'deltae_theta',zpk(FT_long.nofact{4}));

%% Cálculo de propiedades de los modos
[wn, amort, Poles] = damp(FT_long.fact.deltae_u);
FT_long.phugoid.poles = Poles([1:2]);
FT_long.phugoid.wn = min(wn);
FT_long.phugoid.amort = min(amort);
FT_long.phugoid.period = 2*pi/FT_long.phugoid.wn;
FT_long.phugoid.t12 = -log(2)/real(FT_long.phugoid.poles(1));
FT_long.shortperiod.poles = Poles([3:4]);
FT_long.shortperiod.wn = max(wn);
FT_long.shortperiod.amort = max(amort);
FT_long.shortperiod.period = 2*pi/FT_long.shortperiod.wn;
FT_long.shortperiod.t12 = -log(2)/real(FT_long.shortperiod.poles(1));

%% FT de deltae a u 
display(FT_long.fact.deltae_u)
% Diagrama de Bode
Bode1 = figure(1);
bode(FT_long.fact.deltae_u)
grid on; set(gcf,'Color',[1 1 1]); 
% Diagrama de Nichols
Nichols1 = figure(2);
nichols(FT_long.fact.deltae_u)
grid on; set(gcf,'Color',[1 1 1]);
% Respuesta escalón
Step1 = figure(3);
step(FT_long.fact.deltae_u); set(gcf,'Color',[1 1 1]);
% Respuesta rampa unitarla
t = 0:1:900;
Rampa1 = figure(4);
lsim(FT_long.fact.deltae_u,t,t);

%% FT de deltae a alpha
display(FT_long.fact.deltae_alpha)
% Diagrama de Bode
Bode2 = figure(5);
bode(FT_long.fact.deltae_alpha)
grid on; set(gcf,'Color',[1 1 1]); 
% Diagrama de Nichols
Nichols2 = figure(6);
nichols(FT_long.fact.deltae_alpha)
grid on; set(gcf,'Color',[1 1 1]);
% Respuesta escalón
Step2 = figure(7);
step(FT_long.fact.deltae_alpha); set(gcf,'Color',[1 1 1]);
% Respuesta rampa unitarla
Rampa2 = figure(8);
lsim(FT_long.fact.deltae_alpha,t,t);

%% FT de deltae a theta 
display(FT_long.fact.deltae_theta)
% Diagrama de Bode
Bode3 = figure(9);
bode(FT_long.fact.deltae_theta)
grid on; set(gcf,'Color',[1 1 1]); 
% Diagrama de Nichols
Nichols3 = figure(10);
nichols(FT_long.fact.deltae_theta)
grid on; set(gcf,'Color',[1 1 1]);
% Respuesta escalón
Step3 = figure(11);
step(FT_long.fact.deltae_theta); set(gcf,'Color',[1 1 1]);
% Respuesta rampa unitarla
Rampa3 = figure(12);
lsim(FT_long.fact.deltae_theta,t,t);

%% Ft de deltae a q (no haría falta)
display(FT_long.fact.deltae_q)
% Diagrama de Bode
Bode4 = figure(13);
bode(FT_long.fact.deltae_q)
grid on; set(gcf,'Color',[1 1 1]); 
% Diagrama de Nichols
Nichols4 = figure(14);
nichols(FT_long.fact.deltae_q)
grid on; set(gcf,'Color',[1 1 1]);
% Respuesta escalón
Step4 = figure(15);
step(FT_long.fact.deltae_q); set(gcf,'Color',[1 1 1]);
% Respuesta rampa unitarla
Rampa4 = figure(16);
lsim(FT_long.fact.deltae_q,t,t);