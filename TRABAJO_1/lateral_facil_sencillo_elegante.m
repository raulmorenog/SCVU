%% Obtención de las FT para el Beech 99
%Una vez finalizada la ejecución todos los datos se almacenan en la
%estructura FT_lat
clc
clear all 
%close all 

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

% Derivadas de estabilidad lateral-direccionales 
p.Cl_beta = -0.130; p.Cl_p = -0.50; p.Cl_r = 0.14; 
p.Cy_beta = -0.59; p.Cy_p = -0.19; p.Cy_r = 0.39;
p.Cn_beta = 0.08; p.Cn_p = 0.019; p.Cn_r = -0.197;  %Cn_p puesto positivo aunque en los datos era negativo.
p.Cn_Tbeta = 0; 
p.Cl_deltaA = 0.156; p.Cl_deltaR = -0.0106; 
p.Cy_deltaA = 0; p.Cy_deltaR = -0.144; 
p.Cn_deltaA = -0.0012 ; p.Cn_deltaR = 0.0758;

%% Derivadas de estabilidad y adimensionalización --> los guardamos en la estructura s
s.mu = p.m/(p.rhos*p.Sw*p.b/2);
s.I_xxnd = p.I_xx/(p.rhos*p.Sw*(p.b/2)^3);
s.I_zznd = p.I_zz/(p.rhos*p.Sw*(p.b/2)^3);
s.I_xznd = p.I_xz/(p.rhos*p.Sw*(p.b/2)^3);

s.Cxs = p.Cxs;
s.Czs = p.Czs;

s.Cy_beta = p.Cy_beta;    
s.Cy_p = p.Cy_p;    
s.Cy_r = p.Cy_r;    
s.Cy_deltaA = p.Cy_deltaA;    
s.Cy_deltaR = p.Cy_deltaR;   

s.Cl_beta = p.Cl_beta;    
s.Cl_p = p.Cl_p;    
s.Cl_r = p.Cl_r;    
s.Cl_deltaA = p.Cl_deltaA;
s.Cl_deltaR = p.Cl_deltaR;

s.Cn_beta = p.Cn_beta;
s.Cn_p = p.Cn_p;    
s.Cn_r = p.Cn_r;    
s.Cn_deltaA = p.Cn_deltaA;
s.Cn_deltaR = p.Cn_deltaR;

%% Definición del sistema --> Se definen las matrices de la ecuación de estado x_t = F*x + B*u
sys.F(1,1) = s.Cy_beta / (2*s.mu);
sys.F(1,2) = s.Cy_p / (2*s.mu);
sys.F(1,3) = -(2*s.mu - s.Cy_r)/(2*s.mu);
sys.F(1,4) = - s.Czs / (2*s.mu);
sys.F(2,1) = (s.Cl_beta*s.I_zznd + s.Cn_beta*s.I_xznd) / (s.I_xxnd*s.I_zznd - s.I_xznd^2);
sys.F(2,2) = (s.Cl_p*s.I_zznd + s.Cn_p*s.I_xznd) / (s.I_xxnd*s.I_zznd - s.I_xznd^2);
sys.F(2,3) = (s.Cl_r*s.I_zznd + s.Cn_r*s.I_xznd) / (s.I_xxnd*s.I_zznd - s.I_xznd^2);
sys.F(2,4) = 0;
sys.F(3,1) = (s.Cn_beta*s.I_xxnd + s.Cl_beta*s.I_xznd) / (s.I_xxnd*s.I_zznd - s.I_xznd^2);
sys.F(3,2) = (s.Cn_p*s.I_xxnd + s.Cl_p*s.I_xznd) / (s.I_xxnd*s.I_zznd - s.I_xznd^2);
sys.F(3,3) = (s.Cn_r*s.I_xxnd + s.Cl_r*s.I_xznd) / (s.I_xxnd*s.I_zznd - s.I_xznd^2);
sys.F(3,4) = 0;
sys.F(4,1) = 0;
sys.F(4,2) = 1;
sys.F(4,3) = 0;
sys.F(4,4) = 0;

sys.B(1,1) = 0;
sys.B(1,2) = s.Cy_deltaR / (2*s.mu);
sys.B(2,1) = (s.Cl_deltaA*s.I_zznd + s.Cn_deltaA*s.I_xznd) / (s.I_xxnd*s.I_zznd - s.I_xznd^2);
sys.B(2,2) = (s.Cl_deltaR*s.I_zznd + s.Cn_deltaR*s.I_xznd) / (s.I_xxnd*s.I_zznd - s.I_xznd^2);
sys.B(3,1) = (s.Cn_deltaA*s.I_xxnd + s.Cl_deltaA*s.I_xznd) / (s.I_xxnd*s.I_zznd - s.I_xznd^2);
sys.B(3,2) = (s.Cn_deltaR*s.I_xxnd + s.Cl_deltaR*s.I_xznd) / (s.I_xxnd*s.I_zznd - s.I_xznd^2);
sys.B(4,1) = 0;
sys.B(4,2) = 0;

%% Saco las funciones de transferencia adimensionalizadas
sys.sistema = ss(sys.F,sys.B,eye(4),zeros(4,2));
FT_nd = tf(sys.sistema);
% Guardo las tf en un cell para operar mejor
for i = 1:4
    for j = 1:2
        FT_lat_nd{i,j} = tf([FT_nd.Numerator{i,j}],[FT_nd.Denominator{i,j}]); 
    end
end

%% Dimensionalizamos
% El denominador
for i = 1:length(FT_nd.Denominator{1})
    j = i - 1;
   FT.Denominator(length(FT_nd.Denominator{1}) - j) = ...
       FT_nd.Denominator{1}(length(FT_nd.Denominator{1}) - j)*(p.b/(2*p.Us))^j; 
end
% Los numeradores 
for i = 1:4
    for m = 1:2
        for j = 1:length(FT_nd.Numerator{i,m})
            k = j - 1;
            FT.Numerator{i,m}(length(FT_nd.Numerator{i,m}) - k) = ...
                FT_nd.Numerator{i,m}(length(FT_nd.Numerator{i,m}) - k)*(p.b/(2*p.Us))^k; 
        end
        FT_lat.nofact{i,m} = tf(FT.Numerator{i,m},FT.Denominator);
    end
end
%Damos dimensiones de velocidad entre grados
FT_lat.nofact{2,1} = FT_lat.nofact{2,1}*(2*p.Us/p.b);
FT_lat.nofact{3,1} = FT_lat.nofact{3,1}*(2*p.Us/p.b);
FT_lat.nofact{2,2} = FT_lat.nofact{2,2}*(2*p.Us/p.b);
FT_lat.nofact{3,2} = FT_lat.nofact{3,2}*(2*p.Us/p.b);

%% Factorizamos las funciones de transferencia
FT_lat.fact = struct('deltaA_beta',zpk(FT_lat.nofact{1,1}),...
                     'deltaA_p',zpk(FT_lat.nofact{2,1}),...
                     'deltaA_r',zpk(FT_lat.nofact{3,1}),...
                     'deltaA_phi',zpk(FT_lat.nofact{4,1}),...
                     'deltaR_beta',zpk(FT_lat.nofact{1,2}),...
                     'deltaR_p',zpk(FT_lat.nofact{2,2}),...
                     'deltaR_r',zpk(FT_lat.nofact{3,2}),...
                     'deltaR_phi',zpk(FT_lat.nofact{4,2}));
                 
%% Cálculo de propiedades de los modos
[wn, amort, Poles] = damp(FT_lat.fact.deltaA_beta);
FT_lat.espiral.wn = wn(1);
FT_lat.espiral.poles = Poles(1);
FT_lat.espiral.tau_S = 1/wn(1);
FT_lat.espiral.t12 = -log(2)/real(FT_lat.espiral.poles);
FT_lat.balance.wn = wn(4);
FT_lat.balance.poles = Poles(4);
FT_lat.balance.tau_R = 1/wn(4);
FT_lat.balance.t12 = -log(2)/real(FT_lat.balance.poles);
FT_lat.dutchroll.poles = Poles([2:3]);
FT_lat.dutchroll.wn = wn(2);
FT_lat.dutchroll.period = 2*pi/FT_lat.dutchroll.wn;
FT_lat.dutchroll.amort = amort(2);
FT_lat.dutchroll.t12 = -log(2)/real(FT_lat.dutchroll.poles(1)); 

