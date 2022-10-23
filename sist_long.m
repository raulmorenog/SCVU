%% Funciones de Transferencia Longitudinales
function  long = sist_long(p)
        % Adimensionalización longitudinal y coeficientes
    mu_long = p.m/(p.rhos*p.Sw*p.c/2);
    I_yynd = p.I_yy/(p.rhos*p.Sw*(p.c/2)^3);
    
    Cxs = p.Cxs;
    Cx_u = p.CT_xu*cos(p.epsilons) - p.Cd_u;     
    Cx_alpha = p.Cls - p.Cd_alpha;     
    Cx_alphap = 0;                                  
    Cx_deltae = -p.Cd_deltae;                                   
    
    Czs = p.Czs;
    Cz_u = -p.CT_xu*sin(p.epsilons) - p.Cl_u;    
    Cz_alpha = -p.Cl_alpha - p.Cds;    
    Cz_alphap = -p.Cl_alphap;    
    Cz_q = -p.Cl_q;    
    Cz_deltae = -p.Cl_deltae; 
    
    Cm_u = p.Cm_u;
    Cm_alpha = p.Cm_alpha;           
    Cm_alphap = p.Cm_alphap;
    Cm_q = p.Cm_q;     
    Cm_deltae = p.Cm_deltae;
    Cm_deltaep = p.Cm_deltaep;        
 
        % Construcción del Sistema Longitudinal: [A]*{Y}=[b]*{X} 
    syms s z
    sympref('FloatingPointOutput',true)
    
    A = [2*mu_long*s - Cx_u, -Cx_alpha, -Czs;...
        -(Cz_u + 2*Czs), (2*mu_long - Cz_alphap)*s - Cz_alpha, -(2*mu_long + Cz_q)*s;...
        -Cm_u, -(Cm_alphap*s + Cm_alpha), I_yynd*s^2 - Cm_q*s];
    b = [Cx_deltae; Cz_deltae; Cm_deltaep*s + Cm_deltae];
    
    % Cramer y factorización con s ADIMENSIONAL (s gorro) --> FT con s
    % adimensional
    Den = sym2poly(det(A));     % Cuártica de estabilidad Longitudinal
    
    Num_deltae_u = sym2poly(det([b,A(:,2),A(:,3)]));
    TFdeltae_u_nd = tf([Num_deltae_u],[Den]);   % Adimensional!!                                                   
    TFdeltae_u = TFdeltae_u_nd*tf([p.Us]);      % Dimensional en u, todavia no en s
    
    Num_deltae_alpha = sym2poly(det([A(:,1),b,A(:,3)]));
    TFdeltae_alpha = tf([Num_deltae_alpha],[Den]); 
    
    Num_deltae_theta=sym2poly(det([A(:,1),A(:,2),b]));
    TFdeltae_theta=tf([Num_deltae_theta],[Den]); 
    
    %Dimensionalización de s (s = z*c/2Us) --> FT con s dimensional
    DenDim = sym2poly(subs(det(A),s,z*p.c/(2*p.Us)));
    
    Num_deltae_uDim = sym2poly(subs(det([b,A(:,2),A(:,3)]),s,z*p.c/(2*p.Us)));
    TFdeltae_uDim = tf([Num_deltae_uDim],[DenDim])*tf([p.Us]);
    
    Num_deltae_alphaDim = sym2poly(subs(det([A(:,1),b,A(:,3)]),s,z*p.c/(2*p.Us)));
    TFdeltae_alphaDim = tf([Num_deltae_alphaDim],[DenDim]);
    
    Num_deltae_thetaDim = sym2poly(subs(det([A(:,1),A(:,2),b]),s,z*p.c/(2*p.Us)));
    TFdeltae_thetaDim = tf([Num_deltae_thetaDim],[DenDim]);

    % Cálculo de las propiedades de las FT 
        % Zeros, polos y ganancias
    [zeros_u,polos_u,ganancia_u] = tf2zp(Num_deltae_uDim*p.Us,DenDim);
    %F = factorizacion(Num_deltae_uDim*p.Us,DenDim)
    [zeros_alpha,polos_alpha,ganancia_alpha] = tf2zp(Num_deltae_alphaDim,DenDim);
    [zeros_theta,polos_theta,ganancia_theta] = tf2zp(Num_deltae_thetaDim,DenDim);
        
        % Amortiguamientos y frecuencias de los modos
    [w_natural,amort,polos] = damp(TFdeltae_uDim); 
    
    omega_fug = w_natural(1);  
    chi_fug = amort(1);
    t12_fug = -0.693/real(polos(1));
    omega_cp = w_natural(3); 
    chi_cp = amort(3); 
    t12_cp = -0.693/real(polos(3));
    
    % Salida de la función
        % Funciones de transferencia no factorizadas dimensionales
    long.FT_nfact = struct( ...
    'FTdeltaE_u',TFdeltae_uDim, ...
    'FTdeltaE_alpha',TFdeltae_alphaDim, ...
    'FTdeltaE_theta',TFdeltae_thetaDim);
        
        % Funciones de transferencia factorizadas dimensionales
    long.FT_fact = struct( ...
    'FTdeltaE_u',zpk(zeros_u,polos_u,ganancia_u), ...
    'FTdeltaE_alpha',zpk(zeros_alpha,polos_alpha,ganancia_alpha), ...
    'FTdeltaE_theta',zpk(zeros_theta,polos_theta,ganancia_theta));

        % Propiedades de los modos 
    long.Prop_mod = struct(...
        'Omega_Cp',omega_cp,...
        'Amort_Cp',chi_cp,...
        'T12_Cp',t12_cp, ...
        'Omega_Fug',omega_fug,...
        'Amort_Fug',chi_fug,...
        'T12_Fug',t12_fug);

end

% function TF = factorizacion(num,den,a)
%     % Función para la factorización de las FT de la forma que se da en los 
%     % apuntes. Devuelve la función de transferencia, con la entrada del
%     % numerador,num, el denominador, den, y un vector de [zeros, polos,
%     % ganancia], a .
% 
%     [zeros,polos,ganancia] = tf2zp(num,den);
%     [w_n,amort,polos] = damp(tf(num,den));
% 
%     w_fug = w_n(1);  
%     chi_fug = amort(1);
%     w_cp = w_n(3); 
%     chi_cp = amort(3); 
%     
%     for i = 1:length(zeros)
%         if isreal(zeros(i)) == 1 
%             tau(i) = 1/zeros(i)
%             w(i) = 0; 
%             chi(i) = 0;
%         else 
%             tau(i) = 0;
%             [w(i),chi(i)] = fsolve(@(x) w_amor(x,zeros(i)),[w_fug,chi_fug]);
%         end
%     end
%     k = ganancia;
%     for i = 1:length(zeros)
%         if w(i) == 0 & chi(i) == 0 
%             k = k/zeros(i);
%         else 
%             k = k/w(i)^2; 
%         end
%     end 
%     K = k/(w_fug*w_cp)^2
%     
%     FT = tf(K*[tau(1) 1]*[tau(2) 1],[]*[])
%     i=1
% end
% 
% function F = w_amor(x,raiz)
%     % x(1) = frecuencia natural ; x(2) = amortiguamiento;
%     F(1) = x(1)*x(2)-real(raiz);
%     F(2) = x(1)*sqrt(x(2)^2-1)-imag(raiz); 
% end