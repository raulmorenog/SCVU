%% Función cálculo FT del movimiento longitudinal de avión

function [m] = longitudinal(p)
    % Se calculan las funciones de transferencia longitudinales a partir de
    % las derivadas de estabilidad del avión. Se proporcionan como dato en
    % estructuras y se dan como salida las funciones de transferenciano sus
    % numeradores y denominadores correspondientes. Además de ello, algunos
    % datos importantes como las frecuencias, amortiguamientos y tiempos de
    % los modos fugoide y de corto periodo
    
        % Parámetros adimensionales neccesarios 
    mu = p.m/(p.Sw*p.rhos*p.b/2);           % Parámetro másico
    Ix = p.I_xx/(p.Sw*p.rhos*(p.b/2)^3);
    Iy = p.I_yy/(p.Sw*p.rhos*(p.b/2)^3);
    Iz = p.I_zz/(p.Sw*p.rhos*(p.b/2)^3);

        % Se determinan los coeficientes de fuerzas 
    Cx_u = p.CT_xu*cos(p.epsilons)-p.Cd_u; 
    Cz_u = p.CT_xu*sin(p.epsilons)-(p.Ms^2/(1-p.Ms^2))*p.Cls; 
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

        % Se escriben los coeficientes de la cuártica de estabilidad
    A = 2*mu*Iy*(2*mu-Cz_alphap);
    B = -2*mu*Iy*(Cz_alpha+Cx_u)+Iy*Cx_u*Cz_alphap...
        -2*mu*(Cz_q*Cm_alphap-Cm_q*Cz_alphap)...
        -4*mu^2*(Cm_alphap+Cm_q);
    C = Iy*(Cx_u*Cz_alpha-Cx_alpha*Cz_u)...
        +2*mu*(Cz_alpha*Cm_q-Cm_alpha*Cz_q+Cx_u*Cm_q+Cx_u*Cm_alphap)...
        -4*mu^2*Cm_alpha-Cx_u*(Cm_q*Cz_alphap-Cz_q*Cm_alphap)-2*Iy*Czs*Cx_alpha;
    D = -2*Czs^2*Cm_alphap+2*mu*(Cx_u*Cm_alpha-Cx_alpha*Cm_u-Czs*Cm_u)...
        +Cx_u*(Cm_alpha*Cz_q-Cm_q*Cz_alpha)-Cx_alpha*(Cm_u*Cz_q-Cm_q*Cz_u)...
        +Czs*(Cm_u*Cz_alphap-Cz_u*Cm_alphap)+2*Czs*Cm_q*Cx_alpha;
    E = Czs*(Cm_u*Cz_alpha-Cm_alpha*(2*Czs+Cz_u));
    
    den = [A B C D E]
    m.coef_cuartica = den;

        % Se analiza la cuártica --> Frecuencias y amortiguamientos
    a = roots([A B C D E]);
    b = fsolve(@(x) modos(x,a(1),a(3)),[2.8 0.35 0.09 0.11]);
    w_fug = abs(b(3)); chi_fug = abs(b(4)); 
    w_cp = abs(b(1)); chi_cp = abs(b(2));
    m.fugoide.w_fug = w_fug; 
    m.fugoide.chi_fug = chi_fug;
    m.fugoide.t_12 = -0.693/real(a(3));
    m.cortoperiodo.w_cp = w_cp; 
    m.cortoperiodo.chi_cp = chi_cp; 
    m.cortoperiodo.t_12 = -0.693/real(a(1));

        % Numeradores de las funciones de transferencia 
    Bu = Cx_deltae*Iy*(2*mu-Cz_alphap); 
    Cu = Cx_deltae*(Cm_q*(Cz_alphap-2*mu)-Cz_alpha*Iy-Cm_alphap*(2*mu+Cz_q))...
        +Cz_deltae*Cx_alpha*Iy...
        +Cm_deltae*(Cx_alpha*(2*mu+Cz_q)+Czs*(2*mu-Cz_alphap)); 
    Du = Cx_deltae*(Cz_alpha*Cm_q-(2*mu+Cz_q)*Cm_alpha)...
        +Cz_deltae*(Czs*Cm_alphap-Cx_alpha*Cm_q)...
        -Cm_deltaep*Czs+Cm_deltae*(Cx_alpha*(2*mu+Cz_q)+(2*mu-Cz_alpha)*Czs); 
    Eu = Czs*(Cz_deltae*Cm_alpha-Cm_deltae*Cz_alpha); 
    Nu_deltae = [Bu Cu Du Eu];  

    Balpha = 2*mu*(Cz_deltae*Iy+Cm_deltaep*(2*mu+Cz_q)); 
    Calpha = Cx_deltae*Iy*(Cz_u+2*Czs)-Cz_deltae*(2*mu*Cm_q+Iy*Cx_u)...
        -Cm_deltaep*Cx_u*(2*mu+Cz_q)+Cm_deltae*2*mu*(2*mu+Cz_q); 
    Dalpha = -Cx_deltae*(Cm_q*(Cm_u+2*Czs)-Cm_u*(2*mu+Cz_q))...
        +Cz_deltae*Cx_u*Cm_q+Cm_deltaep*Czs*(Cz_u+2*Czs)...
        -Cm_deltae*Cx_u*(2*mu+Cz_q);
    Ealpha = Czs*(-Cm_u*Cz_deltae+Cm_deltae*(Cz_u+2*Czs)); 
    Nalpha_deltae = [Balpha Calpha Dalpha Ealpha];
    
    Btheta = Cm_deltae*2*mu*(2*mu-Cz_alphap);
    Ctheta = Cz_deltae*2*mu*Cm_alphap...
        -Cm_deltaep*(2*mu*Cz_alpha+Cx_u*(2*mu-Cz_alphap))...
        +Cm_deltae*2*mu*(2*mu-Cz_alphap); 
    Dtheta = Cx_deltae*(Cm_alphap*(Cz_u+2*Czs)+Cx_u*(2*mu-Cz_alphap))...
        -Cm_deltae*(2*mu*Cz_alpha+Cx_u*(2*mu-Cz_alphap))...
        -Cz_deltae*(Cx_u*Cm_alphap-2*mu*Cm_alpha)...
        +Cm_deltaep*(Cx_u*Cz_alpha-Cx_alpha*(Cz_u+2*Czs));
    Etheta = Cx_deltae*(Cm_alpha*(Cz_u+2*Czs)-Cm_u*Cz_alpha)...
        -Cz_deltae*(Cx_u*Cm_alpha-Cm_u*Cx_alpha)...
        +Cm_deltae*(Cx_u*Cz_alpha-Cx_alpha*(Cz_u+2*Czs)); 
    Ntheta_deltae = [Btheta Ctheta Dtheta Etheta]; 

    % Funciones de tranferencia del movimiento longitudinal
    m.Gu_deltae = tf(Nu_deltae, den); 
    m.Galpha_deltae = tf(Nalpha_deltae, den);
    m.Gtheta_deltae = tf(Ntheta_deltae, den);




end

function F = modos(x,a,b)       % revisar porque algo no va del todo bien
   F(1) = -x(1)*x(2)-real(a); 
   F(2) = x(1)*sqrt(x(2)^2-1)-imag(a); 
   F(3) = -x(3)*x(4)-real(b); 
   F(4) = x(3)*sqrt(x(4)^2-1)-imag(b); 
end