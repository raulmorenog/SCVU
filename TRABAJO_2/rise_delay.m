% Función para el cálculo del rise time y el time delay de una respuesta
function [rise_time time_delay] = rise_delay(t,signal)
    time = t; 
    y = signal;

    % Rise time (tiempo entre 0.9y_ss y 0.1y_ss)
        % Cálculo del valor máximo (estacionario para phi) de phi y su
        % tiempo asociado
    [y_max index_max] = max(y);
    t_max = time(index_max);
        
        % Cálculo de los valores de referencia y del tiempo asociado a esto
    y_01 = 0.1*y_max; 
    y_09 = 0.9*y_max;
    for i = 1:length(time)
        if abs(y_01-y(i)) <= 0.5     % Impuesta una diferencia de 0.1 entre valores
            index_01 = i; 
        elseif abs(y_09-y(i)) <= 0.1 % Impuesta una diferencia de 0.1 entre valores
            index_09 = i; 
            break           % Cogemos el primer valor y salimos del bucle porque
        else                % la respuesta puede disminuir (caso espiral estable)
        end
    end
    t_01 = time(index_01); 
    t_09 = time(index_09);
    rise_time = t_09-t_01;

    % Time delay (tiempo de reacción)
        % Cálculo de la derivada local máxima en la subida. La recta de
        % esta pendiente cortará al eje de tiempos dando el time delay
    p = polyfit(time,y,20);         % Se ajustan los datos con un polinomio
    time_pol = 0:0.05:t_max;        % de orden 20 (muy alto), y así se toma un
    y_pol = polyval(p,time_pol);    % paso de tiempo uniforme con el que trabajar
%     figure()
%     plot(time,y,'r-',time_pol,y_pol,'b--')

    paso = (t_max-0)/(length(time_pol));
    for i = 2:(length(time_pol)-1)                  % Calculamos la derivada de la señal
        F(i) = (y_pol(i+1)-y_pol(i-1))/(2*paso);   % solo en el tramo de subida (diferencias   
    end                                            % finitas centradas) 
    [F_max index_Fmax] = max(F);

        % Ecuación de la recta para y=0 (0-y0=Fx-Fx0)
    time_delay = time_pol(index_Fmax)-y_pol(index_Fmax)/F_max; 

end