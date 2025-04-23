clearvars
close all
clc

% opciones del usuario
disp("A: Ingresar señales.");
disp("B: Señales aleatorias.");
opcion = input("Elige una opción: ", "s");

% funcion para generar señales aleatorias, magia
generarX = @(largo) -1 + 2*rand(1, largo);

% qof, se puede elegir entre inrgesar señales o generar aleatoriamente
% la generacion de señales es de largo 8 y un filtro de largo 6 cuyos
% valores varian entre -1 y 1
if opcion == 'A'
    x_n = str2num(input("Ingrese la primera señal, separada por comas: ", "s"));
    h_n = str2num(input("Ingrese el filtro, separada por comas: ", "s"));

elseif opcion == 'B'
    disp("Generando señales...");
    % señal largo 8
    x_n = generarX(8);
    % filtro largo 6
    h_n = generarX(6);
    disp("x[n]:"); disp(x_n);
    disp("h[n]:"); disp(h_n);

else 
    error("Intentalo de nuevo...");
end

% convolucion normal
disp("Convolución de Matlab:");
y_matlab = conv(x_n, h_n);
disp('Resultado conv():');
disp(y_matlab);

% convolucion con matriz Toeplitz
disp("Convolución con matriz Toeplitz:");
y_toeplitz = toep_convolucion(x_n, h_n);
disp("Resultado Toeplitz:");
disp(y_toeplitz);

% convolucion usando Winograd
disp("Convolución Winograd:");
y_winograd = winograd(x_n, h_n);
disp("Resultado Winograd:");
disp(y_winograd);

% generar arreglos enormes de 8 a 2^13
disp("Generando gráficos, espere un minuto...");
tamanos = 2.^(3:13);
% comparacion de tiempos
times_winograd = zeros(size(tamanos));
times_toeplitz = zeros(size(tamanos));

% se itera sobre las señales y se obtienen los dato para graficar
% apropiadamente
for i = 1:length(tamanos)
    n = tamanos(i);
    x = generarX(n);
    h = generarX(6);
    
    times_winograd(i) = timeit(@() winograd(x, h));
    times_toeplitz(i) = timeit(@() toep_convolucion(x, h));
end

% grafico de resultados
figure;
loglog(tamanos, times_winograd, 'b-o', tamanos, times_toeplitz, 'r-*');
xlabel("Tamaño de la señal x[n]");
ylabel("Tiempo de ejecución (s)");
legend("Winograd", "Toeplitz");
title("Comparación de tiempos de convolución");
grid on;

% funcion principal de Winograd para filtros de largo multiplo de 3
% division en chunks y luego se aplica por ventana
% cuidando el solapado
function y_n = winograd(x, h)
    N = length(x);
    M = length(h);
    % verificacion de largo de filtro
    if mod(M, 3) ~= 0
        error("La longitud del filtro h debe ser múltiplo de 3.");
    end
    
    % longitud del resultado de convolución
    L = N + M - 1;
    % preinicializar output
    y_n = zeros(1, L);
    
    % numero de chunks de tamaño 3 en el filtro
    num_chunks = M / 3;
    
    for i = 0:num_chunks-1
        % extraer chunk del filtro de largo 3 para aplicar ventana
        % se tiene que dar vuelta ???
        h_chunk = fliplr(h(3*i+1:3*i+3));
        
        % Calcular convolución F(2,3) para este fragmento
        conv_chunk = winograd_ventana(x, h_chunk);
        
        % Calcular shift y agregar al resultado (overlap-add)
        shift = 3*i;
        end_idx = min(L, shift + length(conv_chunk));
        y_n(shift+1:end_idx) = y_n(shift+1:end_idx) + conv_chunk(1:end_idx-shift);
    end
end

% Implementación de F(2,3) para filtros de exactamente 3 elementos
function y_n = winograd_ventana(x, h)
    N = length(x);
    
    % Verificar que el filtro sea de largo 3
    if length(h) ~= 3
        error("winograd_ventana requiere un filtro de largo 3");
    end
    
    % calcular longitud del resultado
    L = N + 2;
    
    % se hace padding a la señal de entrada
    x_pad = [zeros(1, 2), x, zeros(1, 3)];
    
    % preinicializar resultado, si no tira error
    y_n = zeros(1, L);
    
    % Calcular pares de salida usando F(2,3)
    for i = 0:ceil(L/2)-1
        % indices para la ventana de 4 elementos
        idx = 2*i+1:2*i+4;
        
        % asegurar que no nos pasamos del padding
        if idx(end) > length(x_pad)
            break;
        end
        
        % Extraer ventana
        d = x_pad(idx);
        
        % algoritmo F(2,3) de Winograd
        % basado en las formulas entregadas en la tarea
        % se hacen por ventanas
        m1 = (d(1) - d(3)) * h(1);
        m2 = (d(2) + d(3)) * (h(1) + h(2) + h(3)) / 2;
        m3 = (d(3) - d(2)) * (h(1) - h(2) + h(3)) / 2;
        m4 = (d(2) - d(4)) * h(3);
        
        y0 = m1 + m2 + m3;
        y1 = m2 - m3 - m4;
        
        % agregar al resultado si los indices están dentro del rango
        if 2*i+1 <= L
            y_n(2*i+1) = y0;
        end
        if 2*i+2 <= L
            y_n(2*i+2) = y1;
        end
        % si no sencillamente se ignora lol
    end
end

% función para convolucion usando matriz de Toeplitz
function y_n = toep_convolucion(x, h)
    len_x = length(x);
    len_h = length(h);
    
    % crear matriz de Toeplitz
    col = [h(:); zeros(len_x - 1, 1)];
    row = [h(1), zeros(1, len_x - 1)];
    H = toeplitz(col, row);
    
    % realizar convolucion
    y_n = H(1:(len_h + len_x - 1), 1:len_x) * x(:);
    % se entrega en formato fila
    y_n = y_n';
end