clearvars
close all
clc

% Opciones del usuario
disp('A: Ingresar señales.')
disp('B: Señales aleatorias.')
opcion = input("Elige una opción: ", "s");

% Función para generar señales aleatorias
generarX = @(largo) -1 + 2*rand(1, largo);

if opcion == 'A'
    x_n = str2num(input("Ingrese la primera señal, separada por comas: ", "s"));
    h_n = str2num(input("Ingrese la respuesta al impulso, separada por comas: ", "s"));
else
    disp("Generando señales...");
    x_n = generarX(8); % Señal de largo 8
    h_n = generarX(6); % Filtro de largo 6 (múltiplo de 3)
    disp('x_n:'); disp(x_n);
    disp('h_n:'); disp(h_n);
end

% Convolución convencional con MATLAB
disp("Convolución de Matlab:");
tic
y_matlab = conv(x_n, h_n);
toc
disp('Resultado conv():');
disp(y_matlab);

% Convolución usando Winograd F(2,3)
disp("Convolución Winograd:");
tic
y_winograd = winograd_f23_1d_generalized(x_n, h_n);
toc
disp('Resultado Winograd:');
disp(y_winograd);

% Verificar exactitud
error = max(abs(y_matlab - y_winograd));
disp(['Error máximo: ', num2str(error)]);

% Convolución usando matriz Toeplitz
disp("Convolución con matriz Toeplitz:");
tic
y_toeplitz = toep_convolucion(x_n, h_n)';
toc
disp('Resultado Toeplitz:');
disp(y_toeplitz);

% Medición de tiempos para diferentes tamaños
sizes = 2.^(3:13);
times_winograd = zeros(size(sizes));
times_toeplitz = zeros(size(sizes));

for i = 1:length(sizes)
    n = sizes(i);
    x = generarX(n);
    h = generarX(6); % Filtro de largo 6
    
    % Medir Winograd
    f_win = @() winograd_f23_1d_generalized(x, h);
    times_winograd(i) = timeit(f_win);
    
    % Medir Toeplitz
    f_toep = @() toep_convolucion(x, h);
    times_toeplitz(i) = timeit(f_toep);
end

% Graficar resultados
figure;
loglog(sizes, times_winograd, 'b-o', sizes, times_toeplitz, 'r-*');
xlabel('Tamaño de la señal x[n]');
ylabel('Tiempo de ejecución (s)');
legend('Winograd', 'Toeplitz');
title('Comparación de tiempos de convolución');
grid on;



function y_n = toep_convolucion(x, h)
    len_x = length(x);
    len_h = length(h);
    
    % Crear columna base para Toeplitz (h seguido de ceros para alcanzar N+M-1)
    col = [h(:); zeros(len_x - 1, 1)]; 
    
    % Crear matriz Toeplitz con desplazamientos correctos
    H = toeplitz(col, [h(1), zeros(1, len_x - 1)]);
    
    % Asegurar dimensiones correctas (N + M - 1 x N)
    H = H(1:len_h + len_x - 1, 1:len_x);
    
    % Calcular convolución
    y_n = H * x(:);
end