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
y_winograd = winograd_conv(x_n, h_n);
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
    f_win = @() winograd_conv(x, h);
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

function y = winograd_conv(x, h)
    % Verificar que h sea múltiplo de 3
    assert(mod(length(h), 3) == 0, 'h debe ser múltiplo de 3');
    
    % Dividir h en bloques de 3
    num_blocks = length(h) / 3;
    y = zeros(1, length(x) + length(h) - 1);
    
    % Procesar cada bloque
    for k = 0:num_blocks-1
        h_block = h(k*3 + 1 : (k+1)*3);
        y_block = winograd_block_conv(x, h_block);
        % Desplazar y sumar
        y(1:length(y_block) + k*3) = y(1:length(y_block) + k*3) + [y_block, zeros(1, k*3)];
    end
end

function y_block = winograd_block_conv(x, h_block)
    % Parámetros exactos Winograd F(2,3)
    Bt = [1, 0, -1, 0; 
          0, 1, 1, 0; 
          0, -1, 1, 0; 
          0, 1, 0, -1];
    G = [1, 0, 0; 
         1/2, 1/2, 1/2; 
         1/2, -1/2, 1/2; 
         0, 0, 1];
    At = [1, 1, 1, 0; 
          0, 1, -1, -1];
    
    Lx = length(x);
    y_block = zeros(1, Lx + 2); % Lx + 3 -1
    
    % Procesar cada bloque con zero-padding
    for i = 1:2:Lx
        d = x(i:min(i+3, Lx));
        d_padded = [d, zeros(1, 4 - length(d))];
        
        U = G * h_block(:);
        V = Bt * d_padded(:);
        M = U .* V;
        y_win = At * M;
        
        % Escalar y acumular
        y_block(i:i+1) = y_block(i:i+1) + (y_win / 2)';
    end
    % Recortar ceros finales
    y_block = y_block(1:Lx + 2);
end

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