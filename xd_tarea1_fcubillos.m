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
y_matlab = conv(x_n, h_n);
disp('Resultado conv():');
disp(y_matlab);

% Convolución usando Winograd F(2,3)
disp("Convolución Winograd:");
y_winograd = winograd(x_n, h_n);
disp('Resultado Winograd:');
disp(y_winograd);

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
    f_win = @() winograd(x, h);
    times_winograd(i) = timeit(f_win);
    
    % Medir Toeplitz
    f_toep = @() toep_convolucion(x, h);
    times_toeplitz(i) = timeit(f_toep);
end

%Graficar resultados
figure;
loglog(sizes, times_winograd, 'b-o', sizes, times_toeplitz, 'r-*');
xlabel('Tamaño de la señal x[n]');
ylabel('Tiempo de ejecución (s)');
legend('Winograd', 'Toeplitz');
title('Comparación de tiempos de convolución');
grid on;

% --- Ensure this is the function being called ---
function y_n = winograd(x, h)
    N = length(x);
    M = length(h);

    if mod(M, 3) ~= 0
        padding_needed = 3 - mod(M, 3);
        h = [h, zeros(1, padding_needed)];
        M = length(h);
        disp('Advertencia: La longitud del filtro h no era múltiplo de 3. Se ha añadido padding.');
    end

    L = N + length(h) - 1; % Use original length(h) for L? No, use padded length for consistency? Let's use original N + original M - 1.
    L_correct = N + (M - (3-mod(length(h), 3)) ) -1 ; % L based on original M
    if mod(length(h),3) == 0 % M is the original length if multiple of 3
        L_correct = N + M -1;
    end

    y_n = zeros(1, L_correct); % Initialize output with CORRECT final length
    num_h_chunks = M / 3; % M is now padded length

    for i = 0:num_h_chunks-1
        % 1. Extract Filter Chunk (from potentially padded h)
        start_idx_h = 3*i + 1;
        end_idx_h = start_idx_h + 2;
        h_chunk = h(start_idx_h:end_idx_h);
        h_chunk = fliplr(h_chunk);

        % 2. Calculate convolution for this chunk
        % Assuming winograd_ventana is correct for 3-tap filter
        conv_i = winograd_ventana(x, h_chunk); % length(conv_i) = N + 3 - 1 = N + 2

        % 3. Overlap-Add with correct shift and index handling
        shift = 3*i; % Shift for filter chunk 'i'

        % Determine indices in y_n where conv_i should be added
        y_start_index = shift + 1;
        y_end_index   = shift + length(conv_i); % Tentative end index

        % Determine indices in conv_i that are valid to add
        conv_start_index = 1;
        conv_end_index = length(conv_i);

        % Adjust if y_end_index goes beyond the final length L_correct
        if y_end_index > L_correct
            % How many elements are out of bounds?
            overhang = y_end_index - L_correct;
            % Reduce the end index for y_n and the effective length from conv_i
            y_end_index = L_correct;
            conv_end_index = conv_end_index - overhang;
        end

        % Ensure we have a valid range to add
        if y_start_index <= y_end_index && conv_start_index <= conv_end_index
             y_indices_to_update = y_start_index:y_end_index;
             conv_indices_to_use = conv_start_index:conv_end_index;

             % Check if lengths match before adding (should always match here)
             if length(y_indices_to_update) == length(conv_indices_to_use)
                 y_n(y_indices_to_update) = y_n(y_indices_to_update) + conv_i(conv_indices_to_use);
             else
                 error('Internal error: Index length mismatch during overlap-add.');
             end
        end
        % Otherwise, this chunk's contribution falls entirely outside L_correct range
    end
     % Final result should already be size L_correct
     % y_n = y_n(1:L_correct); % This line is redundant if y_n initialized correctly
end

% --- Use the SAME winograd_ventana function provided in the previous answer ---
% --- Make sure it's exactly this one: ---
function y_n = winograd_ventana(x, h)
    N = length(x);
    M = length(h); % Should be 3
    if M ~= 3
        error('winograd_ventana expects h of length 3');
    end

    expected_length = N + M - 1;
    y_acum = [];

    % Pad x for window processing. Needs M-1 zeros before and >=M after for last window.
    pad_before = M-1;
    pad_after = M; % Sufficient padding for the last window
    x_padded = [zeros(1, pad_before), x, zeros(1, pad_after)];

    num_output_pairs = ceil(expected_length / 2);

    for i = 0:num_output_pairs-1
        % Input window d = [d0, d1, d2, d3] corresponding to output y0, y1
        % Indices relative to start of x_padded:
        idx_start = 2*i + 1;
        window_indices = idx_start : idx_start + 3;

        % Extract window safely, padding with zeros if indices go out of bounds
        window = zeros(1, 4);
        valid_mask = (window_indices >= 1) & (window_indices <= length(x_padded));
        actual_indices_in_padded = window_indices(valid_mask);
        % Ensure assignment matches window size if mask is partial
        window(valid_mask(1:length(actual_indices_in_padded))) = x_padded(actual_indices_in_padded);


        % Algoritmo F(2,3) de winograd (g0=h(1), g1=h(2), g2=h(3))
        g0 = h(1); g1 = h(2); g2 = h(3);
        d0 = window(1); d1 = window(2); d2 = window(3); d3 = window(4);

        m1 = (d0 - d2) * g0; % m0 in standard notation
        m2 = (d1 + d2) * (g0 + g1 + g2) / 2; % m1
        m3 = (d2 - d1) * (g0 - g1 + g2) / 2; % m2
        m4 = (d1 - d3) * g2; % m3

        y0 = m1 + m2 + m3;
        y1 = m2 - m3 - m4;

        % Resultado acumulado
        y_acum = [y_acum, y0, y1];
    end

    % Trim to the correct length
    if length(y_acum) >= expected_length
         y_n = y_acum(1:expected_length);
    else
         % Pad if somehow too short (shouldn't happen with ceil)
         y_n = [y_acum, zeros(1, expected_length - length(y_acum))];
    end
end

function y_n = toep_convolucion(x, h)
    % Computes convolution y = x * h using Toeplitz matrix formulation.
    len_x = length(x);
    len_h = length(h);

    % Create the first column for the Toeplitz matrix
    col = [h(:); zeros(len_x - 1, 1)];

    % Create the first row for the Toeplitz matrix
    row = [h(1), zeros(1, len_x - 1)];

    % Construct the full Toeplitz matrix
    H_full = toeplitz(col, row);

    % Extract the submatrix relevant for linear convolution
    H = H_full(1:(len_h + len_x - 1), 1:len_x);

    % Perform convolution via matrix multiplication
    y_n = H * x(:); % Ensure x is a column vector
end