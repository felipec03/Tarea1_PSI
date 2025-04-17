clearvars
close all
clc

% tomar input de usuario
disp('A: Ingresar señales.')
disp('B : Señales aleatorias.')

opcion = input("Elige una opción: ", "s");

function x_senal = generarX(largo)
    largo = largo -1;
    mnval = -1;
    mxval = 1;
    val = mnval + rand*(mxval-mnval);
    x_dumb = [val];
    for i = 1:largo
        val = mnval + rand*(mxval-mnval);
        x_dumb = [x_dumb, val];
    end
    x_senal = x_dumb;
end

if opcion == 'A'
    x_n = input("Ingrese la primera señal, separada por comas: ","s");
    h_n = input("Ingrese la respuesta al impulso, separada por comas: ","s");
else
    % implementacion de señales aleatorias entre -1 y 1
    %x_n = randi([-1, 1],1, 100000);
    disp("Generando señal...");
    tic
    %x_n = generarX(8);
    x_n = [1,0,1,0,1,0,1,0];
    h_n = [0,0,1,0,0,0];
    toc
end

% prueba convencional
disp("Convolución de Matlab");
tic
y_n = conv(x_n, h_n);
disp(y_n);
toc

% usando la matriz de toeplitz
function y_n = toep_convolucion(x, h)
    n = length(x) + length(h) -1;
    h_n_paddeado = [h, zeros(1, n - length(h))];
    disp(h_n_paddeado);
    toeplitz_matriz = toeplitz([h_n_paddeado(1) fliplr(h_n_paddeado(2:n))], h);
    display(toeplitz_matriz);
    %y_n = toeplitz_matriz * x;

    %h_mayusc = transpose(toeplitz_matriz);
    %disp(h_mayusc);

    %y_n = h_mayusc * x;
    y_n = 0;
end
y_n = toep_convolucion(x_n, h_n);
%disp(y_n);