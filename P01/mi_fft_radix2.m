function X = mi_fft_radix2(x)
% mi_fft_radix2 - Implementación educativa de una FFT Radix-2 (DIT)
% x: Vector de entrada de longitud N (debe ser potencia de 2, ej. 1024)

N = length(x);

% 1. Verificar que N sea potencia de 2
if bitand(N, N-1) ~= 0
    error('La longitud de la señal debe ser estrictamente una potencia de 2.');
end

% 2. Reordenamiento inicial por inversión de bits (Bit-Reversal)
bits = log2(N);
indices = 0:(N-1);
% Generamos las cadenas binarias, las invertimos y convertimos a índice MATLAB (+1)
rev_indices = bin2dec(fliplr(dec2bin(indices, bits))) + 1;
X = x(rev_indices);

% 3. Cálculo por etapas (Total: log2(N) etapas)
for stage = 1:bits
    m = 2^stage;        % Tamaño del bloque actual en esta etapa
    m2 = m / 2;         % Mitad del bloque (número de mariposas independientes por bloque)

    % Factor de giro base para esta etapa: W_m = e^(-j * 2*pi / m)
    W_m = exp(-1j * 2 * pi / m);

    % Iterar sobre cada bloque de la etapa actual
    for k = 1:m:N
        W = 1; % Inicializamos el factor de giro W^0 = 1 para el inicio del bloque

        % Iterar sobre cada operación mariposa dentro del bloque
        for j = 0:(m2-1)
            i1 = k + j;         % Índice superior de la mariposa
            i2 = k + j + m2;    % Índice inferior de la mariposa

            % Extracción de valores
            u = X(i1);
            v = X(i2) * W;      % Multiplicación por el factor de giro (Twiddle Factor)

            % Operación mariposa fundamental
            X(i1) = u + v;      % Rama superior
            X(i2) = u - v;      % Rama inferior

            % Actualizar el factor de giro para la siguiente mariposa del bloque
            W = W * W_m;
        end
    end
end
end
