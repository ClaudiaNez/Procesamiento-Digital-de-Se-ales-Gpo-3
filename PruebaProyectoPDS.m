% Parámetros de prueba
N = 1024;
fs = 1000;                  % Frecuencia de muestreo (1 kHz)
t = (0:N-1) / fs;           % Vector de tiempo

% Crear una señal de prueba: suma de dos senoides (ej. 50 Hz y 120 Hz)
x = sin(2 * pi * 50 * t) + 0.5 * sin(2 * pi * 120 * t);

% 1. Ejecutar nuestra implementación personalizada
tic;
X_custom = mi_fft_radix2(x);
tiempo_custom = toc;

% 2. Ejecutar la FFT nativa de MATLAB para comparar
tic;
X_builtin = fft(x);
tiempo_builtin = toc;

% 3. Calcular el error máximo absoluto para verificar precisión
error_max = max(abs(X_custom - X_builtin));

fprintf('--- RESULTADOS DE LA PRUEBA (N = 1024) ---\n');
fprintf('Tiempo función propia: %.6f segundos\n', tiempo_custom);
fprintf('Tiempo función nativa:  %.6f segundos\n', tiempo_builtin);
fprintf('Error máximo absoluto:  %.2e (debe ser cercano a cero)\n', error_max);

% 4. Gráfica comparativa de magnitudes
f = (0:(N/2)) * (fs / N);
figure;
subplot(2,1,1);
plot(f, abs(X_custom(1:N/2+1)), 'LineWidth', 1.5);
title('Espectro de Amplitud obtenido con la FFT Propia (Radix-2)');
xlabel('Frecuencia (Hz)'); ylabel('|X(f)|');
grid on;

subplot(2,1,2);
plot(f, abs(X_builtin(1:N/2+1)), 'r--', 'LineWidth', 1.5);
title('Espectro de Amplitud obtenido con la FFT Nativa de MATLAB');
xlabel('Frecuencia (Hz)'); ylabel('|X(f)|');
grid on;