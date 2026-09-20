% =========================================================================
% Universidad Nacional Autónoma de México
% Asigantura: Procesamiento Digital de Señales
% TEMA: Ejercicios de Señales en Tiempo Continuo y Discreto
% 
% INTEGRANTES:
%   - Cruces Nuñez Claudia Meera (320277038)
%   - Téllez Olmedo Eduardo Iván (321134873)
% =========================================================================
%% Limpieza del entorno
clear; clc; close all;

%% ------------------------------------------------------------------------
% EJERCICIO 1: Visualización de las señales
% ------------------------------------------------------------------------
figure('Name', 'Ejercicio 1');

% a) x(t) = sin(2*pi*10*t)
t = -0.2:0.001:0.2; % Vector de tiempo continuo (-0.2 a 0.2 segundos)
xa = sin(2*pi*10*t);
subplot(3,1,1);
plot(t, xa, 'LineWidth', 1.5);
title('1a) x(t) = sin(2\pi \cdot 10 t) [Continua, Periódica, No Causal]');
xlabel('Tiempo t (s)'); ylabel('Amplitud'); grid on;

% b) x[n] = u[n-3]
n = -5:10; % Vector de tiempo discreto (-5 a 10)
xb = double(n >= 3); % Función escalón desplazada 3 unidades a la derecha
subplot(3,1,2);
stem(n, xb, 'filled', 'LineWidth', 1.5);
title('1b) x[n] = u[n-3] [Discreta, Aperiódica, Causal]');
xlabel('Muestras n'); ylabel('Amplitud'); grid on;

% c) x(t) = e^(-t)*u(t)
tc = -1:0.01:4;
xc = exp(-tc) .* (tc >= 0); % Multiplicación por u(t)
subplot(3,1,3);
plot(tc, xc, 'LineWidth', 1.5);
title('1c) x(t) = e^{-t} u(t) [Continua, Aperiódica, Causal]');
xlabel('Tiempo t (s)'); ylabel('Amplitud'); grid on;

%% ------------------------------------------------------------------------
% EJERCICIO 2: Periodicidad de x[n] = cos((5*pi/6)*n)
% ------------------------------------------------------------------------
figure('Name', 'Ejercicio 2');
n2 = 0:35; % Mostramos 3 periodos completos (3 * 12 = 36 muestras)
x2 = cos((5*pi/6) * n2);

stem(n2, x2, 'filled', 'LineWidth', 1.5);
title('Ejercicio 2: x[n] = cos((5\pi/6)n) - Período N = 12');
xlabel('Muestras n'); ylabel('Amplitud'); grid on;

%% ------------------------------------------------------------------------
% EJERCICIO 3: Muestreo y Aliasing
% ------------------------------------------------------------------------
figure('Name', 'Ejercicio 3');

f0 = 100;                 % Frecuencia original de la señal (200*pi / 2*pi)
fs = 150;                 % Frecuencia de muestreo dada
f_app = abs(f0 - fs);     % Frecuencia aparente por aliasing (50 Hz)

t_cont = 0:0.0001:0.03;   % Tiempo denso para graficar señal continua
x_original = cos(2*pi*f0*t_cont);
x_aliased = cos(2*pi*f_app*t_cont);

% Muestras tomadas a fs = 150 Hz
n_sample = 0:(0.03*fs);
t_sample = n_sample / fs;
x_sampled = cos(2*pi*f0*t_sample);

plot(t_cont, x_original, 'b--', 'LineWidth', 1.2); hold on;
plot(t_cont, x_aliased, 'r', 'LineWidth', 1.5);
stem(t_sample, x_sampled, 'k', 'filled', 'LineWidth', 1.5);
hold off;

title('Ejercicio 3: Muestreo con Aliasing (f_0 = 100 Hz, f_s = 150 Hz \rightarrow f_{app} = 50 Hz)');
xlabel('Tiempo (s)'); ylabel('Amplitud');
legend('Señal Original (100 Hz)', 'Señal Percibida (50 Hz)', 'Muestras Tomadas');
grid on;

%% ------------------------------------------------------------------------
% EJERCICIO 4: Gráficas de señales básicas para -5 <= n <= 5
% ------------------------------------------------------------------------
figure('Name', 'Ejercicio 4');
n4 = -5:5;

% a) Escalón unitario u[n]
u = double(n4 >= 0); %La señal inicia en 0
subplot(3,1,1);
stem(n4, u, 'filled', 'LineWidth', 1.5);
title('4a) Escalón Unitario u[n]');
xlabel('n'); ylabel('Amplitud'); grid on;

% b) Impulso unitario delta[n]
delta = double(n4 == 0); %La señal inicia en 0
subplot(3,1,2);
stem(n4, delta, 'filled', 'LineWidth', 1.5);
title('4b) Impulso Unitario \delta[n]');
xlabel('n'); ylabel('Amplitud'); grid on;

% c) Rampa unitaria r[n] = n * u[n]
r = n4 .* (n4 >= 0); %La señal n4 inicia 0
subplot(3,1,3);
stem(n4, r, 'filled', 'LineWidth', 1.5);
title('4c) Rampa Unitaria r[n] = n \cdot u[n]');
xlabel('n'); ylabel('Amplitud'); grid on;

%% ------------------------------------------------------------------------
% EJERCICIOS 5, 6 y 7: Transformaciones de tiempo sobre x[n] = {1, 2, 3, 4}
% ------------------------------------------------------------------------
% Definimos la señal original y su soporte de índices
x_val = [1, 2, 3, 4];
n_x = 0:3; % x[0]=1, x[1]=2, x[2]=3, x[3]=4

figure('Name', 'Ejercicios 5, 6 y 7 - Transformaciones de tiempo');

% Señal Original x[n]
subplot(3,3,1);
stem(n_x, x_val, 'filled', 'LineWidth', 1.5);
title('Original: x[n]');
xlabel('n'); xlim([-4 6]); grid on;

% Ejercicio 5: Inversión temporal x[-n]
% Propiedad: x[-n] se obtiene invirtiendo el eje de tiempo n -> -n
n_5 = -fliplr(n_x);
x_5 = fliplr(x_val);
subplot(3,3,2);
stem(n_5, x_5, 'filled', 'r', 'LineWidth', 1.5);
title('5) x[-n]');
xlabel('n'); xlim([-4 6]); grid on;

% Ejercicio 6a: Desplazamiento a la derecha x[n-2]
% Propiedad: x[n - k] desplaza la señal k unidades a la derecha
n_6a = n_x + 2;
subplot(3,3,4);
stem(n_6a, x_val, 'filled', 'g', 'LineWidth', 1.5);
title('6a) x[n-2]');
xlabel('n'); xlim([-4 6]); grid on;

% Ejercicio 6b: Desplazamiento a la izquierda x[n+1]
% Propiedad: x[n + k] desplaza la señal k unidades a la izquierda
n_6b = n_x - 1;
subplot(3,3,5);
stem(n_6b, x_val, 'filled', 'g', 'LineWidth', 1.5);
title('6b) x[n+1]');
xlabel('n'); xlim([-4 6]); grid on;

% Ejercicio 7a: Inversión y desplazamiento x[-n+2] = x[-(n-2)]
% Método de sustitución de índice: n_orig = -n_nuevo + 2 => n_nuevo = 2 - n_orig
n_7a = 2 - n_x;
% Reordenamos para graficar en orden creciente de tiempo
[n_7a_sorted, idx_a] = sort(n_7a);
x_7a_sorted = x_val(idx_a);

subplot(3,3,7);
stem(n_7a_sorted, x_7a_sorted, 'filled', 'm', 'LineWidth', 1.5);
title('7a) x[-n+2]');
xlabel('n'); xlim([-4 6]); grid on;

% Ejercicio 7b: Inversión y desplazamiento x[3-n] = x[-(n-3)]
% Método de sustitución de índice: n_orig = 3 - n_nuevo => n_nuevo = 3 - n_orig
n_7b = 3 - n_x;
[n_7b_sorted, idx_b] = sort(n_7b);
x_7b_sorted = x_val(idx_b);

subplot(3,3,8);
stem(n_7b_sorted, x_7b_sorted, 'filled', 'm', 'LineWidth', 1.5);
title('7b) x[3-n]');
xlabel('n'); xlim([-4 6]); grid on;
