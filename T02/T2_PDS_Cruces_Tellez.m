%% ========================================================
%  Tarea 2: DSP Exercise Set 2 - Correlation & Cross-Correlation
%  CRUCES NUÑEZ CLAUDIA MEERA					 
%  TELLEZ OLMEDO EDUARDO IVAN					
%  ========================================================
clear; clc; close all;

%% ===== EJERCICIO 1: Autocorrelación de una señal discreta =====
% x[n] = {1, 2, 3, 2, 1},  n = 0..4
x1 = [1 2 3 2 1];
n1 = 0:length(x1)-1;
N1 = length(x1);

% (a) y (b) Cálculo de la autocorrelación (y = x)
[Rxx, kRxx] = corr_pp(x1, x1);
fprintf('--- Ejercicio 1: Resultado ---\n');
fprintf('k      :'); fprintf('%5d', kRxx); fprintf('\n');
fprintf('Rxx[k] :'); fprintf('%5d', Rxx);  fprintf('\n\n');

%% Desglose elemento a elemento (procedimiento manual)
fprintf('--- Desarrollo del cálculo ---\n');
fprintf('n          :');  fprintf('%5d', n1);  fprintf('\n');
fprintf('x[n]       :');  fprintf('%5d', x1);  fprintf('\n');
fprintf('%s\n', repmat('-', 1, 70));
for k = kRxx
    % Construir x[n-k] dentro de la ventana n = 0..N1-1 (ceros fuera)
    xd = zeros(1, N1);
    for n = 0:N1-1
        if (n-k >= 0) && (n-k <= N1-1)
            xd(n+1) = x1(n-k+1);
        end
    end
    
    prod = x1 .* xd;   % productos columna a columna
    
    fprintf('k=%2d x[n-k]:', k);   fprintf('%5d', xd);
    fprintf('   prod:');           fprintf('%4d', prod);
    fprintf('   suma = %3d\n', sum(prod));
end
fprintf('\n');

%% Resultados de los incisos (d) y (e)
[maxRxx, idx] = max(Rxx);
fprintf('(d) Máximo: Rxx[%d] = %d\n', kRxx(idx), maxRxx);
fprintf('(e) ¿Simetría par (Rxx[k] = Rxx[-k])? %d\n\n', isequal(Rxx, fliplr(Rxx)));

%% (c) Gráficas
figure('Name','Ejercicio 1','Position',[100 100 700 800]);
subplot(3,1,1);
stem(n1, x1, 'filled', 'LineWidth', 1.2);
xlabel('n [muestras]'); ylabel('x[n]');
title('Señal x[n]'); grid on; ylim([0 4]);

subplot(3,1,2);
stem(kRxx, Rxx, 'filled', 'LineWidth', 1.2, 'Color', [0.20 0.40 0.70]);
xlabel('Lag k [muestras]'); ylabel('R_{xx}[k]');
title('Autocorrelación R_{xx}[k]'); grid on; ylim([0 21]);
xline(0, '--', 'Color', [0.5 0.5 0.5]);

subplot(3,1,3); hold on;
colores  = lines(numel(kRxx));
leyendas = cell(1, numel(kRxx));
for i = 1:numel(kRxx)
    k = kRxx(i);
    n_full = k : k+N1-1;                  % soporte de x[n-k]
    stem(n_full, x1, 'filled', 'LineWidth', 1, 'Color', colores(i,:));
    leyendas{i} = sprintf('x[n-(%d)]', k);
end
hold off;
xlabel('n [muestras]'); ylabel('Amplitud');
title('Señales desplazadas x[n-k] utilizadas en el cálculo');
grid on; ylim([0 4]);
legend(leyendas, 'Location','eastoutside','FontSize',7);

%% ===== EJERCICIO 2: Correlación cruzada y desplazamiento temporal =====
x2 = [0 1 2 3 2 1 0];
y2 = [0 0 0 1 2 3 2 1 0];
nx2 = 0:length(x2)-1;
ny2 = 0:length(y2)-1;
Nx2 = length(x2);
Ny2 = length(y2);

[rxy, krxy] = corr_pp(x2, y2);
fprintf('--- Ejercicio 2: Resultado ---\n');
fprintf('k      :'); fprintf('%5d', krxy); fprintf('\n');
fprintf('rxy[k] :'); fprintf('%5d', rxy);  fprintf('\n\n');

fprintf('--- Desarrollo del cálculo ---\n');
fprintf('n            :');  fprintf('%4d', ny2);  fprintf('\n');
fprintf('y[n]         :');  fprintf('%4d', y2);   fprintf('\n');
fprintf('%s\n', repmat('-', 1, 80));
for k = krxy
    xd = zeros(1, Ny2);
    for n = 0:Ny2-1
        if (n+k >= 0) && (n+k <= Nx2-1)
            xd(n+1) = x2(n+k+1);
        end
    end
    prod = xd .* y2;
    fprintf('k=%3d x[n+k]:', k);  fprintf('%4d', xd);
    fprintf('  prod:');           fprintf('%4d', prod);
    fprintf('  suma = %3d\n', dot(xd, y2));
end
fprintf('\n');

[maxr, idx] = max(rxy);
kmax = krxy(idx);
fprintf('(d) Máximo: rxy[%d] = %d\n', kmax, maxr);
fprintf('(e) y[n] = x[n%+d]  -->  y está retrasada %d muestras respecto a x\n\n', kmax, abs(kmax));

figure('Name','Ejercicio 2','Position',[150 50 700 900]);
subplot(4,1,1);
stem(nx2, x2, 'filled', 'LineWidth', 1.2);
xlabel('n [muestras]'); ylabel('x[n]'); title('Señal x[n]'); grid on;
xlim([-0.5 8.5]); ylim([0 3.5]);

subplot(4,1,2);
stem(ny2, y2, 'filled', 'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
xlabel('n [muestras]'); ylabel('y[n]'); title('Señal y[n]'); grid on;
xlim([-0.5 8.5]); ylim([0 3.5]);

subplot(4,1,3);
stem(krxy, rxy, 'filled', 'LineWidth', 1.2, 'Color', [0.20 0.40 0.70]);
xlabel('Lag k [muestras]'); ylabel('r_{xy}[k]');
title('Correlación cruzada r_{xy}[k]'); grid on; ylim([0 21]);
xline(kmax, '--r', sprintf('k = %d', kmax), 'LabelVerticalAlignment','bottom');

subplot(4,1,4); hold on;
colores2  = lines(numel(krxy));
leyendas2 = cell(1, numel(krxy));
for i = 1:numel(krxy)
    k = krxy(i);
    n_full = (0:Nx2-1) - k;
    stem(n_full, x2, 'filled', 'LineWidth', 1, 'Color', colores2(i,:));
    leyendas2{i} = sprintf('x[n+(%d)]', k);
end
hold off;
xlabel('n [muestras]'); ylabel('Amplitud');
title('Señales desplazadas x[n+k] utilizadas en el cálculo');
grid on; ylim([0 3.5]);
legend(leyendas2, 'Location','eastoutside','FontSize',6);

%% ===== EJERCICIO 3: Detección de un patrón conocido =====
x3 = [0 0 1 -1 2 1 -1 0 0];
p3 = [1 -1 2 1 -1];
nx3 = 0:length(x3)-1;
np3 = 0:length(p3)-1;
Nx3 = length(x3);
Np3 = length(p3);

[rxp, krxp] = corr_pp(x3, p3);
fprintf('--- Ejercicio 3: Resultado ---\n');
fprintf('k      :'); fprintf('%5d', krxp); fprintf('\n');
fprintf('rxp[k] :'); fprintf('%5d', rxp);  fprintf('\n\n');

fprintf('--- Desarrollo del cálculo ---\n');
fprintf('n            :');  fprintf('%4d', np3);  fprintf('\n');
fprintf('p[n]         :');  fprintf('%4d', p3);   fprintf('\n');
fprintf('%s\n', repmat('-', 1, 70));
for k = krxp
    xd = zeros(1, Np3);
    for n = 0:Np3-1
        if (n+k >= 0) && (n+k <= Nx3-1)
            xd(n+1) = x3(n+k+1);
        end
    end
    prod = xd .* p3;
    fprintf('k=%3d x[n+k]:', k);  fprintf('%4d', xd);
    fprintf('  prod:');           fprintf('%4d', prod);
    fprintf('  suma = %3d\n', dot(xd, p3));
end
fprintf('\n');

magnitudes = abs(rxp);
[maxmag, idx3] = max(magnitudes);
kmax3   = krxp(idx3);
valmax3 = rxp(idx3);
fprintf('(c) Máximo: rxp[%d] = %d\n', kmax3, valmax3);
fprintf('    El patrón inicia en la muestra n = %d\n', kmax3);
fprintf('    x[%d..%d] = ', kmax3, kmax3+Np3-1);
fprintf('%d ', x3(kmax3+1 : kmax3+Np3));
fprintf('  =  p[n]\n\n');

figure('Name','Ejercicio 3','Position',[200 50 700 900]);
subplot(4,1,1);
stem(nx3, x3, 'filled', 'LineWidth', 1.2);
xlabel('n [muestras]'); ylabel('x[n]'); title('Señal observada x[n]'); grid on;
xlim([-0.5 8.5]); ylim([-2 3]); yline(0, 'k-');

subplot(4,1,2);
stem(np3, p3, 'filled', 'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
xlabel('n [muestras]'); ylabel('p[n]'); title('Patrón buscado p[n]'); grid on;
xlim([-0.5 8.5]); ylim([-2 3]); yline(0, 'k-');

subplot(4,1,3);
stem(krxp, rxp, 'filled', 'LineWidth', 1.2, 'Color', [0.20 0.40 0.70]);
xlabel('Lag k [muestras]'); ylabel('r_{xp}[k]');
title('Correlación cruzada r_{xp}[k]'); grid on; ylim([-3 9]); yline(0, 'k-');
xline(kmax3, '--r', sprintf('k = %d', kmax3), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment','left');

subplot(4,1,4); hold on;
colores3  = lines(numel(krxp));
leyendas3 = cell(1, numel(krxp));
for i = 1:numel(krxp)
    k = krxp(i);
    n_full = (0:Nx3-1) - k;
    stem(n_full, x3, 'filled', 'LineWidth', 1, 'Color', colores3(i,:));
    leyendas3{i} = sprintf('x[n+(%d)]', k);
end
hold off;
xlabel('n [muestras]'); ylabel('Amplitud');
title('Señales desplazadas x[n+k] utilizadas en el cálculo');
grid on; ylim([-2 3]); yline(0, 'k-');
legend(leyendas3, 'Location','eastoutside','FontSize',6);

%% ===== EJERCICIO 4: Autocorrelación de una senoide =====
fs4 = 100;              % frecuencia de muestreo [Hz]
f0  = 5;                % frecuencia de la señal [Hz]
Ts4 = 1/fs4;            % periodo de muestreo [s]
t4  = 0:Ts4:2;          % vector de tiempo [s]
x4  = sin(2*pi*f0*t4);  % señal muestreada
N4  = length(x4);
fprintf('--- Ejercicio 4 ---\n');
fprintf('Número de muestras: N = %d\n', N4);

[Rxx4, k4] = corr_pp(x4, x4);

picos = [];
for i = 2:numel(Rxx4)-1
    if Rxx4(i) > Rxx4(i-1) && Rxx4(i) > Rxx4(i+1)
        picos(end+1) = k4(i);    
    end
end
fprintf('(c) Lags de los máximos: ');
fprintf('%d  ', picos); fprintf('\n');
separacion = diff(picos);
N0 = separacion(1);
fprintf('    Separación entre máximos: %d muestras\n', N0);

dt = N0 * Ts4;
fprintf('(d) Separación en tiempo: %.2f s\n', dt);

f_est = 1/dt;
fprintf('(e) Frecuencia estimada: 1/%.2f = %.1f Hz  (original: %d Hz)\n\n', dt, f_est, f0);

figure('Name','Ejercicio 4','Position',[250 100 700 600]);
subplot(2,1,1);
plot(t4, x4, 'LineWidth', 1.2);
xlabel('t [s]'); ylabel('x(t)'); title('Señal x(t) = sin(2\pi \cdot 5t)'); grid on;
xlim([0 2]); ylim([-1.2 1.2]);

subplot(2,1,2);
plot(k4, Rxx4, 'LineWidth', 1.2, 'Color', [0.20 0.40 0.70]); hold on;
plot(picos, Rxx4(ismember(k4,picos)), 'ro', 'MarkerFaceColor','r'); hold off;
xlabel('Lag k [muestras]'); ylabel('R_{xx}[k]');
title('Autocorrelación R_{xx}[k]  (máximos marcados en rojo)');
grid on; yline(0,'k-');

%% ===== EJERCICIO 5: Comparación de señales con distintas frecuencias =====
fs5 = 100;
Ts5 = 1/fs5;
t5  = 0:Ts5:2;
x5  = sin(2*pi*5*t5);
r1  = sin(2*pi*5*t5);
r2  = sin(2*pi*8*t5);

[rx_r1, k_r1] = corr_pp(x5, r1);
[rx_r2, k_r2] = corr_pp(x5, r2);

mag1 = abs(rx_r1); [max1, i1] = max(mag1);
mag2 = abs(rx_r2); [max2, i2] = max(mag2);
fprintf('--- Ejercicio 5 ---\n');
fprintf('(c) max|r_x,r1| = %.2f   en k = %d\n', max1, k_r1(i1));
fprintf('    max|r_x,r2| = %.2f   en k = %d\n', max2, k_r2(i2));
fprintf('    Razón entre máximos: %.1f\n', max1/max2);
fprintf('(d) La referencia r1 (5 Hz) da la respuesta más fuerte.\n\n');

figure('Name','Ejercicio 5 - Señales','Position',[300 100 700 600]);
subplot(3,1,1); plot(t5, x5, 'LineWidth', 1.2);
xlabel('t [s]'); ylabel('x(t)'); title('x(t) = sin(2\pi \cdot 5t)'); grid on; xlim([0 1]); ylim([-1.2 1.2]);

subplot(3,1,2); plot(t5, r1, 'LineWidth', 1.2, 'Color', [0.20 0.60 0.30]);
xlabel('t [s]'); ylabel('r_1(t)'); title('r_1(t) = sin(2\pi \cdot 5t)   (misma frecuencia)'); grid on; xlim([0 1]); ylim([-1.2 1.2]);

subplot(3,1,3); plot(t5, r2, 'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
xlabel('t [s]'); ylabel('r_2(t)'); title('r_2(t) = sin(2\pi \cdot 8t)   (frecuencia distinta)'); grid on; xlim([0 1]); ylim([-1.2 1.2]);

figure('Name','Ejercicio 5 - Correlaciones','Position',[350 100 700 500]);
subplot(2,1,1); plot(k_r1, rx_r1, 'LineWidth', 1.2, 'Color', [0.20 0.60 0.30]);
xlabel('Lag k [muestras]'); ylabel('r_{x,r1}[k]');
title(sprintf('Correlación con r_1 (5 Hz)   -   máximo = %.1f', max1)); grid on; yline(0,'k-'); ylim([-110 110]);

subplot(2,1,2); plot(k_r2, rx_r2, 'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
xlabel('Lag k [muestras]'); ylabel('r_{x,r2}[k]');
title(sprintf('Correlación con r_2 (8 Hz)   -   máximo = %.1f', max2)); grid on; yline(0,'k-'); ylim([-110 110]);

%% =========================================================================
%  Parámetros de Tiempo Continuo/Muestreado para Ejercicios 6, 7 y 10
% =========================================================================
fs = 100;              % Frecuencia de muestreo [Hz]
t  = 0 : 1/fs : 2;      % Vector de tiempo de 0 a 2 segundos

%% =========================================================================
% EJERCICIO 6: Encontrar una frecuencia oculta
% =========================================================================
fprintf('=== EJERCICIO 6: Frecuencia Oculta ===\n');
x6 = 2 * sin(2*pi*7*t);
f_candidatas6 = [3, 5, 7, 9, 11];
resp6 = zeros(size(f_candidatas6));

for i = 1:length(f_candidatas6)
    f_cand = f_candidatas6(i);
    ref_s = sin(2*pi*f_cand*t);
    ref_c = cos(2*pi*f_cand*t);
    
    I = dot(x6, ref_s);
    Q = dot(x6, ref_c);
    resp6(i) = sqrt(I^2 + Q^2);
end

figure('Name', 'Ejercicio 6 - Frecuencia Oculta', 'NumberTitle', 'off');
stem(f_candidatas6, resp6, 'filled', 'LineWidth', 1.5);
title('Detección de Frecuencia Oculta (x(t) = 2 sin(2\pi \cdot 7t))');
xlabel('Frecuencia Candidata (Hz)'); ylabel('Magnitud del Producto Punto'); grid on;

%% =========================================================================
% EJERCICIO 7: Detección de dos frecuencias
% =========================================================================
fprintf('=== EJERCICIO 7: Dos Frecuencias Dominantes ===\n');
x7 = 2*sin(2*pi*5*t) + 0.8*sin(2*pi*12*t);
f_candidatas7 = 2:20;
resp7 = zeros(size(f_candidatas7));

for i = 1:length(f_candidatas7)
    f_cand = f_candidatas7(i);
    ref_s = sin(2*pi*f_cand*t);
    ref_c = cos(2*pi*f_cand*t);
    
    I = dot(x7, ref_s);
    Q = dot(x7, ref_c);
    resp7(i) = sqrt(I^2 + Q^2);
end

figure('Name', 'Ejercicio 7 - Dos Frecuencias', 'NumberTitle', 'off');
stem(f_candidatas7, resp7, 'filled', 'LineWidth', 1.5, 'Color', [0.494 0.184 0.556]);
title('Respuesta del Detector para Dos Frecuencias (5 Hz y 12 Hz)');
xlabel('Frecuencia Candidata (Hz)'); ylabel('Magnitud de Correlación'); grid on;

%% =========================================================================
%  EJERCICIO 8: Fase y Detección Sinusoidal
% =========================================================================
fprintf('=== EJERCICIO 8: Efecto de la Fase en la Detección ===\n');

% Señales de entrada y referencia a 8 Hz
x1_8 = sin(2*pi*8*t);                   % Fase 0
x2_8 = sin(2*pi*8*t + pi/3);            % Fase pi/3
rs_8 = sin(2*pi*8*t);                   % Referencia Seno
rc_8 = cos(2*pi*8*t);                   % Referencia Coseno

% (a) Gráfica de las señales x1, x2 y rs
figure('Name', 'Ejercicio 8 - Señales y Fase', 'NumberTitle', 'off');
plot(t, x1_8, 'b', 'LineWidth', 1.2); hold on;
plot(t, x2_8, 'r--', 'LineWidth', 1.2);
plot(t, rs_8, 'k:', 'LineWidth', 1.5); hold off;
title('Ejercicio 8: Señales x_1(t), x_2(t) y referencia r_s(t)');
xlabel('Tiempo (s)'); ylabel('Amplitud'); grid on; xlim([0 0.5]);
legend('x_1(t) = sin(2\pi\cdot8t)', 'x_2(t) = sin(2\pi\cdot8t + \pi/3)', 'r_s(t) = sin(2\pi\cdot8t)');

% (b) Correlación de x1 con rs[cite: 3]
corr_x1_rs = dot(x1_8, rs_8);

% (c) Correlación de x2 con rs[cite: 3]
corr_x2_rs = dot(x2_8, rs_8);

% (d) Imprimir y comparar[cite: 3]
fprintf('(b) Correlación <x1, rs> = %.2f\n', corr_x1_rs);
fprintf('(c) Correlación <x2, rs> = %.2f\n', corr_x2_rs);
fprintf('(d) Explicación: Ambas señales tienen la misma frecuencia (8 Hz), pero x2 está desfasada.\n');
fprintf('    La correlación con seno depende de cos(fase), reduciendo la respuesta en cos(pi/3) = 0.5.\n');

% (e) Correlación con la referencia Coseno rc[cite: 3]
corr_x1_rc = dot(x1_8, rc_8);
corr_x2_rc = dot(x2_8, rc_8);
fprintf('(e) Correlación <x1, rc> = %.2f  |  <x2, rc> = %.2f\n', corr_x1_rc, corr_x2_rc);

% (f) Magnitud total usando cuadratura[cite: 3]
mag_x1 = sqrt(corr_x1_rs^2 + corr_x1_rc^2);
mag_x2 = sqrt(corr_x2_rs^2 + corr_x2_rc^2);
fprintf('(f) Magnitud combinada x1: %.2f  |  x2: %.2f (¡Iguales e independientes de la fase!)\n\n', mag_x1, mag_x2);

%% =========================================================================
%  EJERCICIO 9: Detección de Frecuencias en Ruido Blanco Gaussiano
% =========================================================================
fprintf('=== EJERCICIO 9: Detección con Ruido Blanco Gaussiano ===\n');

% Parámetros[cite: 3]
sigmas = [0.5, 2.0, 5.0];               % Diferentes niveles de ruido (inciso f)[cite: 3]
freqs9 = 1:20;                          % Candidatas de 1 a 20 Hz (inciso b)[cite: 3]
s_clean = sin(2*pi*6*t) + 0.7*sin(2*pi*15*t); % Señal limpia (6 Hz y 15 Hz)[cite: 3]

figure('Name', 'Ejercicio 9 - Detección en Ruido', 'NumberTitle', 'off', 'Position', [100 100 1000 700]);

for idx = 1:length(sigmas)
    sigma = sigmas(idx);

    % Generar ruido blanco gaussiano N(0, sigma^2)[cite: 3]
    rng(42);                            % Semilla para reproducibilidad
    w = sigma * randn(size(t));
    x9 = s_clean + w;                   % Señal ruidosa[cite: 3]

    % (c) Calcular la respuesta de correlación para cada frecuencia candidata[cite: 3]
    R_f = zeros(size(freqs9));
    for i = 1:length(freqs9)
        f = freqs9(i);
        r_s = sin(2*pi*f*t);
        r_c = cos(2*pi*f*t);

        I = dot(x9, r_s);
        Q = dot(x9, r_c);
        R_f(i) = sqrt(I^2 + Q^2);        % Magnitud en cuadratura[cite: 3]
    end

    % (a) Graficar señal limpia vs ruidosa[cite: 3]
    subplot(length(sigmas), 2, 2*idx - 1);
    plot(t, x9, 'Color', [0.85 0.33 0.10], 'LineWidth', 0.8); hold on;
    plot(t, s_clean, 'b', 'LineWidth', 1.5); hold off;
    title(sprintf('a) Señal en el Tiempo (\\sigma = %.1f)', sigma));
    xlabel('Tiempo (s)'); ylabel('Amplitud'); grid on; xlim([0 1]);
    legend('Ruidosa x[n]', 'Limpia s[n]');

    % (d) Graficar Magnitud vs Frecuencia[cite: 3]
    subplot(length(sigmas), 2, 2*idx);
    stem(freqs9, R_f, 'filled', 'LineWidth', 1.2, 'Color', [0.20 0.60 0.30]); hold on;
    xline(6, '--r', '6 Hz');
    xline(15, '--m', '15 Hz'); hold off;
    title(sprintf('d) Espectro de Correlación (\\sigma = %.1f)', sigma));
    xlabel('Frecuencia (Hz)'); ylabel('Magnitud'); grid on;
end

fprintf('(e) Con \\sigma = 0.5 y 2.0, los picos en 6 Hz y 15 Hz aún destacan claramente.\n');
fprintf('(f) Al aumentar \\sigma (ej. 5.0), el piso de ruido sube y dificulta la detección de 15 Hz.\n\n');

%% =========================================================================
% EJERCICIO 10: Detector de Frecuencias Robusto (Fase y Ruido)
% =========================================================================
fprintf('=== EJERCICIO 10: Detector con Cuadratura ===\n');
x10 = 1.5*sin(2*pi*4*t) + 0.7*sin(2*pi*9*t + pi/4) + 2*sin(2*pi*16*t - pi/3);
freqs10 = 1:20;
mags10 = zeros(size(freqs10));

for i = 1:length(freqs10)
    f = freqs10(i);
    r_sin = sin(2*pi*f*t);
    r_cos = cos(2*pi*f*t);
    I = dot(x10, r_sin);
    Q = dot(x10, r_cos);
    mags10(i) = sqrt(I^2 + Q^2);
end

figure('Name', 'Ejercicio 10 - Detector Completo', 'NumberTitle', 'off');
subplot(2,1,1); plot(t, x10, 'LineWidth', 1.2);
title('Señal Compuesta x(t)'); xlabel('Tiempo (s)'); grid on;
subplot(2,1,2); stem(freqs10, mags10, 'filled', 'LineWidth', 1.5, 'Color', [0.929 0.694 0.125]);
title('Respuesta del Detector de Frecuencias (Producto Punto Seno/Coseno)');
xlabel('Frecuencia (Hz)'); ylabel('Magnitud Detectada'); grid on;

fprintf('\n=== Todos los ejercicios han sido ejecutados correctamente ===\n');

% =========================================================================
% FUNCIONES LOCALES (Deben ubicarse al final del archivo)
% =========================================================================
function [r, lags] = corr_pp(x, y) 
    % CORR_PP  Correlación cruzada por producto punto.
    % Formula: r[k] = sum_n x[n+k]*y[n]
    x = x(:).';  y = y(:).';
    Nx = length(x);  Ny = length(y); 
    lags = -(Ny-1) : (Nx-1);     
    r = zeros(1, numel(lags));
    
    for i = 1:numel(lags)
        k = lags(i);
        n_ini = max(0, -k);
        n_fin = min(Ny-1, Nx-1-k);
        if n_ini <= n_fin
            n = n_ini:n_fin;
            r(i) = dot( x(n+k+1), y(n+1) );
        else
            r(i) = 0;          
        end
    end
end

