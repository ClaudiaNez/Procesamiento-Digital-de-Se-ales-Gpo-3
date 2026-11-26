clear; close all; clc;

%% Parametros
archA = 'PDS_voz_60Hz_4s.wav';            
archB = 'PDS_voz_fuente_simulada_4s.wav';   

% Filtro pasa altas, lo usaremos para quitar el tono de 60 Hz ---

fcA = 95;       % frecuencia de corte [Hz]
NA  = 1001;     % número de coeficientes 

% Filtro pasa bajas, lo usaremos para quitar el ruido arriba de la voz 
% ya que arriba de ~4 kHz casi todo es ruido.
fcB = 3500;     % [Hz] frecuencia de corte
NB  = 101;      % número de coeficientes

t1 = 3.80;  t2 = 3.95;      % tramo SIN voz 

[xA, fs]  = audioread(archA);
[xB, fsB] = audioread(archB);
assert(fs == fsB, 'Los audios deben tener la misma fs');

%% Diseño por ventanas

hA = zeros(NA, 1);   
hB = zeros(NB, 1);

for caso = 1:2
    % Seleccion de parametros segun el filtro a diseñar
    if caso == 1, N = NA; fc = fcA; %Longitud y frecuecnia corte del filtro A
        else, N = NB; fc = fcB;  %Longitud y frecuecnia corte del filtro B
    end

    M  = N - 1; %Orden del filtro
    wc = 2*pi*fc/fs; % frecuencia de corte
    n  = (0:M)'; %Vector de Indices
    w  = 0.5*(1 - cos(2*pi*(n+1)/(N+1)));   % ventana Hann de diseño para suavizar el filtro

    hd = zeros(N, 1); % Inicializacion de la respuesta impulsiva ideal
    for i = 1:N
        d = n(i) - M/2; % Distancia respecto al centro del filtro
        if d == 0 % Caso especial para evitar division entre cero
            hd(i) = wc/pi;
        else
            hd(i) = sin(wc*d) / (pi*d);
        end
    end
   % Aqui generamos cada filtro
    if caso == 1  % A pasa-altas 
        hd = -hd;
        hd(M/2 + 1) = hd(M/2 + 1) + 1;
        hA = hd .* w;
    else  % B Se conserva el pasa-bajas ideal y solo se ventanea
        hB = hd .* w;
    end
end

%% Ecuaciòn en diferencia y respuesta en frecuencia

Nfft = 2^16;   % Número de puntos para la FFT  
f    = (0:Nfft/2)' * fs/Nfft;
% Transformada rápida de Fourier de ambos filtros
HA   = fft(hA, Nfft);   HA = HA(1:Nfft/2+1);
HB   = fft(hB, Nfft);   HB = HB(1:Nfft/2+1);
% Funciones anónimas para calcular la ganancia en una frecuencia específica
ganA = @(fx) 20*log10(abs(sum(hA .* exp(-1j*2*pi*fx/fs*(0:NA-1)'))));  
ganB = @(fx) 20*log10(abs(sum(hB .* exp(-1j*2*pi*fx/fs*(0:NB-1)'))));

fprintf('Filtro pasa-altas, fc = %d Hz, N = %d \n', fcA, NA);
c = (NA-1)/2;                       % índice del centro (h(c+1) = h[M/2])
fprintf('H(z) = %+.2e %+.2e z^-1 + ... %+.5f z^-%d + ... %+.2e z^-%d\n', hA(1), hA(2), hA(c+1), c, hA(end), NA-1);
fprintf('y[n] = %+.2e x[n] %+.2e x[n-1] + ... %+.5f x[n-%d] + ... %+.2e x[n-%d]\n', hA(1), hA(2), hA(c+1), c, hA(end), NA-1);
fprintf('h[n] simétrica: max|h[n]-h[M-n]| = %.1e  -> fase lineal\n', max(abs(hA - flipud(hA))));
fprintf('Retardo = M/2 = %d muestras = %.2f ms\n', (NA-1)/2, (NA-1)/2/fs*1000);
fprintf('|H| en  60 Hz = %6.1f dB   (tono)\n', ganA(60));
fprintf('|H| en 120 Hz = %6.2f dB   |H| en 200 Hz = %6.2f dB\n\n', ganA(120), ganA(200));

fprintf('Filtro pasa-bajas, fc = %d Hz, N = %d \n', fcB, NB);
c = (NB-1)/2;                       % índice del centro (h(c+1) = h[M/2])
fprintf('H(z) = %+.2e %+.2e z^-1 + ... %+.5f z^-%d + ... %+.2e z^-%d\n', hB(1), hB(2), hB(c+1), c, hB(end), NB-1);
fprintf('y[n] = %+.2e x[n] %+.2e x[n-1] + ... %+.5f x[n-%d] + ... %+.2e x[n-%d]\n', hB(1), hB(2), hB(c+1), c, hB(end), NB-1);
fprintf('h[n] simétrica: max|h[n]-h[M-n]| = %.1e  -> fase lineal\n', max(abs(hB - flipud(hB))));
fprintf('Retardo = M/2 = %d muestras = %.2f ms\n', (NB-1)/2, (NB-1)/2/fs*1000);
fprintf('|H| en 1000 Hz = %6.2f dB   |H| en 3000 Hz = %6.2f dB\n', ganB(1000), ganB(3000));
fprintf('|H| en 4000 Hz = %6.1f dB   |H| en 5000 Hz = %6.1f dB\n\n', ganB(4000), ganB(5000));

%% Grafica de los ceros de la funciòn


zA = roots(hA);
zB = roots(hB);

% Se buscan los ceros localizados sobre el círculo unitario.
% Cada cero sobre el círculo unitario anula una frecuencia.
for caso = 1:2
    if caso == 1, z = zA; nom = 'A'; else, z = zB; nom = 'B'; end
    k  = abs(abs(z) - 1) < 1e-3 & imag(z) >= 0;     % Selección de los ceros con magnitud aproximadamente igual a 1
    fz = sort(angle(z(k)) * fs/(2*pi));   % Conversión del ángulo del cero a frecuencia
    fprintf('%s: %d ceros sobre el círculo unitario. Frecuencias anuladas [Hz]:\n', nom, 2*sum(k));
    fprintf('   %s\n', sprintf('%.0f  ', fz(1:min(8, end))));
end

th = linspace(0, 2*pi, 400);     % para dibujar el círculo unitario

for caso = 1:2
    if caso == 1
        h = hA; H = HA; z = zA; fc = fcA; nom = 'A: pasa-altas';
        fmaxZoom = 300;                         % zoom de |H| en bajas frecuencias
    else
        h = hB; H = HB; z = zB; fc = fcB; nom = 'B: pasa-bajas';
        fmaxZoom = fs/2;
    end
    M  = length(h) - 1;
    wc = 2*pi*fc/fs;

    figure(caso);
    
    % Respuesta al impulso h[n]
    subplot(2,3,1); stem(0:M, h, 'Marker', 'none'); grid on;
    xlabel('n'); title(['h[n]  ', nom]);
    % (A: pico de ~1 en M/2 = la delta, menos una sinc ancha y pequeña)

    % Magnitud completa del filtro en dB
    subplot(2,3,2); plot(f, 20*log10(abs(H) + eps)); grid on;
    xlim([0 fs/2]); ylim([-100 5]); xlabel('f [Hz]'); ylabel('dB'); title('|H(e^{j\omega})|');

    % Acercamiento alrededor de la frecuencia de corte
    subplot(2,3,3); plot(f, 20*log10(abs(H) + eps)); grid on; hold on;
    if caso == 1
        plot([60 60], [-100 5], 'r--'); xlim([0 fmaxZoom]);
    else
        xlim(fc + [-1500 1500]);
    end
    plot([fc fc], [-100 5], 'k:'); ylim([-100 5]); xlabel('f [Hz]'); ylabel('dB'); title('|H| ampliada');
    
    % Respuesta de fase
    % Se compara la fase real con la recta ideal
    % Característica de un filtro de fase lineal
    subplot(2,3,4);
    if caso == 1, k = f >= 1.5*fc & f <= 1000; else, k = f <= fc; end
    fase = unwrap(angle(H(k)));
    recta = -2*pi*f(k)/fs*M/2;
    recta = recta + (fase(1) - recta(1));
    plot(f(k), fase, f(k), recta, 'r--'); grid on;
    xlabel('f [Hz]'); ylabel('rad'); title('fase (banda de paso)'); legend('unwrap(\angle H)', '-\omega M/2 + cte');

    % Plano Z completo mostrando ceros y polos
    subplot(2,3,5); plot(cos(th), sin(th), 'k:'); hold on;
    plot(real(z), imag(z), 'bo', 'MarkerSize', 3);
    plot(0, 0, 'rx', 'MarkerSize', 10);
    text(0.05, 0.12, sprintf('%d polos', M));
    axis equal; grid on; xlim([-1.6 1.6]); ylim([-1.6 1.6]);
    xlabel('Re\{z\}'); ylabel('Im\{z\}'); title('ceros (o) y polos (x)');

    % Plano Z ampliado alrededor del ángulo de corte wc
    subplot(2,3,6); plot(cos(th), sin(th), 'k:'); hold on;
    plot(real(z), imag(z), 'bo', 'MarkerSize', 4);
    plot([0 2*cos(wc)], [0 2*sin(wc)], 'r--');      % rayo en el ángulo wc
    if caso == 1
        xlim([0.97 1.03]); ylim([-0.035 0.035]);       % A: cerca de z = 1 (f = 0)
    else
        xlim(cos(wc) + [-0.4 0.4]); ylim(sin(wc) + [-0.4 0.4]);
    end
    axis equal; grid on; xlabel('Re\{z\}'); ylabel('Im\{z\}'); title('ceros cerca de \omega_c');
end

%% Filtrado y Wav

% Se calcula el retardo introducido por el filtro FIR.
% Debido a la simetría de sus coeficientes, el retardo es constante e igual a M/2.
DA = (NA-1)/2;

% Se aplica el filtro mediante la función filter().
% Se añaden ceros al final para compensar el retardo y después
% se eliminan las primeras muestras retrasadas para realinear la señal.
yA = filter(hA, 1, [xA; zeros(DA,1)]);   yA = yA(DA+1 : end);
DB = (NB-1)/2;
yB = filter(hB, 1, [xB; zeros(DB,1)]);   yB = yB(DB+1 : end);

% Verificación de que la señal resultante no exceda ±1,
% evitando saturación al guardar el archivo de audio.
fprintf('Pico de salida: A = %.3f, B = %.3f (debe ser < 1 para no saturar el WAV)\n', ...
        max(abs(yA)), max(abs(yB)));

% Almacenamiento de las señales filtradas en formato WAV.
audiowrite('A_FIR.wav', yA, fs);
audiowrite('B_FIR.wav', yB, fs);

% Para escuchar antes/después:
 sound(xA, fs); pause(4.5); sound(yA, fs); pause(4.5);
 sound(xB, fs); pause(4.5); sound(yB, fs);

%% Comparaciòn
% Se calcula el espectro promedio utilizando ventanas Hann
% y solapamiento del 50 %, tal como se realizó en el análisis inicial.
Na    = 1024;  salto = Na/2;
wa    = 0.5*(1 - cos(2*pi*(0:Na-1)'/(Na-1)));  % Ventana Hann empleada para el análisis espectral.
fa    = (0:Na/2)' * fs/Na; % Vector de frecuencias del espectro.
X4    = [xA yA xB yB]; % Matriz que contiene las señales originales y filtradas.
P4    = zeros(Na/2+1, 4);  % promedio de todos los segmentos
S4    = zeros(Na/2+1, 4);  % promedio solo de segmentos sin voz
nSeg  = floor((length(xA) - Na)/salto) + 1; % Número total de segmentos analizados.
nSin  = 0;

% Procesamiento de cada segmento temporal.
for m = 1:nSeg
    i0 = (m-1)*salto + 1;
    sinVoz = (i0-1)/fs >= t1;  % el segmento empieza en el tramo sin voz
    nSin = nSin + sinVoz;
    for c = 1:4
        Xs = fft(X4(i0:i0+Na-1, c) .* wa);
        Pm = (2*abs(Xs(1:Na/2+1))/sum(wa)).^2;
        P4(:, c) = P4(:, c) + Pm;
        if sinVoz, S4(:, c) = S4(:, c) + Pm; end
    end
end
P4 = sqrt(P4 / nSeg);  % promedio (RMS) de los segmentos
S4 = sqrt(S4 / nSin);

% Graficas de comparacion
% Comparación entre las señales originales y filtradas.
figure(3);
subplot(2,2,1); plot(fa, 20*log10(P4(:,1)), fa, 20*log10(P4(:,2))); grid on;
xlim([0 500]); ylim([-110 -20]); title('A: ampliación 0-500 Hz'); xlabel('f [Hz]'); ylabel('dB');
legend('original', 'FIR');
subplot(2,2,2); plot(fa, 20*log10(P4(:,1)), fa, 20*log10(P4(:,2))); grid on;
xlim([0 fs/2]); ylim([-110 -20]); title('A: completo'); xlabel('f [Hz]');
subplot(2,2,3); plot(fa, 20*log10(P4(:,3)), fa, 20*log10(P4(:,4))); grid on;
xlim([0 fs/2]); ylim([-110 -20]); title('B: completo'); xlabel('f [Hz]'); ylabel('dB');
legend('original', 'FIR');
subplot(2,2,4); plot(fa, 20*log10(S4(:,3)), fa, 20*log10(S4(:,4))); grid on;
xlim([0 fs/2]); ylim([-110 -20]); title('B: solo ruido (tramo sin voz)'); xlabel('f [Hz]');

t  = (0:length(xB)-1)'/fs;
tr = (t >= t1) & (t <= t2); % tramo sin voz
rmsdB = @(v) 20*log10(sqrt(mean(v.^2))); % Función auxiliar para calcular el valor RMS en dB.
fprintf('\nRMS DEL TRAMO SIN VOZ (%.2f-%.2f s) \n', t1, t2);
fprintf('A: antes %6.1f dB  después %6.1f dB  -> interferencia atenuada %5.1f dB\n', ...
        rmsdB(xA(tr)), rmsdB(yA(tr)), rmsdB(xA(tr)) - rmsdB(yA(tr)));
fprintf('B: antes %6.1f dB  después %6.1f dB  -> interferencia atenuada %5.1f dB\n', ...
        rmsdB(xB(tr)), rmsdB(yB(tr)), rmsdB(xB(tr)) - rmsdB(yB(tr)));
fprintf('RMS DE TODO EL AUDIO (voz + interferencia) \n');
fprintf('A: antes %6.1f dB  después %6.1f dB\n', rmsdB(xA), rmsdB(yA));
fprintf('B: antes %6.1f dB  después %6.1f dB\n', rmsdB(xB), rmsdB(yB));

%% Ajuste de la frecuencia de corte

fprintf('\n PRUEBA DE fc PARA B (N = %d) \n', NB);
fprintf('  fc [Hz]   ruido atenuado [dB]   audio completo atenuado [dB]\n');
M = NB - 1;  n = (0:M)';
w = 0.5*(1 - cos(2*pi*(n+1)/(NB+1)));
for fcP = [2000 2500 3000 3500 4000]
    wc = 2*pi*fcP/fs;
    hP = zeros(NB, 1);
    for i = 1:NB
        d = n(i) - M/2;
        if d == 0, hP(i) = wc/pi; else, hP(i) = sin(wc*d)/(pi*d); end
    end
    hP = hP .* w;
    yP = filter(hP, 1, [xB; zeros(M/2,1)]);   yP = yP(M/2+1 : end);
    fprintf('  %6d        %6.1f                 %6.1f\n', fcP, ...
            rmsdB(xB(tr)) - rmsdB(yP(tr)), rmsdB(xB) - rmsdB(yP));
end
