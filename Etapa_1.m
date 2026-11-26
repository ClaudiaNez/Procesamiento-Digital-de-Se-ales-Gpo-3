clear; close all; clc;

%Archivos de entrada
archA = 'PDS_voz_60Hz_4s.wav';              
archB = 'PDS_voz_fuente_simulada_4s.wav';    

N_principal = 1024;  % --- Longitud de segmento N                 
N_lista = [256 512 1024 2048]; %Lista de longitu  des que se comparan entre sí
f_zoom = 500; % Límite superior del zoom de bajas frecuencias (Hz)
sin_voz   = 3.80;  % [s] desde aquí ya no hay voz 

%%  Carga y tiempo
% Se cargan ambos archivos de audio mediante audioread, obteniendo la señal digital normalizada y su frecuencia de muestreo.
[xA, fs]  = audioread(archA);
[xB, fsB] = audioread(archB);

% assert(condición, mensaje): si la condición es falsa, se detiene el script
assert(size(xA,2) == 1 && size(xB,2) == 1, 'Los audios deben ser mono');
assert(fs == fsB, 'Los audios deben tener la misma fs');

% Reproduciomos los audios
 sound(xA, fs);  pause(4.5);  sound(xB, fs);

% Generamos un vector de tiempo a partir de la frecuencia de muestreo, para poder representar las señales en función del tiempo.
t = (0:length(xA)-1)' / fs;          


% Se selecciona una ventana de análisis entre 3.85 s y 3.95 s, ya que dentro de esa region, no hay voz.
% Considerando que la interferencia esperada es de 60 Hz, un intervalo de
% 100 ms contiene aproximadamente seis periodos completos de dicha señal.
zoom = (t >= 3.85) & (t <= 3.95);

%Graficamos la señal completa
figure(1);
subplot(2,2,1); plot(t, xA); grid on; title('A completo'); xlabel('t [s]');
subplot(2,2,3); plot(t, xB); grid on; title('B completo'); xlabel('t [s]');
%Graficmos la ampliacion temporal
subplot(2,2,2); plot(t(zoom)*1000, xA(zoom)); grid on; title('A sin voz (zoom)'); xlabel('t [ms]');
subplot(2,2,4); plot(t(zoom)*1000, xB(zoom)); grid on; title('B sin voz (zoom)'); xlabel('t [ms]');

%% Sesgmentos, HANN, FFT

N    = N_principal;
salto = N/2;                                 % 50 % de solape
w    = 0.5*(1 - cos(2*pi*(0:N-1)'/(N-1)));   % ventana Hann (vector columna)
f    = (0:N/2)' * fs/N;                      % frecuencias de 0 a fs/2
nSeg = floor((length(xA) - N)/salto) + 1;    % cuántos segmentos caben

SA = zeros(N/2+1, nSeg);    % espectros de A
SB = zeros(N/2+1, nSeg);    % espectros de B
tSeg = zeros(1, nSeg);      % instante central de cada segmento

for m = 1:nSeg
    i0 = (m-1)*salto + 1;                    % primera muestra del segmento
    i1 = i0 + N - 1;                         % última muestra del segmento

    XA = fft(xA(i0:i1) .* w);                % segmento x ventana -> FFT
    XB = fft(xB(i0:i1) .* w);

    SA(:, m) = 2*abs(XA(1:N/2+1)) / sum(w);  % solo tomamos la parte positiva
    SB(:, m) = 2*abs(XB(1:N/2+1)) / sum(w);
    tSeg(m) = (i0 - 1 + N/2) / fs;
end

% Con los espectros obtenidos para cada segmento se construyen los espectrogramas de las grabaciones A y B. 
% %Para ello se utiliza la función imagesc, que representa mediante colores la amplitud espectral en función del tiempo y la frecuencia.
figure(2);
subplot(2,1,1);
imagesc(tSeg, f, 20*log10(SA + eps)); axis xy; colorbar; caxis([-110 -20]);
title(sprintf('A: espectrograma, Hann, N = %d', N)); ylabel('f [Hz]');
subplot(2,1,2);
imagesc(tSeg, f, 20*log10(SB + eps)); axis xy; colorbar; caxis([-110 -20]);% (eps es un número muy pequeño; se suma para no calcular log10(0).)
title(sprintf('B: espectrograma, Hann, N = %d', N)); ylabel('f [Hz]'); xlabel('t [s]');

%% Espectro pormedio para cada N

figure(3);
fprintf('    N   duración [ms]   Δf [Hz]   pico de A entre 20 y 200 Hz\n');

for q = 1:length(N_lista)
    Nq     = N_lista(q);
    saltoq = Nq/2;
    wq     = 0.5*(1 - cos(2*pi*(0:Nq-1)'/(Nq-1)));
    fq     = (0:Nq/2)' * fs/Nq;
    nSegq  = floor((length(xA) - Nq)/saltoq) + 1;

    PA = zeros(Nq/2+1, 1);   % suma de potencias de A
    PB = zeros(Nq/2+1, 1);   % suma de potencias de B
    for m = 1:nSegq
        i0 = (m-1)*saltoq + 1;
        i1 = i0 + Nq - 1;
        XA = fft(xA(i0:i1) .* wq);
        XB = fft(xB(i0:i1) .* wq);
        PA = PA + (2*abs(XA(1:Nq/2+1)) / sum(wq)).^2;
        PB = PB + (2*abs(XB(1:Nq/2+1)) / sum(wq)).^2;
    end
    promA = sqrt(PA / nSegq);    % espectro promedio de A
    promB = sqrt(PB / nSegq);    % espectro promedio de B

    subplot(2,2,1); hold on; plot(fq, 20*log10(promA));   % A completo
    subplot(2,2,2); hold on; plot(fq, 20*log10(promB));   % B completo
    subplot(2,2,3); hold on; plot(fq, 20*log10(promA));   % A zoom
    subplot(2,2,4); hold on; plot(fq, 20*log10(promB));   % B zoom

    % ¿En qué frecuencia aparece el tono? Solo puede caer en un múltiplo de Δf.
    banda = (fq >= 20) & (fq <= 200);
    fb = fq(banda);   pb = promA(banda);
    [~, k] = max(pb);
    fprintf('%5d   %10.1f   %8.2f   %8.2f Hz\n', Nq, Nq/fs*1000, fs/Nq, fb(k));
end

subplot(2,2,1); grid on; xlim([0 fs/2]);   ylim([-120 0]); title('A: promedio vs N');
subplot(2,2,2); grid on; xlim([0 fs/2]);   ylim([-120 0]); title('B: promedio vs N');
legend('N = 256', 'N = 512', 'N = 1024', 'N = 2048', 'Location', 'southwest');
subplot(2,2,3); grid on; xlim([0 f_zoom]); ylim([-100 0]); title('A: zoom'); xlabel('f [Hz]');
subplot(2,2,4); grid on; xlim([0 f_zoom]); ylim([-100 0]); title('B: zoom'); xlabel('f [Hz]');

%% Ventana Rectangular vs Hann

Nr   = 256;
Nfft = 16*Nr;
seg  = xA(end-Nr+1 : end);                       % últimas 256 muestras
wh   = 0.5*(1 - cos(2*pi*(0:Nr-1)'/(Nr-1)));     % Hann
fz   = (0:Nfft/2)' * fs/Nfft;

Xr = fft(seg, Nfft);           Ar = 2*abs(Xr(1:Nfft/2+1)) / Nr;        % rectangular
Xh = fft(seg .* wh, Nfft);     Ah = 2*abs(Xh(1:Nfft/2+1)) / sum(wh);   % Hann

figure(4);
plot(fz, 20*log10(Ar), 'b', fz, 20*log10(Ah), 'r'); grid on;
xlim([0 1000]); ylim([-110 0]); xlabel('f [Hz]'); ylabel('dB');
title('A sin voz, N = 256: rectangular vs Hann');
legend('rectangular', 'Hann');

%% Comparación entre el audio completo y el tramo sin voz
%  Los segmentos después de sinVoz contienen solo la interferencia

sinVoz = (tSeg - (N/2)/fs) >= sin_voz;    % segmentos que empiezan después de t_sin_voz

promA_todo   = sqrt(mean(SA.^2, 2));            % promedio de todos los segmentos
promA_sinVoz = sqrt(mean(SA(:, sinVoz).^2, 2)); % promedio solo sin voz
promB_todo   = sqrt(mean(SB.^2, 2));
promB_sinVoz = sqrt(mean(SB(:, sinVoz).^2, 2));

figure(5);
subplot(2,1,1);
plot(f, 20*log10(promA_todo + eps), f, 20*log10(promA_sinVoz + eps)); grid on;
xlim([0 fs/2]); ylim([-120 0]); ylabel('dB'); title('A: audio completo vs sin voz');
legend('todo el audio', 'solo interferencia');
subplot(2,1,2);
plot(f, 20*log10(promB_todo + eps), f, 20*log10(promB_sinVoz + eps)); grid on;
xlim([0 fs/2]); ylim([-120 0]); ylabel('dB'); xlabel('f [Hz]'); title('B: audio completo vs sin voz');
legend('todo el audio', 'solo interferencia');
