clear; close all; clc;
if exist('OCTAVE_VERSION', 'builtin'), pkg load signal; end   % solo en Octave

[xA, fs] = audioread('PDS_voz_60Hz_4s.wav');              % A: voz + 60 Hz
xB       = audioread('PDS_voz_fuente_simulada_4s.wav');   % B: voz + agua
fA       = audioread('A_FIR.wav');                        % salidas FIR (etapa 2,
fB       = audioread('B_FIR.wav');                        %  ya sin retardo M/2)

%% DISEÑO IIR
% Frecuencias normalizadas respecto a fs/2.
% Filtro A: Butterworth rechaza-banda de 55 a 65 Hz.
% Filtro B: Butterworth pasa-bajas con fc = 3500 Hz.

[bA, aA] = butter(2, [55 65]/(fs/2), 'stop');
[bB, aB] = butter(8, 3500/(fs/2));

%% POLOS, CEROS Y |H| 
% Polos y ceros obtenidos a partir de los coeficientes del filtro.
% El sistema es estable si todos los polos están dentro del círculo unitario.
Nfft = 64000;   % Resolución en frecuencia
f    = (0:Nfft/2)' * fs/Nfft;
HA   = fft(bA(:), Nfft) ./ fft(aA(:), Nfft);   HA = HA(1:Nfft/2+1);
HB   = fft(bB(:), Nfft) ./ fft(aB(:), Nfft);   HB = HB(1:Nfft/2+1);
dB   = @(H, fx) 20*log10(abs(H(fx*Nfft/fs + 1)));          % |H| en una f [dB]

% Cálculo del retardo de grupo.
gd = @(H, fx) -diff(unwrap(angle(H(fx*Nfft/fs + [1 2])))) / (2*pi/Nfft);

fprintf('A IIR: orden %d | max|polo| = %.4f | |H(60)| = %.1f dB | |H(120)| = %.2f dB\n', ...
        length(aA)-1, max(abs(roots(aA))), dB(HA,60), dB(HA,120));
fprintf('       retardo: %.1f muestras a 120 Hz, %.2f a 1 kHz   (FIR: 500 constante)\n', ...
        gd(HA,120), gd(HA,1000));
fprintf('B IIR: orden %d | max|polo| = %.4f | |H(3k)| = %.2f dB | |H(4k)| = %.1f dB\n', ...
        length(aB)-1, max(abs(roots(aB))), dB(HB,3000), dB(HB,4000));
fprintf('       retardo: %.1f muestras a 500 Hz, %.1f a 3 kHz   (FIR: 50 constante)\n\n', ...
        gd(HB,500), gd(HB,3000));

th = linspace(0, 2*pi, 400);
figure(1);
% Plano Z del filtro rechaza-banda.
subplot(2,2,1); plot(f, 20*log10(abs(HA))); grid on; xlim([0 300]); ylim([-80 5]);
xlabel('f [Hz]'); ylabel('dB'); title('A: |H| IIR rechaza-banda');
subplot(2,2,2); plot(f, 20*log10(abs(HB))); grid on; ylim([-80 5]);
xlabel('f [Hz]'); ylabel('dB'); title('B: |H| IIR pasa-bajas');
% Plano Z del filtro pasa-bajas.
subplot(2,2,3); plot(cos(th), sin(th), 'k:', real(roots(bA)), imag(roots(bA)), 'bo', ...
                     real(roots(aA)), imag(roots(aA)), 'rx'); grid on;
xlim([0.9975 1.0005]); ylim([0.021 0.026]); xlabel('Re\{z\}'); ylabel('Im\{z\}');
title('A: plano Z (zoom en 60 Hz)');
subplot(2,2,4); plot(cos(th), sin(th), 'k:', real(roots(bB)), imag(roots(bB)), 'bo', ...
                     real(roots(aB)), imag(roots(aB)), 'rx'); grid on; axis equal;
xlabel('Re\{z\}'); ylabel('Im\{z\}'); title('B: plano Z (8 ceros en z = -1)');

%% FILTRADO CAUSAL Y EXPORTACIÓN
% Aplicación de los filtros IIR y almacenamiento de los audios resultantes.

yA = filter(bA, aA, xA);
yB = filter(bB, aB, xB);
audiowrite('A_IIR.wav', yA, fs);
audiowrite('B_IIR.wav', yB, fs);

%% COMPARACIÓN ENTRE FILTROS FIR E IIR
% Comparación mediante atenuación de interferencia y análisis espectral.
t   = (0:length(xA)-1)'/fs;
tr  = t >= 3.80 & t <= 3.95;
% Cálculo del nivel RMS en un tramo sin voz.
rms = @(v) 20*log10(sqrt(mean(v(tr).^2)));
fprintf('Interferencia atenuada (tramo sin voz):\n');
fprintf('  A:  FIR %.1f dB   IIR %.1f dB\n', rms(xA)-rms(fA), rms(xA)-rms(yA));
fprintf('  B:  FIR %.1f dB   IIR %.1f dB\n', rms(xB)-rms(fB), rms(xB)-rms(yB));

% Análisis espectral entre 3 y 4 s.
Na = 4096;  wa = 0.5*(1 - cos(2*pi*(0:Na-1)'/(Na-1)));  fa = (0:Na/2)'*fs/Na;
X  = [xA fA yA xB fB yB];   P = zeros(Na/2+1, 6);   inicios = 48001:Na/2:64001-Na;
for i0 = inicios
    S = fft(X(i0:i0+Na-1, :) .* wa);                 % las 6 señales a la vez
    P = P + abs(S(1:Na/2+1, :)).^2;
end
P = 10*log10(P / length(inicios) * (2/sum(wa))^2);

figure(2);
% Comparación espectral del caso A.
subplot(1,2,1); plot(fa, P(:,1:3)); grid on; xlim([0 300]); ylim([-120 -20]);
xlabel('f [Hz]'); ylabel('dB'); title('A: espectro 3-4 s'); legend('original', 'FIR', 'IIR');
% Comparación espectral del caso B.
subplot(1,2,2); plot(fa, P(:,4:6)); grid on; xlim([0 fs/2]); ylim([-120 -20]);
xlabel('f [Hz]'); title('B: espectro 3-4 s'); legend('original', 'FIR', 'IIR');
