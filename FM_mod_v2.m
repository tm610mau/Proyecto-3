% Caracterización del mensaje que se desea transmitir
f0 = 300; % frecuencia inicial
f1 = 1000; % frecuencia final
fs1 = 3000;
fs2 = 1e+05;
%fs= 80000; % tasa de muestreo de 80 kHz
A = 3/2; % amplitud peak to peak

tf = 0.2; % tiempo final
%t=0:1/fs:tf; % vector de tiempo
t = (0:1/fs2:tf)';

signal = A*chirp(t,f0,tf,f1); % up chirp

figure(1)
plot(t, signal)
title('Señal up chirp')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')

% transformada de Fourier
l = length(t);
Y = fft(signal)/l;
Y2 = abs(fftshift(Y)); % correr la frecuencia cero al centro y aplicarle valor absoluto
fr = fs2/2 * linspace(-1,1-2/l,l); %vector de frecuencias

%fr = Fs*(0:L-1);
figure(2)
plot(fr, Y2);
title('fft de la señal up chirp')
xlabel('Frecuencia (Hz)', 'FontSize', 12, 'FontWeight', 'Bold')
xlim([-fs1/2 fs1/2])

%Modulación de la señal.

fdev1 = 100; %desviación de frecuencia de 100 Hz
fdev2 = 500; %desviación de frecuencia de 500 Hz
fc = 40000; %frecuencia de carrier
fd = 3e+03;

carrier = A*cos(2*pi*fc*t); % portadora 
phi1 = 2*pi*fdev1*cumsum(signal)/fs2; %phi1(t)
phi2 = 2*pi*fdev2*cumsum(signal)/fs2; %phi2(t)

fm_signal1 = A*cos(2*pi*fc*t + phi1); % señal modulada (respuesta directa)  
%fm_signal = 3*cos(2*pi*fc*t).*cos(phi)-3*sin(2*pi*fc*t).*sin(phi);
fm_signal2 = A*cos(2*pi*fc*t + phi2); % señal modulada (respuesta directa)

figure(3)
plot(t, fm_signal1)
title('Señal up chirp modulada, \Deltaf = 100 Hz')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')

figure(5)
plot(t, fm_signal2)
title('Señal up chirp modulada, \Deltaf = 500 Hz')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')

% transformada de fourier
Yfm1 = fft(fm_signal1)/l;
Yfm11 =abs(fftshift(Yfm1));
figure(4)
plot(fr, Yfm11);
title('fft de la señal up chirp modulada, \Deltaf = 100 Hz')
xlabel('Frecuencia (Hz)', 'FontSize', 12, 'FontWeight', 'Bold')
xlim([fc - fd, fc + fd])

Yfm2 = fft(fm_signal2)/l;
Yfm22 =abs(fftshift(Yfm2));
figure(6)
plot(fr, Yfm22);
title('fft de la señal up chirp modulada, \Deltaf = 500 Hz')
xlabel('Frecuencia (Hz)', 'FontSize', 12, 'FontWeight', 'Bold')
xlim([fc - fd, fc + fd])


% Demo ejemplo para agregar ruido

test_sine = A*cos(2*pi*100*t); % sinusoide 100 Hz
y_noise_15 = awgn(test_sine,15); % ruido de SNR = 15 dB a test_sine
y_noise_30 = awgn(test_sine,60); % ruido de SNR = 30 dB a test_sine

figure(90) % comparacion ruidos
plot(t,y_noise_15)
hold on
plot(t,y_noise_30, 'LineWidth', 2)