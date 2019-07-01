% Caracterización del mensaje que se desea transmitir
f0 = 300; % frecuencia inicial
f1 = 1000; % frecuencia final
fs1 = 3000;
fs2 = 1e+04;
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



%-----------------------------------------------------------------------------------
%Demodulación de la señal.

%hilb = hilbert(fm_signal).*exp(-1i*2*pi*fc*t); % aplicar hilbert
%demod_signal = diff(unwrap(angle(hilb))); %diferencial del angulo de hilbert
%demod_signal = demod_signal*fs2/(2*pi*fdev1); % cambiar amplitud

demod_signal1 = hilbert(fm_signal1).*exp(-1i*2*pi*fc*t); % aplicar hilbert
I1 = real(demod_signal1); % I
Q1 = imag(demod_signal1); % Q
demod_signal1 = diff(atan(Q1/I1))*fs2/(2*pi*fdev1) ; % señal demodulada

demod_signal2 = hilbert(fm_signal2).*exp(-1i*2*pi*fc*t); % aplicar hilbert
I2 = real(demod_signal2); % I
Q2 = imag(demod_signal2); % Q
demod_signal2 = diff(atan(Q2/I2))*fs2/(2*pi*fdev2) ;  % señal demodulada



figure(7)
plot(t(1:l-1), demod_signal1)
title('Señal demodulada up chirp \Deltaf = 100 Hz')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')
xlim([0,0.2])
figure(8)
plot(t(1:l-1), demod_signal2)
title('Señal demodulada up chirp \Deltaf = 500 Hz')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')
xlim([0,0.2])
% transformada de fourier
Ydemod1 = fft(demod_signal1)/l; % Fourier
Ydemod1 = abs(fftshift(Ydemod1)); % correr la frecuencia cero al centro y aplicarle valor absoluto

Ydemod2 = fft(demod_signal1)/l; % Fourier
Ydemod2 = abs(fftshift(Ydemod2)); % correr la frecuencia cero al centro y aplicarle valor absoluto

figure(9)
plot(fr(1:l-1), Ydemod1);
title('fft de la señal up chirp modulada, \Deltaf = 100 Hz')
xlabel('Frecuencia (Hz)', 'FontSize', 12, 'FontWeight', 'Bold')

figure(10)
plot(fr(1:l-1), Ydemod2);
title('fft de la señal up chirp demodulada, \Deltaf = 500 Hz')
xlabel('Frecuencia (Hz)', 'FontSize', 12, 'FontWeight', 'Bold')
