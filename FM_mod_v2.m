% Caracterización del mensaje que se desea transmitir

f0 = 300; % frecuencia inicial
f1 = 1000; % frecuencia final
fs1 = 3400;
fs2 = 1e+05;
%fs= 80000; % tasa de muestreo de 80 kHz
A = 3/2; % amplitud peak to peak

tf = 0.2; % tiempo final
%t=0:1/fs:tf; % vector de tiempo
t = (0:1/fs2:tf)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Digital LowPass Filter %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Frecuencia de corte de 1500 Hz

rf = 1000; % Cantidad de retardos 
nf = rf/2;
Tf = .01; % Intervalo de muestreo, en segundos
t_filter = 0:(1/fs2):Tf;

mf = (-rf/2):1:rf/2;

c_m = sin(13*pi*mf/500)./(pi*mf); % Formula de coeficientes (sinc)
c_m((rf/2)+1) = 13/500; % c_0, se indefine con la formula

w = hamming(rf+1);

% figure(100)
% plot(c_m)
% hold on
% plot(w)

cH_m = zeros(1, rf+1);

for i = 1:(rf+1)
    cH_m(i) = c_m(i)*w(i);
end

% figure(100) test plot coeficientes filtro
% plot(cH_m)
% hold on
% plot(c_m)

% lf = length(c_m);
% LP = fft(cH_m)/lf;
% LP2 = abs(fftshift(LP)); % correr la frecuencia cero al centro y aplicarle valor absoluto
% frLP = fs2/2 * linspace(-1,1-2/lf,lf); %vector de frecuencias
% 
% figure(200) % Test plot fourier del filtro
% plot(frLP,LP2, '.')


% Test de accion del filtro
% f10 = 900;
% f11 = 1700;
% 
% carrier_900 = cos(2*pi*f10*t_filter);
% carrier_2000 = cos(2*pi*f11*t_filter);
% 
% c900_filter = conv(carrier_900,cH_m,'same');
% c2000_filter = conv(carrier_2000,cH_m,'same');
% 
% figure(300) % Test plot accion del filtro
% plot(t_filter,c900_filter)
% hold on
% plot(t_filter,c2000_filter)
% legend('sinuoside de 900 Hz', 'sinusoide de 3000 Hz')
% %xlim([0 .1])
% title('accion del filtro que corta en 1500 Hz')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Digital LowPass Filter %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
ylabel('Amplitud')
xlim([-fs1/2 fs1/2])

%Modulación de la señal.

fdev1 = 100; %desviación de frecuencia de 100 Hz
fdev2 = 500; %desviación de frecuencia de 500 Hz
fc = 40000; %frecuencia de carrier
fd = 3e+03;

carrier = A*cos(2*pi*fc*t); % portadora 
carrier_90 = -A*sin(2*pi*fc*t); % portadora adelantada en 90°
phi1 = 2*pi*fdev1*cumsum(signal)/fs2; %phi1(t)
phi2 = 2*pi*fdev2*cumsum(signal)/fs2; %phi2(t)

fm_signal1 = A*cos(2*pi*fc*t + phi1); % señal modulada (respuesta directa)  
%fm_signal = 3*cos(2*pi*fc*t).*cos(phi)-3*sin(2*pi*fc*t).*sin(phi);
fm_signal2 = A*cos(2*pi*fc*t + phi2); % señal modulada (respuesta directa)

figure(3)
subplot(2,1,1);
plot(t, fm_signal1)
title('Señal up chirp modulada, \Deltaf = 100 Hz')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')

subplot(2,1,2);
plot(t, fm_signal2)
title('Señal up chirp modulada, \Deltaf = 500 Hz')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')

% transformada de fourier
Yfm1 = fft(fm_signal1)/l;
Yfm11 =abs(fftshift(Yfm1));
figure(4)
subplot(2,1,1);
plot(fr, Yfm11);
title('fft de la señal up chirp modulada, \Deltaf = 100 Hz')
xlabel('Frecuencia (Hz)', 'FontSize', 12, 'FontWeight', 'Bold')
xlim([fc - fd, fc + fd])
ylabel('Amplitud')

Yfm2 = fft(fm_signal2)/l;
Yfm22 =abs(fftshift(Yfm2));

subplot(2,1,2);
plot(fr, Yfm22);
title('fft de la señal up chirp modulada, \Deltaf = 500 Hz')
xlabel('Frecuencia (Hz)', 'FontSize', 12, 'FontWeight', 'Bold')
xlim([fc - fd, fc + fd])
ylabel('Amplitud')


% Demo ejemplo para agregar ruido

test_sine = A*cos(2*pi*100*t); % sinusoide 100 Hz
y_noise_15 = awgn(test_sine,15); % ruido de SNR = 15 dB a test_sine
y_noise_30 = awgn(test_sine,30); % ruido de SNR = 30 dB a test_sine

% figure(90) % comparacion ruidos
% plot(t,y_noise_15)
% hold on
% plot(t,y_noise_30, 'LineWidth', 2)

s1_15dB = awgn(fm_signal1,15);
s1_30dB = awgn(fm_signal1,30);
s2_15dB = awgn(fm_signal2,15);
s2_30dB = awgn(fm_signal2,30);

figure(5)
subplot(2,1,1);
plot(t, s1_15dB)
title('Señal up chirp modulada, \Deltaf = 100 Hz, SNR = 15 dB')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')

subplot(2,1,2);
plot(t, s1_30dB)
title('Señal up chirp modulada, \Deltaf = 100 Hz, SNR = 30 dB')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')

figure(6)
subplot(2,1,1);
plot(t, s2_15dB)
title('Señal up chirp modulada, \Deltaf = 500 Hz, SNR = 15 dB')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')

subplot(2,1,2);
plot(t, s2_30dB)
title('Señal up chirp modulada, \Deltaf = 500 Hz, SNR = 30 dB')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')

%Demodulación de la señal.

% s1_15dB_hilb = hilbert(s1_15dB);
% s1_30dB_hilb = hilbert(s1_30dB);
% s2_15dB_hilb = hilbert(s2_15dB);
% s2_30dB_hilb = hilbert(s2_30dB);


%%%%%%%%%%%%%%%%%%% Delta_f = 100 Hz, SNR = 15 dB %%%%%%%%%%%%%%%%%%%%%%%%%

s1_15dB_cos = carrier.*s1_15dB;
s1_15dB_sin = carrier_90.*s1_15dB;

I1_15 = conv(s1_15dB_cos,cH_m,'same');
Q1_15 = conv(s1_15dB_sin,cH_m,'same');

R1_15 = I1_15.*I1_15 + Q1_15.*Q1_15;

phi1_15 = fs2*(diff(Q1_15).*I1_15(1:end-1) - diff(I1_15).*Q1_15(1:end-1))...
    ./(2*pi*fdev1*R1_15(1:end-1));

L1_15 = length(phi1_15); % Transformada de Fourier
Y1_15 = fft(phi1_15)/L1_15;
Y1_15 = abs(fftshift(Y1_15)); % correr la frecuencia cero al centro y aplicarle valor absoluto
fr1_15 = (fs2/2)*linspace(-1,1-2/L1_15,L1_15); %vector de frecuencias

figure(7)
plot(fr1_15,Y1_15)
hold on
plot(fr,Y2)
xlim([-fs1/2 fs1/2])
xlabel('Frecuencia (Hz)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Amplitud', 'FontSize', 12, 'FontWeight', 'Bold')
title('Demodulación en espacio de Fourier, \Deltaf = 100 Hz, SNR = 15 dB ')
legend('Demodulación', 'Original')

figure(8)
plot(t(1:end-1),phi1_15)
hold on
plot(t, signal)
title('Señal up chirp, \Deltaf = 100 Hz, SNR = 15 dB')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')
legend('Demodulación', 'Original')

%%%%%%%%%%%%%%%%%%% Delta_f = 100 Hz, SNR = 30 dB %%%%%%%%%%%%%%%%%%%%%%%%%

s1_30dB_cos = carrier.*s1_30dB;
s1_30dB_sin = carrier_90.*s1_30dB;

I1_30 = conv(s1_30dB_cos,cH_m,'same');
Q1_30 = conv(s1_30dB_sin,cH_m,'same');

R1_30 = I1_30.*I1_30 + Q1_30.*Q1_30;

phi1_30 = fs2*(diff(Q1_30).*I1_30(1:end-1) - diff(I1_30).*Q1_30(1:end-1))...
    ./(2*pi*fdev1*R1_30(1:end-1));

L1_30 = length(phi1_30); % Transformada de Fourier
Y1_30 = fft(phi1_30)/L1_30;
Y1_30 = abs(fftshift(Y1_30)); % correr la frecuencia cero al centro y aplicarle valor absoluto
fr1_30 = (fs2/2)*linspace(-1,1-2/L1_30,L1_30); %vector de frecuencias

figure(9)
plot(fr1_30,Y1_30)
hold on
plot(fr,Y2)
xlim([-fs1/2 fs1/2])
xlabel('Frecuencia (Hz)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Amplitud', 'FontSize', 12, 'FontWeight', 'Bold')
title('Demodulación en espacio de Fourier, \Deltaf = 100 Hz, SNR = 30 dB ')
legend('Demodulación', 'Original')

figure(10)
plot(t(1:end-1),phi1_30)
hold on
plot(t, signal)
title('Señal up chirp, \Deltaf = 100 Hz, SNR = 30 dB')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')
legend('Demodulación', 'Original')

%%%%%%%%%%%%%%%%%%% Delta_f = 500 Hz, SNR = 15 dB %%%%%%%%%%%%%%%%%%%%%%%%%

s2_15dB_cos = carrier.*s2_15dB;
s2_15dB_sin = carrier_90.*s2_15dB;

I2_15 = conv(s2_15dB_cos,cH_m,'same');
Q2_15 = conv(s2_15dB_sin,cH_m,'same');

R2_15 = I2_15.*I2_15 + Q2_15.*Q2_15;

phi2_15 = fs2*(diff(Q2_15).*I2_15(1:end-1) - diff(I2_15).*Q2_15(1:end-1))...
    ./(2*pi*fdev2*R2_15(1:end-1));

L2_15 = length(phi2_15); % Transformada de Fourier
Y2_15 = fft(phi2_15)/L2_15;
Y2_15 = abs(fftshift(Y2_15)); % correr la frecuencia cero al centro y aplicarle valor absoluto
fr2_15 = (fs2/2)*linspace(-1,1-2/L2_15,L2_15); %vector de frecuencias

figure(11)
plot(t(1:end-1),phi2_15)
hold on
plot(t, signal)
title('Señal up chirp, \Deltaf = 500 Hz, SNR = 15 dB')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')
legend('Demodulación', 'Original')

figure(12)
plot(fr2_15,Y2_15)
hold on
plot(fr,Y2)
xlim([-fs1/2 fs1/2])
xlabel('Frecuencia (Hz)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Amplitud', 'FontSize', 12, 'FontWeight', 'Bold')
title('Demodulación en espacio de Fourier, \Deltaf = 500 Hz, SNR = 15 dB ')
legend('Demodulación', 'Original')

%%%%%%%%%%%%%%%%%%% Delta_f = 500 Hz, SNR = 30 dB %%%%%%%%%%%%%%%%%%%%%%%%%

s2_30dB_cos = carrier.*s2_30dB;
s2_30dB_sin = carrier_90.*s2_30dB;

I2_30 = conv(s2_30dB_cos,cH_m,'same');
Q2_30 = conv(s2_30dB_sin,cH_m,'same');

R2_30 = I2_30.*I2_30 + Q2_30.*Q2_30;

phi2_30 = fs2*(diff(Q2_30).*I2_30(1:end-1) - diff(I2_30).*Q2_30(1:end-1))...
    ./(2*pi*fdev2*R2_30(1:end-1));

L2_30 = length(phi2_30); % Transformada de Fourier
Y2_30 = fft(phi2_30)/L2_30;
Y2_30 = abs(fftshift(Y2_30)); % correr la frecuencia cero al centro y aplicarle valor absoluto
fr2_30 = (fs2/2)*linspace(-1,1-2/L2_30,L2_30); %vector de frecuencias

figure(13)
plot(fr2_30,Y2_30)
hold on
plot(fr,Y2)
xlim([-fs1/2 fs1/2])
xlabel('Frecuencia (Hz)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Amplitud', 'FontSize', 12, 'FontWeight', 'Bold')
title('Demodulación en espacio de Fourier, \Deltaf = 500 Hz, SNR = 30 dB ')
legend('Demodulación', 'Original')

figure(14)
plot(t(1:end-1),phi2_30)
hold on
plot(t, signal)
title('Señal up chirp, \Deltaf = 500 Hz, SNR = 30 dB')
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'Bold')
ylabel('Voltaje (V)', 'FontSize', 12, 'FontWeight', 'Bold')
legend('Demodulación', 'Original')

% % Graficos de prueba
% 
% figure(555)
% plot(phi2_30)
% 
% % transformada de Fourier
% L2_30 = length(phi2_30);
% Y2_30 = fft(phi2_30)/L2_30;
% Y2_30 = abs(fftshift(Y2_30)); % correr la frecuencia cero al centro y aplicarle valor absoluto
% fr2_30 = (fs2/2)*linspace(-1,1-2/L2_30,L2_30); %vector de frecuencias
% 
% figure(556)
% plot(fr2_30,Y2_30)
% xlim([-fs1/2 fs1/2])
