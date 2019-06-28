% Caracterización del mensaje que se desea transmitir
f0 = 300; % frecuencia inicial
f1 = 1000; % frecuencia final

Fs=2100; % tasa de muestreo
tf=0.4; % tiempo final
t=0:1/Fs:tf; % vector de tiempo
t1 = t(end);

signal = 3*chirp(t,f0,tf,f1); % up chirp

%sound(y, Fs)
figure(1)
plot(t, signal)
title('Señal up chirp')



l = length(t);
Y = fft(signal)/l;
Y2 = abs(fftshift(Y)); % correr la frecuencia cero al centro y aplicarle valor absoluto
fr = Fs/2 * linspace(-1,1-2/l,l); %vector de frecuencias

%fr = Fs*(0:L-1);
figure(2)
plot(fr, Y2);
title('fft de la señal up chirp')

%-----------------------------------------------------------------------------------
%Modulación de la señal.

