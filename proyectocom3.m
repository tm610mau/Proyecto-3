% Caracterización del mensaje que se desea transmitir
f0 = 300; % frecuencia inicial
f1 = 1000; % frecuencia final

fs=2000; % tasa de muestreo

%t=0:1/fs:tf; % vector de tiempo
t = (0:1/fs:tf)';
tf = t(end); % tiempo final

signal = 3*chirp(t,f0,tf,f1); % up chirp

%sound(y, Fs)
figure(1)
plot(t, signal)
title('Señal up chirp')

l = length(t);
Y = fft(signal)/l;
Y2 = abs(fftshift(Y)); % correr la frecuencia cero al centro y aplicarle valor absoluto
fr = fs/2 * linspace(-1,1-2/l,l); %vector de frecuencias

%fr = Fs*(0:L-1);
figure(2)
plot(fr, Y2);
title('fft de la señal up chirp')

%-----------------------------------------------------------------------------------
%Modulación de la señal.

fdev1 = 100; %desviación de frecuencia de 100 Hz
fdev2 = 500; %desviación de frecuencia de 500 Hz
fc = 40000;
carrier = 3*cos(2*pi*fc*t); % portadora 
phi = 2*pi*fdev1*cumsum(signal)/fs; %phi(t)

fm_signal = cos(2*pi*fc*t + phi); % señal modulada (respuesta directa)  
%fm_signal = 3*cos(2*pi*fc*t).*cos(phi)-3*sin(2*pi*fc*t).*sin(phi);

figure(3)
plot(t, fm_signal)
title('Señal up chirp modulada')

%Ycarrier = fft(carrier)/l;
%figure(3)
%plot(fr, Ycarrier); %grafico de la portadora (solo para chequear)

Yfm = fft(fm_signal)/l;
figure(4)
plot(fr, Yfm);
title('fft de la señal up chirp modulada')

%-----------------------------------------------------------------------------------
%Agregar Ruido.



%-----------------------------------------------------------------------------------
%Demodulación de la señal.

t2 = (0:1/fs:((size(fm_signal,1)-1)/fs))';
t2 = t2(:,ones(1,size(fm_signal,2)));
hilb = hilbert(fm_signal).*exp(-1i*2*pi*fc*t2); % aplicar hilbert
demod_signal = (1/(2*pi*fdev1))*[zeros(1,size(hilb,2)); diff(unwrap(angle(hilb)))*fs];

figure(5)
plot(t, demod_signal)
title('Señal demodulada up chirp')

Ydemod = fft(demod_signal)/l;
Y2demod = abs(fftshift(Ydemod)); % correr la frecuencia cero al centro y aplicarle valor absoluto

figure(6)
plot(fr, Y2demod);
title('fft de la señal up chirp demodulada')
