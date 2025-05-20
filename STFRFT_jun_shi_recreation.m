clc; clear all;

%modifiable params
fs=1000;        %sampling freq
k=3;            %chirp rate
alpha=-acot(k); %rotation angle for frft

%fixed params
t= 0:1/fs:12;
t1= 0:1/fs:4-1/fs;
t2= 4:1/fs:8-1/fs;
t3= 8:1/fs:12;

%input signal
signal=complex(zeros(1,length(t)));
signal(1:length(t1)) = 0.9*exp(50*t1*1j+1j*k*(t1.^2)/2);
signal(length(t1)+1:length(t2)+length(t1))= 0.7*exp(10*t2*1j+1j*k*(t2.^2)/2);
signal(length(t1)+length(t2)+1:length(t1)+length(t2)+length(t3)) = 0.8*exp(30*t3*1j+1j*k*(t3.^2)/2);

%input signal plot
subplot(2,2,1);
plot(t,real(signal));
xlim([0,12]);
ylim([-1,1]);
xticks(0:2:12);
yticks(-1:0.2:1);
grid on;
xlabel("time");
ylabel("amplitude");
title("Time Domain Waveform s(t)");

%FRFT Kernel default (verification needed)
u = linspace(-30, 30, length(t));  % same number of points as t

[T, U] = meshgrid(t, u);
A = sqrt((1 - 1j*cot(alpha)) / (2*pi));
K = A * exp(1j * (U.^2 + T.^2)/2 * cot(alpha) - 1j * U .* T * csc(alpha));

F_alpha = K * signal.' * (1/fs);
subplot(2,2,2);
plot(u, abs(F_alpha));
xlabel('u'); ylabel('|F_\alpha(u)|');
title('Matched-Angle FRFT');