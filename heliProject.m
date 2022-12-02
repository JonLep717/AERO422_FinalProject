clear; clc; close all;

%% Helicopter Model
s = tf('s');

Om1 = 0.24; z1 = 0.2; % Phugoid
Om2 = 6.5; z2 = 0.05; % Short period
Om3 = 50; z3 = 0.001; % Flex mode

G1 = tf(Om1^2,[1 2*z1*Om1 Om1^2]);
G2 = tf(Om2^2,[1 2*z2*Om2 Om2^2]);
G3 = tf(Om3^2,[1 2*z3*Om3 Om3^2]);
G = 0.5*G1*G2*G3;

bodemag(G);
print -depsc plantBode.eps

%% Create disturbance and noise signals -- You can comment this section out and use the saved.mat file.
dt = 0.01; Fs = 1/dt;
F = 1/(s/5+1);
T = 0:dt:200;

% Create disturbance
D = 2*randn(length(T),1);
wc = (20/2/pi)/(Fs/2); % Normalized cutoff frequency (10 rqd/s)
[fb,fa] = butter(10,wc,'low'); % butter worth filter.
d = filter(fb,fa,D); % Now filter the disturbance. This generates colored noise
L = length(T)+1;

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Z = 2*fft(d,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1); % Only interested in signals upto Fs/2 Hz.
absZ = abs(Z(1:NFFT/2+1));
figure(3);
subplot(2,2,1);plot(T,d); 
xlabel('t'); ylabel('|d(t)|');

subplot(2,2,2);plot(2*pi*f,absZ); axis([0 100 0 .5]);
title('Spectrum of dist(t)');
xlabel('rad/s');
ylabel('|X(f)|');

% Create sensor noise
N = randn(length(T),1);
wc = (50/2/pi)/(Fs/2); % Normalized cutoff frequency (10 Hz)
[fb,fa] = butter(10,wc,'high'); % butter worth filter.
n = filter(fb,fa,N); % Now filter the disturbance. This generates colored noise
L = length(T)+1;

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Z = fft(n,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1); % Only interested in signals upto Fs/2 Hz.
absZ = 2*abs(Z(1:NFFT/2+1));
figure(3);
subplot(2,2,3);plot(T,n);
xlabel('t');
ylabel('n(t)');

subplot(2,2,4);plot(2*pi*f,absZ); axis([0 100 0 .5]);
title('Spectrum of noise(t)');
xlabel('rad/s');
ylabel('|X(f)|');

% Save Data
P = G;
distTime = d;
noiseTime = n;

save heli.mat P distTime noiseTime T
%% Try several control architectures

% First: Only proportional;
K = 10;

% Then try PID, lead/lag. Use rlocus(...) to determine the gain such that
% tr, Mp, ts are achieved.

C = K;
%% Analyze the controller's performance.

Gyr = C*G/(1+C*G);
Gyd = G/(1+C*G);
Gyn = -Gyr;
Gur = C/(1+C*G);


%% Time Simulation

[Y1,T1] = step(Gyr); 
Y2 = lsim(Gyd,d,T);
Y3 = lsim(Gyn,n,T);

u = step(Gur,T);

%%
figure(2); 
subplot(2,2,1); plot(T1,Y1); title('Step Response')
subplot(2,2,2); plot(T,u); title('Control')
subplot(2,2,3); bodemag(Gyr); title('Frequency Response Gyr')
subplot(2,2,4); bodemag(Gur); title('Frequency Response Gur')
