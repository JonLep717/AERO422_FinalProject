clear; clc; close all;

%% Helicopter Model

s = tf('s');

Om1 = 0.1; z1 = 0.2; % Phugoid
Om2 = 6.5; z2 = 0.7; % Short period
Om3 = 50; z3 = 0.01; % Flex mode

G1 = tf(Om1^2,[1 2*z1*Om1 Om1^2]);
G2 = tf(Om2^2,[1 2*z2*Om2 Om2^2]);
G3 = tf(Om3^2,[1 2*z3*Om3 Om3^2]);
G = 0.5*G1*G2*G3;

figure(1)
bodemag(G); title('Open Loop Frequency Response')
% print -depsc plantBode.eps

dt = 0.01; Fs = 1/dt;
F = 1/(s/5+1);
T = 0:dt:200;

%% Lead/Lag Controller

K = 30; %With only proportional control the system oscillates. Needs dampening

%Introduce Lead and Lag compensators

z1 = 0.1;
p1 = 11;
z2 = 0.3;
p2 = 0.25;

Lead = ((s/z1)+1)/((s/p1)+1);  %where p>>z>0

Lag = ((s+z2)/(s+p2)); %z>p>0 but z is close to p

%Having lots of issues with steady state error with just
%Lead/Lag

Ki = 5;

C = K*Lead*Lag + Ki/s;

Gyr = C*G/(1+C*G);

%% Time Simulation

Y1 = step(Gyr,T);

S = stepinfo(Y1)

figure(2); clf;
subplot(2,2,1); plot(T,Y1,'Linewidth',1); title('Gyr: Step Response'); xlabel('Time (s)'); grid on;

subplot(2,2,3); bodemag(Gyr); title('Gyr: Frequency Response'); grid on;