% Toerija estimacije
% 5. laboratorijska vježba 2022./2023.

% Zadatak 1

%% 
clear
% Parameters 1
m=3;
c=1;
k=100;

% Parameters 2
m2=2;
c2=1;
k2=100;

% Continuous models
Gc=tf([1],[m c  k]);
Gc2=tf([1],[m2 c2  k2]);

%% PRBS signal parameters -> a)
%%%%%%%%%%%%%%%%%%%%%%% FILL VALUES %%%%%%%%%%%%%%%%%%%%%%%
figure(1)
bodeplot(Gc)
grid on;
title('Prijenosna funkcija Gc')

figure(2)
step(Gc)
grid on;
%iz step je vidljivo da je stacionarno stanje 0.01 = -40 db
stacionarno_stanje = mag2db(0.01)
% Graničnu frekvenciju odredite kao frekvenciju na kojoj sustav ima 
% prigušenje 3dB u odnosu na pojačanje sustava u stacionarnom stanju
granicno_pojacanje = stacionarno_stanje - 3;

[mag,phase,wout] = bode(Gc);
mag = squeeze(mag);
phase = squeeze(phase);

w_g = interp1(20*log10(mag), wout, granicno_pojacanje)

dt = 2.77 / w_g;

G = stepinfo(Gc, 'SettlingTimeTreshold', 0.05);
settling_time = G.SettlingTime 

T = 1.5*settling_time;

N = T/dt;
n = ceil(log2(N+1));

N = (2^n)-1
dt = T/N
%%%%%%%%%%%%%%%%%%%%%%% FILL VALUES %%%%%%%%%%%%%%%%%%%%%%%
%% Diskretizacija -> b)

Ts=dt;
Tsim = 10*Ts*N;
% PRBS signal
ui=idinput(N,'prbs');

% Discrete models
Gd=c2d(Gc,Ts)
Gd2=c2d(Gc2,Ts)

%
theta=[Gd.den{1}(2:3) Gd.num{1}(2:3)]'
theta2=[Gd2.den{1}(2:3) Gd2.num{1}(2:3)]';

a1_ = theta(1);
a2_ = theta(2);
b1_ = theta(3);
b2_ = theta(4);

