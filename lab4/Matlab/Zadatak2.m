% Teorija estimacije
% 4. laboratorijska vježba 2022./2023.

% Zadatak 2

%% 
close all
clear all;
tmass_ss;
Am = 1;
syms x

%% a)
%% 
[num,den] = ss2tf(A,B,C(1,:),D(1),1);

G_ideal = tf(num, den)
G = stepinfo(G_ideal, 'SettlingTimeTreshold', 0.05);
%vrijeme smirivanja ozdia t_95
settling_time = G.SettlingTime 
%korelacijsko vrijeme
T = 1.5*settling_time;

%zadano
Suu = 0.00691;

%slide 39 i 44
wg = 60; %zadano u Pripremi za vježbu
w_max = wg;
dt = 2.77 / w_max
Ts = dt;

% N?
N = T/dt;
n = ceil(log2(N+1));
N = (2^n)-1

%c?
w0 = (2*pi)/(N*dt);
ni = round(wg/w0);

c = solve(Suu == x^2*((N+1)/N^2)*(sin(ni*pi/N)/(ni*pi/N))^2, x);
c = double(abs(c(1)))


%za simulaciju
signal = idinput(N,'PRBS', [0, 1],[-c,c]);
Tsim = T*10;

sim('tromaseni_model', Tsim);

%% b) 
%% 
figure(1);
hold on;
stairs(data.signals.values(:,1))
title('Ulazni signal')
xlabel('\itN')
ylabel('\itu')
xlim([0 255])
grid on
hold off;

figure(2);
hold on
plot(data.signals.values(:,2))
title('Izlazni signal');
xlabel('\itN')
ylabel('\ity')
xlim([0 255])
grid on
hold off;

% zadana cra funkcija 
cra =cra([data.signals.values(:,2),data.signals.values(:,1)],255);

figure(3)
hold on;
plot(linspace(0, dt*255,255), cra(1:255)*dt^-1);
impulse(G_ideal,'r');
title('Težinska funkcija idealnog i estimiranog modela')
legend('Estimirani odziv','Idealni odziv')
grid on
hold off;