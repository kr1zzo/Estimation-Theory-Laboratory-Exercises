% Teorija estimacije
% 5. laboratorijska vježba 2022./2023.

% Zadatak 3.

%% 
clear all;
close all;
inicijalizacija;

%% a1

figure(1);

global phi1;
phi1 = 1;
sim("model_zad3")
time = data.time;

a1_real = realP.signals.values(:, 1);
a1_estim = rivP.signals.values(:, 1);
plot(time, a1_real, 'b');
hold on;
plot(time, a1_estim, 'r');
hold on;

%mijenjamo faktor zaboravljanja
global phi1;
phi1 = 0.975;
sim("model_zad3")

a1_est_faktor = rivP.signals.values(:, 1);
plot(time, a1_est_faktor, 'g');
xlabel('t[s]');
ylabel('a_1');
legend('Stvarni a1', 'Estimirani a1 bez faktora zaboravljanja', strcat('Estimirani a1 s faktorom zaboravljanja ',num2str(phi1)));
title("a1")
grid on;
hold off;

%%  a2

figure(2);

global phi1;
phi1 = 1;
sim("model_zad3")
time = data.time;

a1_real = realP.signals.values(:, 2);
a1_estim = rivP.signals.values(:, 2);
plot(time, a1_real, 'b');
hold on;
plot(time, a1_estim, 'r');
hold on;

%mijenjamo faktor zaboravljanja
global phi1;
phi1 = 0.975;
sim("model_zad3")

a1_est_faktor = rivP.signals.values(:,2);
plot(time, a1_est_faktor, 'g');
xlabel('t[s]');
ylabel('a_2');
legend('Stvarni a2', 'Estimirani a1 bez faktora zaboravljanja', strcat('Estimirani a1 s faktorom zaboravljanja ',num2str(phi1)));
title("a2")
grid on;
hold off;
