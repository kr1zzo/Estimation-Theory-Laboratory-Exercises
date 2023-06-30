% Teorija estimacije
% 5. laboratorijska vježba 2022./2023.

% Zadatak 2.

%% 
clear all;
close all;
inicijalizacija;
model_zad2;
%% 
Nstm = 100; %broj simulacija

for i = 1:Nstm
   
    warning off
    sim model_zad2;

    [hat_theta{i} MSE(i)] = LSmetoda(u.signals.values,y.signals.values,2);

    a1(i) = hat_theta{i}(1);
    a2(i) = hat_theta{i}(2);
    b1(i) = hat_theta{i}(3);
    b2(i) = hat_theta{i}(4);
end

theta_ = ([mean(a1);mean(a2);mean(b1);mean(b2)])

bias = norm(theta - theta_)

%% Plot

figure(3)

plot(theta(1), theta(2), 'rO', 'LineWidth', 2)
hold on;
plot(a1, a2, 'bO', 'LineWidth', 2)
hold on;
%a1_ i a2_ iz inicijalizacija.m
plot(theta_(1), theta_(2), 'gO', 'LineWidth', 2)
hold on;
title("Parametri a1 i a2")
xlabel('a1')
ylabel('a2')
legend(strcat("Stvarni parametri: (",num2str(theta(1)) ," , ", num2str(theta(2)), ")"), "Estimirani parametri diskretiziranog sustava", strcat("Očekivani parametri(",num2str(theta_(1)) ," , ", num2str(theta_(2)), ")"))
grid on;
hold off;



