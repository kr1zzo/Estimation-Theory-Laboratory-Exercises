% Teorija estimacije
% 4. laboratorijska vježba 2022./2023.

% Zadatak 1.

%% 
clear all;
close all;
signal = 1;
tmass_ss;
%% a)
%% 
Ts = 0.01;
Tsim = 2*pi/0.1;

Am = 1;
sim('tromaseni_model');

FFT_chirp = fft(data.signals.values(:,1));
FFT_sum = FFT_chirp(1:628);
FFT_avg = sum(abs(FFT_sum))/length(FFT_sum)

%srednja vrijednost amplitude izlaza dobivenog izvršavanjem fft naredbe u zanimljivom rasponu
%frekvencija bude približno 188.
Am = 188 / FFT_avg;

sim('tromaseni_model');

%% Plot

figure(1)
hold on
grid on
plot(data.time, data.signals.values(:,1),'r');
title('(i) Pobudni chirp signal')
ylabel('\itu(t)')
xlabel('\it{t}[s]')


figure(2)
hold on
grid on
plot(data.time, data.signals.values(:,2))
title('(ii) Tromaseni model(pomak prve mase x1) sa šumom')
ylabel('\itx_1(t)')
xlabel('\it{t}[s]')
hold on

figure(3)
hold on
grid on
plot(data.time, data.signals.values(:,3))
title('(iii) Idealni odziv')
ylabel('\itx_1(t)')
xlabel('\it{t}[s]')
hold on;

figure(4)
hold on
grid on
plot(data.time, data.signals.values(:,2))
plot(data.time, data.signals.values(:,3))
title('(iii) Idealni odziv vs (ii) Odziv modela')
ylabel('\itx_1(t)')
xlabel('\it{t}[s]')
legend('Odziv sa šumom', 'Idealni odziv')
hold off

%% Uklanjanje DC komponente
without_DC_1 = detrend(data.signals.values(:,1));
without_DC_2 = detrend(data.signals.values(:,2));
without_DC_3 = detrend(data.signals.values(:,3));

vector = linspace(0.1, Tsim, length(without_DC_1));
%% Plot

figure(5)
hold on
grid on
plot(data.time, without_DC_1,'r');
title('Pobudni chirp signal bez DC komponente')
ylabel('\itu(t)')
xlabel('\it{t}[s]')

figure(6)
hold on
grid on
plot(data.time, without_DC_2)
title('(ii) Tromaseni model(pomak prve mase x1) sa šumom bez DC komponente')
ylabel('\itx_1(t)')
xlabel('\it{t}[s]')
hold on

figure(7)
hold on
grid on
plot(data.time, without_DC_3)
title('(iii) Idealni odziv bez DC komponente')
ylabel('\itx_1(t)')
xlabel('\it{t}[s]')
hold on;

figure(8)
hold on
grid on
plot(data.time, without_DC_2)
plot(data.time, without_DC_3)
title('(iii) Idealni odziv vs (ii) Odziv modela (bez DC komponente)')
ylabel('\itx_1(t)')
xlabel('\it{t}[s]')
legend('Odziv sa šumom', 'Idealni odziv')
hold off

%% b)
%% Izracun spektra autokorelacijske i medukorelacijske funkcije

auto_corr = xcorr(data.signals.values(:,1));
inter_corr = xcorr(data.signals.values(:,1),data.signals.values(:,3));

%% Plot

figure(9)
hold on;
plot(auto_corr)
title('Autokorelacijska funkcija')
grid on


figure(10)
hold on;
plot(inter_corr)
title('Međukorelacijska funkcija')
grid on

%% iz dobivenih grafova:
H1 = 500;
H2 = 100;
H3 = 3000;

%% Bodeovi dijagrami

figure(11)
hold on
tf_1 = spektralnaAnaliza(without_DC_2, data.signals.values(:,1), H1, vector, Ts);
bodeplot(tf(num,den), 'r')
bodeplot(tf_1, 'b')
xlim([0.1 1000])
grid on
title(['H = ',num2str(H1)])
legend('Idealno', 'Estimirano')

figure(12)
hold on
tf_2 = spektralnaAnaliza(without_DC_2, data.signals.values(:,1), H2, vector, Ts);
bodeplot(tf(num,den), 'r')
bodeplot(tf_2, 'b')
grid on
title(['H = ',num2str(H2)])
legend('Idealno', 'Estimirano')

figure(13)
hold on
tf_3 = spektralnaAnaliza(without_DC_2, data.signals.values(:,1), H3, vector, Ts);
grid on
bodeplot(tf(num,den), 'r')
bodeplot(tf_3, 'b')
grid on
title(['H = ',num2str(H3)])
legend('Idealno', 'Estimirano')

figure(14)
hold on
bodeplot(tf(num,den),tf_1,tf_2,tf_3)
grid on
title('Bodeov Dijagram')
legend('Idealno',['H = ',num2str(H1)],['H = ',num2str(H2)], ['H = ',num2str(H3)])

%% c)
%% 

%ulaz
figure(15)
subplot(2,1,1)
plot(abs(fft(without_DC_1)));
title('Ulaz DFT-a')
ylabel('u')
xlabel('N')
grid on

%izlaz
subplot(2,1,2)
plot(abs(fft(without_DC_2)));
title('Izlaz DFT-a')
ylabel('y')
xlabel('N')
grid on

%po formuli
G_jw= fft(without_DC_2)./fft(without_DC_1);
G_jw = idfrd(G_jw, linspace(0.01*2*pi,length(without_DC_1)*0.1,length(without_DC_1)), Ts);

figure(16)
hold on
bodeplot(G_jw,tf(num,den))
title('Frekvencijska karakteristika idealnog modela i estimiranog modela pomoću DFT-a')
legend('DTF', 'Idealno')
grid on


