[xhat6, P6, t6, err6] = Satelit('EKF', 5);

[xhat0, P0, t0, err0] = Satelit('IEKF', 0);
[xhat1, P1, t1, err1] = Satelit('IEKF', 1);
[xhat2, P2, t2, err2] = Satelit('IEKF', 2);
[xhat3, P3, t3, err3] = Satelit('IEKF', 3);
[xhat4, P4, t4, err4] = Satelit('IEKF', 4);
[xhat5, P5, t5, err5] = Satelit('IEKF', 5);

set(gca,'FontSize',12); set(gcf,'Color','White');
grid on
xlabel('[sekunde]'); ylabel('Greska estimacije radiusa [m]');
grid on;

figure;
plot(t0, err0);
hold on;
plot(t1, err1);
hold on;
plot(t2, err2);
hold on;
plot(t3, err3);
hold on;
plot(t4, err4);
hold on;
plot(t5, err5);
hold on;
plot(t6, err6);
hold on;
legend('IEKF(Niter = 0)','IEKF(Niter = 1)','IEKF(Niter = 2)','IEKF(Niter = 3)','IEKF(Niter = 4)','IEKF(N_iter = 5)','EKF')
hold off;

figure;
plot(t6, err6);
hold on;
plot(t0, err0);
legend('EKF','IEKF');
title('Niter = 0');
hold off;

figure;
plot(t6, err6);
hold on;
plot(t5, err5);
legend('EKF','IEKF');
title('Niter = 5');
hold off;