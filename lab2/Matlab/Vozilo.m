function [xhat, P, xRMSE] = Vozilo(gps_type)

% Pracenje vozila u 2D
%
% Ulazi su:
%   gps_type = tip koristenih mjerenja
%          'gps_easy' - osnovni skup mjerenja
%          'gps_hard' - prošireni skup mjerenja s velikim pogreškama

% Izlazi su:
%   Estimirano stanje xhat = [x;, y;, xdot; ydot]
%   Kovarijanca pogreske estimacije P
%   XRMSE = srednja kvadraticna pogreska estimacije lokacije vozila

%% 

% ovdje opisati sustav u prostoru stanja: postaviti matrice Phi, Gamma, H, L, M
syms T
T = 1/ 10

A = [1, 0, T, 0; 
    0, 1, 0, T; 
    0, 0, 1, 0; 
    0, 0, 0, 1];

B = [0;0;0;0];

H = [1, 0, 0, 0; 
    0, 1, 0, 0];

L = eye(4);
M = eye(2);


% ovdje postaviti matice kovarijanci procesnog i mjernog suma Q i R
Q = 0.1*eye(4);
R = 4*eye(2);

% ucitaj odgovarajuci skup mjerenja
y = [];
if strcmp(gps_type, 'gps_easy')
    load gps_easy gps_y_easy x_gt
    y = gps_y_easy;
elseif strcmp (gps_type, 'gps_hard')
    load gps_hard gps_y_hard x_gt
    y = gps_y_hard;
else
    error(['Skup mjerenja ' gps_type ' nije definiran.']);
end


xhat = [y(1,1); y(1, 2); 0; 0]; 

P = [R(1,1),0, 0, 0; 
    0, R(1, 2), 0, 0; 
    0, 0, 1500, 0; 
    0, 0, 0, 1500]; 

xhatArray = [];
PtraceArray = [];
xhatErrArray = [];
%% SIMULACIJA KALMANOVOG FILTRA
for k = 1 : size(y,1)

   % Ovdje napisati kod odgovarajuceg  Kalmanovog filtra
   
   %prediction equations:
   P = A*P *transpose(A)+Q;
   
   K_k =P*transpose(H) * inv(H*P*transpose(H)+R);
   
   xhat = A*xhat;
   
   err =  transpose(xhat(1:2)) - x_gt(k,1:2) ; 
    
   if ~(isnan(y(k,1)))
        xhat = xhat + K_k*(y(k,:)' - H * xhat);
   end
   
   %update equation:
   P = (eye(4) - K_k*H) * P *transpose(eye(4) - K_k*H) + K_k*R*transpose(K_k);
   
   %Spremi podatke za iscrtavanje
   xhatArray = [xhatArray; xhat(1:2)'];
   PtraceArray = [PtraceArray; trace(P)] ;
   xhatErrArray = [xhatErrArray; err];
end

xRMSE = rms(xhatErrArray);
if strcmp(gps_type, 'gps_easy')
    yRMSE = rms(y-x_gt);
    disp(['RMSE mjerenja lokacije vozila  ',  num2str(yRMSE) ' m']);
end
disp(['RMSE estimacije lokacije vozila  ',  num2str(xRMSE) ' m']);
%% Iscrtavanje
figure; hold on;
title('Trajektorija vozila i mjerenja senzora');
hold on;
grid on;
plot(xhatArray(:,1), xhatArray(:,2), 'g');
plot(x_gt(:,1), x_gt(:,2));
plot(y(:,1), y(:,2), 'xr');
legend('Estimirana trajektorija vozila', 'Točna trajektorija vozila', 'Mjerenje senzora')

figure; 
title('Trag matrice P');
hold on; 
grid on;
plot(PtraceArray);
figure;
title('Korijen kvadrata pogreske estimacije');
hold on; 
grid on;
plot(xhatErrArray);
