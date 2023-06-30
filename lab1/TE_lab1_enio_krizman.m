%1.a)
close all

syms t
%TO-DO. Generate 3x9 matrix H_ti. This matrix can be defined parametrically
%with regards to parameter t, for example H_ti = [cos(t), 0; 0, sin(t)];

H_ti = [1, 0, 0, cos(t), 0, 0, sin(t), 0, 0;
        0, 1, 0, 0, cos(t), 0, 0, sin(t), 0;
        0, 0, 1, 0, 0, cos(t), 0, 0, sin(t)];
p0 = [200; 300; 300];
d1 = [0; 200; 0];
d2 = [100; 0; 200];
x = [p0; d1; d2];

% Ellipse visualization
figure;
%Generates measurements without noise, i.e. covariance matrix has 0
%variance
[~, ~, z] = genMeasurements(H_ti, x, linspace(0, 2*pi, 50), t, diag([0, 0, 0]));
plot3(z(1:3:end), z(2:3:end), z(3:3:end))
hold on;
%Generates measrements with noise as defined in 1.a), i.e. Sigma_ti =
%diag([4, 4, 4]). Note that Sigma_ti can also be defined parametrically,
%for example diag([t, 2*t, t])
[~, ~, z] = genMeasurements(H_ti, x, linspace(0, 2*pi, 50), t, diag([4, 4, 4]));
scatter3(z(1:3:end), z(2:3:end), z(3:3:end))
xlabel('x')
ylabel('y')
zlabel('z')
title('Measured and real trajectory')
legend('Real trajectory', 'Measured trajectory')

%% 1.b)
close all
syms t

N = linspace(5,50,10);
H_ti = [1, 0, 0, cos(t), 0, 0, sin(t), 0, 0;
        0, 1, 0, 0, cos(t), 0, 0, sin(t), 0;
        0, 0, 1, 0, 0, cos(t), 0, 0, sin(t)];
p0 = [200; 300; 300];
d1 = [0; 200; 0];
d2 = [100; 0; 200];
x = [p0; d1; d2];

x_least_square = zeros(9, 50);
varianca= zeros(1, length(N));
CRLB = zeros(1, length(N));
%counter = 0;
%m_counter = 0;

for M=[5,50]
    
    for i=1:length(N)
        %counter = counter+1;
        P = 0;
        for j = 1:M
            %repeat least square estimaton process for M times
            [H, sigma, z] = genMeasurements(H_ti, x, linspace(0, 2*pi, N(i)), t, diag([2, 2, 2]));
            %least square optimal x
            x_least_square = pinv(transpose(H) * H) * transpose(H) * z;
            %variance P
            P = P + transpose(x_least_square - x)*(x_least_square - x);
        end
        
        varianca(i) = P/M;
        
        %CRAMER-RAO LOWER BOUND - smallest possible uncertainity
        %if estimator variance iz equal to CRLB then estimator is efficient
        CRLB(i) = trace(pinv(transpose(H) * pinv(sigma) * H));

        %display(counter)
         
    end

    %fprintf('izaso')

    plot(N, varianca);
    hold on

    %m_counter = m_counter+1;
    %display(m_counter)
    
end

plot(N, CRLB);
hold on;
grid on;
xlabel('N')
ylabel('P')
legend('M=5', 'M=50', 'CLRB')
hold off;
%% 1. C)
close all
syms t

H_ti = [1, 0, 0, cos(t), 0, 0, sin(t), 0, 0;
        0, 1, 0, 0, cos(t), 0, 0, sin(t), 0;
        0, 0, 1, 0, 0, cos(t), 0, 0, sin(t)];
p0 = [200; 300; 300];
d1 = [0; 200; 0];
d2 = [100; 0; 200];
x = [p0; d1; d2];

x_least_square = zeros(9, 50);
varianca= zeros(1, length(N));
CRLB = zeros(1, length(N));
%counter = 0;
%m_counter = 0;

sigma_z = linspace(1,50,10);

for i=1:10
    %counter = counter+1;
    P = 0;
    for j = 1:10
        %repeat least square estimaton process for M times
        [H, sigma, z] = genMeasurements(H_ti, x, linspace(0, 2*pi, 10), t, diag([sigma_z(i), sigma_z(i), sigma_z(i)]));
        %least square optimal x
        x_least_square = pinv(transpose(H) * H) * transpose(H) * z;
        %variance P
        P = P + transpose(x_least_square - x)*(x_least_square - x);
    end
        
    varianca(i) = P/10;
    
    %CRAMER-RAO LOWER BOUND - smallest possible uncertainity
    %if estimator variance iz equal to CRLB then estimator is efficient
    CRLB(i) = trace(pinv(transpose(H) * pinv(sigma) * H));

    %display(counter)
end

plot(sigma_z, CRLB);
hold on;
plot(sigma_z, varianca);
hold on;
grid on;
xlabel('σ_z')
ylabel('P')
legend('CRLB', 'M = 10') 
hold off;
%% 2. a) first part
close all
syms t

% Ellipse visualization
figure;
%Generates measurements without noise, i.e. covariance matrix has 0
%variance
[~, ~, z] = genMeasurements(H_ti, x, linspace(0, 2*pi, 50), t, diag([0, 0, 0]));
plot3(z(1:3:end), z(2:3:end), z(3:3:end))
hold on;
%Generates measrements with noise as defined in 1.a), i.e. Sigma_ti =
%diag([4, 4, 4]). Note that Sigma_ti can also be defined parametrically,
%for example diag([t, 2*t, t])
[~, ~, z] = genMeasurements(H_ti, x, linspace(0, 2*pi, 50), t, diag([exp(t), exp(t), exp(t)]));
scatter3(z(1:3:end), z(2:3:end), z(3:3:end))
xlabel('x')
ylabel('y')
zlabel('z')
title('Measured and real trajectory')
legend('Real trajectory', 'Measured trajectory')

%% 2 a) second part
close all
syms t

H_ti = [1, 0, 0, cos(t), 0, 0, sin(t), 0, 0;
        0, 1, 0, 0, cos(t), 0, 0, sin(t), 0;
        0, 0, 1, 0, 0, cos(t), 0, 0, sin(t)];
p0 = [200; 300; 300];
d1 = [0; 200; 0];
d2 = [100; 0; 200];
x = [p0; d1; d2];

N = linspace(5, 50, 10);
x_least_square = zeros(9, 50);
varianca= zeros(1, length(N));
CRLB = zeros(1, length(N));

for i=1:length(N)
    P = 0;
    for j=1:20
        %repeat least square estimaton process for M times
        [H, sigma, z] = genMeasurements(H_ti, x, linspace(0, 2*pi, N(i)), t, diag([exp(t), exp(t), exp(t)]));
        %least square optimal x
        x_least_square = pinv(transpose(H) * H) * transpose(H) * z;
        %variance P
        P = P + transpose(x_least_square - x)*(x_least_square - x);
    end

    varianca(i) = P/20;
    
    %CRAMER-RAO LOWER BOUND - smallest possible uncertainity
    %if estimator variance iz equal to CRLB then estimator is efficient
    CRLB(i) = trace(pinv(transpose(H) * pinv(sigma) * H));
end

plot(N, CRLB);
hold on;
plot(N, varianca);
hold on;
grid on;
xlabel('N')
ylabel('P')
legend('CRLB', 'M = 20') 
hold off;
%% 2. b)
close all
syms t

H_ti = [1, 0, 0, cos(t), 0, 0, sin(t), 0, 0;
        0, 1, 0, 0, cos(t), 0, 0, sin(t), 0;
        0, 0, 1, 0, 0, cos(t), 0, 0, sin(t)];
p0 = [200; 300; 300];
d1 = [0; 200; 0];
d2 = [100; 0; 200];
x = [p0; d1; d2];

N = linspace(5, 50, 10);
x_weighted_least_square = zeros(9, 50);
varianca= zeros(1, length(N));
CRLB = zeros(1, length(N));

for i=1:length(N)
    P = 0;
    for j=1:20
        [H, sigma, z] = genMeasurements(H_ti, x, linspace(0, 2*pi, N(i)), t, diag([exp(t), exp(t), exp(t)]));
        x_weighted_least_square = pinv(transpose(H)*pinv(sigma) * H) * transpose(H)*pinv(sigma) * z;
        %variance P
        P = P + transpose(x_weighted_least_square - x)*(x_weighted_least_square - x);
    end

    varianca(i) = P/20;
    
    %CRAMER-RAO LOWER BOUND - smallest possible uncertainity
    %if estimator variance iz equal to CRLB then estimator is efficient
    CRLB(i) = trace(pinv(transpose(H) * pinv(sigma) * H));
end

plot(N, CRLB);
hold on;
plot(N, varianca);
hold on;
grid on;
xlabel('N')
ylabel('P')
legend('CRLB', 'M = 20') 
hold off;
%% 3. a)
close all
syms t

H_ti = [1, 0, 0, cos(t), 0, 0, sin(t), 0, 0;
        0, 1, 0, 0, cos(t), 0, 0, sin(t), 0;
        0, 0, 1, 0, 0, cos(t), 0, 0, sin(t)];
p0 = [200; 300; 300];
d1 = [0; 200; 0];
d2 = [100; 0; 200];
x_mean = [p0; d1; d2];
s = linspace(0.01, 5, 40);
N = 10;
M = 50;

varianca_LS = zeros(1, length(s));
varianca_MMSE = zeros(1, length(s));

for i=1:length(s)
    Px = s(i)*eye(9);
    P1 = 0;
    P2 = 0;
    for j=1:M
        % X = N(x_mean, Px)
        x(:,i) =  mvnrnd(x_mean,Px);

        [H, sigma, z] = genMeasurements(H_ti, x(:,i), linspace(0, 2*pi, N), t, diag([2,2,2]));
        %MMSE estimator
        x_MMSE = x_mean + Px*transpose(H)*pinv(H*Px*transpose(H)+sigma)*(z-H*x_mean);
        P1 = P1 + transpose(x_MMSE - x(:,i))*(x_MMSE - x(:,i));
        %least square optimal x
        x_least_square = pinv(transpose(H) * H) * transpose(H) * z;
        %variance P
        P2 = P2 + transpose(x_least_square - x(:,i))*(x_least_square - x(:,i));
    end

    varianca_MMSE(i) = P1/50;
    varianca_LS(i) = P2/50;
end

plot(s, varianca_MMSE);
hold on;
plot(s, varianca_LS);
hold on;
grid on;
xlabel('N')
ylabel('P')
legend('MMSE', 'LS') 
hold off;
%% 3. b)
close all
syms t

H_ti = [1, 0, 0, cos(t), 0, 0, sin(t), 0, 0;
        0, 1, 0, 0, cos(t), 0, 0, sin(t), 0;
        0, 0, 1, 0, 0, cos(t), 0, 0, sin(t)];
p0 = [200; 300; 300];
d1 = [0; 200; 0];
d2 = [100; 0; 200];
x_mean = [p0; d1; d2];
s_z = linspace(0.01, 10, 40);
N = 10;
M = 50;

varianca_LS = zeros(1, length(s));
varianca_MMSE = zeros(1, length(s));

for i=1:length(s)
    Px = 2*eye(length(x_mean));
    P1 = 0;
    P2 = 0;
    for j=1:M
        % X = N(x_mean, Px)
        x(:,i) =  mvnrnd(x_mean,Px);

        [H, sigma, z] = genMeasurements(H_ti, x(:,i), linspace(0, 2*pi, N), t, diag([s_z(i), s_z(i), s_z(i)]));
        %MMSE estimator
        x_MMSE = x_mean + Px*transpose(H)*pinv(H*Px*transpose(H)+sigma)*(z-H*x_mean);
        P1 = P1 + transpose(x_MMSE - x(:,i))*(x_MMSE - x(:,i));
        %least square optimal x
        x_least_square = pinv(transpose(H) * H) * transpose(H) * z;
        %variance P
        P2 = P2 + transpose(x_least_square - x(:,i))*(x_least_square - x(:,i));
    end

    varianca_MMSE(i) = P1/50;
    varianca_LS(i) = P2/50;
end

plot(s, varianca_MMSE);
hold on;
plot(s, varianca_LS);
hold on;
grid on;
xlabel('N')
ylabel('P')
legend('MMSE', 'LS') 
hold off;
%% 

function [H, Sigma, z] = genMeasurements(H_ti, x, values, t, Sigma_ti)
%GENMEASUREMENTS Generates measurements and appropriate matrices given in equation (3) of task 1

%Function inputs are: H_ti as in equation (2). This is a parametrically defined linear system matrix FOR A SINGLE MEASUREMENT which can depend on parameter t.
%                     Sigma_ti as in equation (2): Parametrically defined normal distribution covariance matrix FOR A SINGLE MEASUREMENT.
%                     t: Symbolic matlab variable which can be defined with function syms
%                     values: Vector of values [t1, ..., tN] for N measurements. 
%                     x: Parameter vector being estimated
%
%Function outputs are: z, H, and Sigma as defined in (3)
rows = size(H_ti, 1);
cols = size(H_ti, 2);
N = length(values);
H = zeros(N * rows, cols);
Sigma = zeros(N * rows, N * rows);

for i=1:length(values)
   H(rows*(i-1)+1:rows*(i), 1:cols) = subs(H_ti, t, values(i));
   Sigma(rows*(i-1)+1:rows*(i), rows*(i-1)+1:rows*(i)) = subs(Sigma_ti, t, values(i));
end

%C_w_stacked = diag(repmat(diag(C_w), N, 1));
w = mvnrnd(zeros(1, N*rows), Sigma);
z = H * x + w.';
end