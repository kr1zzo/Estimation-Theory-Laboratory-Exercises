
function [sys,x0,str,ts] = riv_sfun(t,x,u,flag,Ts)
% Dispatch the flag. The switch function controls the calls to 
% S-function routines at each simulation stage.
switch flag,
case 0
[sys,x0,str,ts] = mdlInitializeSizes(Ts); % Initialization
case 3
sys = mdlOutputs(t,x,u,flag); % Calculate outputs
case { 1, 2, 4, 9 }
   sys = []; % Unused flags
otherwise
error(['Unhandled flag = ',num2str(flag)]); % Error handling
end;
% End of function timestwo

%============================================================== 
% Function mdlInitializeSizes initializes the states, sample 
% times, state ordering strings (str), and sizes structure.
%==============================================================
function [sys,x0,str,ts] = mdlInitializeSizes(Ts)
% Call function simsizes to create the sizes structure.
sizes = simsizes;
% Load the sizes structure with the initialization information.
sizes.NumContStates= 0;
sizes.NumDiscStates= 0;
sizes.NumOutputs= 4;
sizes.NumInputs= 2;
sizes.DirFeedthrough=1;
sizes.NumSampleTimes=1;
% Load the sys vector with the sizes information.
sys = simsizes(sizes);
%
x0 = []; % No continuous states
% 
str = []; % No state ordering
% 
ts = [Ts 0]; % Inherited sample time
% End of mdlInitializeSizes.
%==============================================================
% Function mdlOutputs performs the calculations.
%==============================================================
function sys = mdlOutputs(t,x,u,flag)

% Vektori theta, fi, thetaH, w i matrica P moraju se pamtiti iz koraka u korak, stoga
% su definirani kao persistent varijable
% thetaK_ - vektor estimiranih podataka u k-1 koraku
% fiK - regresijski vektor u k-tom koraku
% PK_ - matrica kovarijance procjene parametara u k-1 koraku
% thetaHK_ - vektor parametara pomoænog modela u k-1 koraku
% wK_ - vektor pomoænih varijabli u k-1 koraku
% thetaKvec - ovu matricu možete koristiti za spremanje prošligh
% vrijednosti estimiranih parametera
persistent PK_ thetaK_ fiK wK_ thetaHK_ thetaKvec 

%% 
%Inicijalizacija - potrebno je upisati poèetne vrijednosti matrica i
%vektora
if (t==0)
    PK_ = 500*eye(4,4);
    thetaK_ = [0 0 0 0]';
    thetaHK_ = [1 0 0 0]';
    fiK = [0 0 0 0]';
    wK_ = [0 0 0 0]';
    thetaKvec = zeros(4,10);
   
end
%% 


% Mjereni podaci ulaza i izlaza:
% yK - mjerenje izlaza procesa u k-tom trenutku
% uK - mjerenje ulaza u proces u k-tom trenutku
yK=u(2);
uK=u(1);

%_____________________________________________________
% Nadopunite funkciju tako da ostvaruje RIV algoritam
%_____________________________________________________

global phi1;

phi = phi1;
% brzina filtra
gamma = 0.01; 

%korak
y_ = wK_' * thetaHK_;
%w matrica
wK_ = [-y_ wK_(1) fiK(3) fiK(4)]';
%poajcanje
K = (PK_ * wK_)*pinv(1 + (fiK') * PK_ * wK_);
v = yK - (fiK') * thetaK_;
thetaK_ = thetaK_ + K * v;
PK_ = (PK_ - K * (fiK') * PK_)/phi;
fiK = [-yK fiK(1) uK fiK(3)]';
thetaHK_ = (1 - gamma) * thetaHK_ + gamma * thetaKvec(:,10);
thetaKvec = [thetaK_ thetaKvec(:,1:end-1)];

sys = [thetaK_];


% End of mdlOutputs