function [sys,x0,str,ts] = model_sfun(t,x,u,flag,theta,theta2)



switch flag,

  case 0
    [sys,x0,str,ts]=mdlInitializeSizes(); % Initialization

  case 1
    sys = mdlDerivatives(t,x,u); % Calculate derivatives

  case 3
    sys = mdlOutputs(t,x,u,theta,theta2); % Calculate outputs
  
  case { 2, 4, 9 } % Unused flags
    sys = [];
  otherwise
    error(['Unhandled flag = ',num2str(flag)]); % Error handling
end
% End of csfunc.
%==============================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the 
% S-function.
%==============================================================
%
function [sys,x0,str,ts] = mdlInitializeSizes()
%
% Call simsizes for a sizes structure, fill it in and convert it 
% to a sizes array.
%
sizes = simsizes;
sizes.NumContStates  = 2;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 5;
sizes.NumInputs      = 1;
sizes.DirFeedthrough = 1;     % Matrix D is nonempty. 
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
%
% Initialize the initial conditions.
%
x0 = zeros(2,1);
%
% str is an empty matrix.
%
str = [];
%
% Initialize the array of sample times; in this example the sample 
% time is continuous, so set ts to 0 and its offset to 0.
%
ts = [0 0];
% End of mdlInitializeSizes.
%
%==============================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%==============================================================
function sys = mdlDerivatives(t,x,u)

if(t<300)
m=3;
c=1;
k=100;
A=[0 1; -k/m -c/m];
B=[0; 1/m];
else
    m=2;
c=1;
k=100;
A=[0 1; -k/m -c/m];
B=[0; 1/m];
end

sys = A*x + B*u;
% End of mdlDerivatives.
%
%==============================================================
% mdlOutputs
% Return the block outputs.
%==============================================================
%
function sys = mdlOutputs(t,x,u,theta,theta2)

if(t<300)
par=theta;
else
par=theta2;
end

C=[1 0];
D=[0];
sys = [C*x + D*u, par'];
% End of mdlOutputs.