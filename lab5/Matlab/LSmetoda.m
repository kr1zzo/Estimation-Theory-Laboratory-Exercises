function [hat_theta,MSE]=LSmetoda(u,y,n)
%% 
% u-ulazni signal (Nx1)
% y-izlazni signal (Nx1)
% n-pretpostavljeni red modela

% hat_theta-vektor estimiranih parametara
% MSE-srednja kvadratna pogreska estimiranog modela
%% 
N = length(y);

for i= 1:(N-n)
    for j = 1:n
        phi_y(i,j) = -y(n+i-j);
        phi_u(i,j) = u(n+i-j);
    end
end

%matrica poadataka
phi = [phi_y phi_u]; 
%izlazni podaci
y_ = y((n+1):end); 

hat_theta = pinv(phi.' * phi)* phi.' * y_;

MSE = ((y_- phi * hat_theta).' * (y_- phi * hat_theta))/(N-n);

end


