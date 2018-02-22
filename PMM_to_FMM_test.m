function [La, alpha_ref, beta_ref, b_x1, b_x2, N_intervals_x, N_intervals_y,...
    N_basis_x, N_basis_y, Nx, nx, Ny, ny] = PMM_to_FMM_test

N_intervals_x = 2;
N_intervals_y = 2;
N_b = 9;
N_basis_x = N_b*ones(N_intervals_x,1);
N_basis_y = N_b*ones(N_intervals_y,1);

%boundaries in the periodic layer
%for ellipse here will be matched coordinates

lambda = 1000;

b_x1 = [0 0.25*lambda 0.9*lambda];    
b_x2 = [0 0.25*lambda 0.9*lambda];

[Nxx, NNxx] = size(b_x1);
[Nyy, NNyy] = size(b_x2);
periodx = b_x1(NNxx)-b_x1(1);
periody = b_x2(NNyy)-b_x2(1);

nx = N_basis_x-1;   %because functions are p(0)...p(n(k)) on interval k
Nx = zeros(N_intervals_x,1);
for k=1:N_intervals_x   
    Nx(k) = -nx(k);  
    for p=1:k
        Nx(k) = Nx(k)+nx(p);
    end
end

ny = N_basis_y-1;
Ny = zeros(N_intervals_y,1);
for k=1:N_intervals_y   
    Ny(k) = -ny(k);  
    for p=1:k
        Ny(k) = Ny(k)+ny(p);
    end
end

alpha_ref = pi/(4*periodx);
beta_ref = pi/(4*periody);

%tau_x = exp(1j*alpha_ref*periodx);
%tau_y = exp(1j*beta_ref*periody);

%La is lambda in Gegenbauer polynomials; La>-1/2

La = 0.5;
