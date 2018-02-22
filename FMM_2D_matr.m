function [M, G] = FMM_2D_matr(N, alpha, beta, eps, k0)

miden = eye((2*N+1)*(2*N+1),(2*N+1)*(2*N+1));

eps_inv = miden/eps;
F11 = alpha*eps_inv*beta;
F12 = miden*k0^2 - alpha*eps_inv*alpha;
F21 = beta*eps_inv*beta - miden*k0^2;
F22 = -beta*eps_inv*alpha;

F_up = cat(2,F11,F12);
F_down = cat(2,F21,F22);
F = cat(1,F_up, F_down);

G11 = -alpha*beta;
G12 = alpha^2 - (k0^2)*eps;
G21 = (k0^2)*eps - beta^2;
G22 = beta*alpha;
G_up = cat(2,G11,G12);
G_down = cat(2,G21,G22);
G = cat(1,G_up, G_down);

M = F*G/(k0^2);
%{
fcond = rcond(F)
gcond = rcond(G)
Mcond = rcond(M)
%}