function [eps] = FMM_happy_epsilon_new(epsilon, P,Q,R)

%P = 2*N+1;
%Q = 2*N+1;
%R = 1;
[Nx,Ny] = size(epsilon);
Nz = 1;
NH = P*Q*R;
p = [-floor(P/2):+floor(P/2)];
q = [-floor(Q/2):+floor(Q/2)];
r = [-floor(R/2):+floor(R/2)];

A = fftshift(fftn(epsilon))/(Nx*Ny*Nz);
p0 = 1+floor(Nx/2);
q0 = 1+floor(Ny/2);
r0 = 1+floor(Nz/2);

eps = zeros(P*Q,P*Q);
for rrow = 1:R
for qrow = 1:Q
for prow = 1:P
    row = (rrow-1)*Q*P + (qrow-1)*P + prow;
    for rcol = 1:R
    for qcol = 1:Q
    for pcol = 1:P
        col = (rcol-1)*Q*P + (qcol-1)*P + pcol;
        pfft = p(prow) - p(pcol);
        qfft = q(qrow) - q(qcol);
        rfft = r(rrow) - r(rcol);
        eps(row,col) = A(p0+pfft, q0+qfft, r0+rfft);
    end
    end
    end
end
end
end