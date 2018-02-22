function [eps11, eps22, eps33] = FMM_eps123(epsilon,N,M)

P = 2*N+1;
Q = 2*N+1;
R = 1;

%eps33
eps33 = FMM_happy_epsilon_new(epsilon, P,Q,R);

%eps22, eps11

epsinv = ones(M,M)./epsilon;
epsinv1 = zeros(2*N+1,2*N+1,M);
Z1 = zeros(2*N+1,2*N+1,M);

eps1 = zeros(2*N+1,2*N+1,M);
X1 = zeros(2*N+1,2*N+1,M);
for i=1:M
    eps1(:,:,i) = FMM_happy_epsilon_new(transpose(epsilon(i,:)), P,R,R);
    X1(:,:,i) = eye(2*N+1,2*N+1)/eps1(:,:,i);
    
    epsinv1(:,:,i) = FMM_happy_epsilon_new(transpose(epsinv(i,:)), P,R,R);
    Z1(:,:,i) = eye(2*N+1,2*N+1)/epsinv1(:,:,i);
end

X1line = zeros(4*N+1,M);%zeros(M,M);
Z1line = zeros(4*N+1,M);
for i =1:(2*N+1)
    for j=1:M
    X1line(j,2*N+i) = X1(i,1,j);   %for eps22
    X1line(j,2*N-i+2) = X1(1,i,j);
    
    Z1line(j,2*N+i) = Z1(i,1,j);   %for eps11
    Z1line(j,2*N-i+2) = Z1(1,i,j);
    end
end

eps1inv2=zeros(M,4*N+1);    %[inv([eps]1)]2
epsinv1inv2=zeros(M,4*N+1); %[inv([inv(eps)]1)]2

for i=1:(4*N+1)
    eps1inv2(:,i) = fftshift(fftn(X1line(:,i)))/M;
    epsinv1inv2(:,i) = fftshift(fftn(Z1line(:,i)))/M;
end
n0 = 1+floor(M/2);
n = [-2*N:2*N];
eps1inv2n = zeros(4*N+1,4*N+1);
epsinv1inv2n = zeros(4*N+1,4*N+1);
for i=1:(4*N+1)
    eps1inv2n(i,:)=eps1inv2(n0+n(i),:);
    epsinv1inv2n(i,:) = epsinv1inv2(n0+n(i),:);
end

NH = P*Q*R;
P = 2*N+1;
Q = 2*N+1;
R = 1;
p = [-floor(P/2):+floor(P/2)];
q = [-floor(Q/2):+floor(Q/2)];
r = [-floor(R/2):+floor(R/2)];

[Nx, Ny] = size(eps1inv2n);
Nz = 1;
p0 = 1+floor(Nx/2);
q0 = 1+floor(Ny/2);
r0 = 1+floor(Nz/2);

eps22inv = zeros(P*Q,P*Q);
eps11 = zeros(P*Q,P*Q);
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
        eps22inv(row,col) = eps1inv2n(p0+pfft, q0+qfft, r0+rfft);
        eps11(row,col) = epsinv1inv2n(p0+pfft, q0+qfft, r0+rfft);
    end
    end
    end
end
end
end
eps22 = eye(P*Q,P*Q)/eps22inv;