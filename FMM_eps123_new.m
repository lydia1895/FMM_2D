function [eps11, eps22, eps33] = FMM_eps123_new(epsilon,N,M)

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

NN = 2*N+1;
eps1inv2=zeros(NN,NN,NN,NN);    %[inv([eps]1)]2
epsinv1inv2=zeros(NN,NN,NN,NN); %[inv([inv(eps)]1)]2
X1t=zeros(M,1);
Z1t=zeros(M,1);
for i=1:(2*N+1)
    for j=1:(2*N+1)
        for k=1:M
        X1t(k) = X1(i,j,k);
        Z1t(k) = Z1(i,j,k);
        end
        n = size(X1t);
        eps1inv2(i,j,:,:) = FMM_happy_epsilon_new(X1t, P,R,R);
        epsinv1inv2(i,j,:,:) = FMM_happy_epsilon_new(Z1t, P,R,R);
    end
end

NH = P*Q*R;
P = 2*N+1;
Q = 2*N+1;
R = 1;

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
        eps22inv(row,col) = eps1inv2(qrow,qcol,prow,pcol);
        eps11(row,col) = epsinv1inv2(qrow,qcol,prow,pcol);
    end
    end
    end
end
end
end
eps22 = eye(P*Q,P*Q)/eps22inv;

