function FMM_graph_output(R_T, alpha0, beta0, N, periodx, periody)

NN = (2*N+1)*(2*N+1);
R1 = R_T(1:NN);
R2 = R_T((NN+1):2*NN);
T1 = R_T((2*NN+1):3*NN);  %d0
T2 = R_T((3*NN+1):4*NN);

alpha_mm = zeros(2*N+1,1);
beta_mm = zeros(2*N+1,1);

for m=1:(2*N+1)
    alpha_mm(m) = (m-N-1)*2*pi/periodx;
    beta_mm(m) =  (m-N-1)*2*pi/periody;
end

alpha_p = alpha0 + alpha_mm;
beta_p =  beta0 + beta_mm;

n_points_FMM_x = 200;
n_points_FMM_y = 200;

reflected_field = zeros(n_points_FMM_x,n_points_FMM_y);
x_new = linspace(0, periodx, n_points_FMM_x);
y_new = linspace(0, periody, n_points_FMM_y);


for i = 1:n_points_FMM_x
for j = 1:n_points_FMM_y
    for m=1:(2*N+1)
    for n=1:(2*N+1)
        row = n + (m-1)*(2*N+1);
        reflected_field(i,j) = reflected_field(i,j)+...
            real(R1(row)*exp(1j*alpha_p(m)*x_new(i))*exp(1j*beta_p(n)*y_new(j)));
    end
    end
end
end

figure(3);
pcolor(x_new, y_new, transpose(reflected_field))
shading interp
%caxis([-1.0 1.0])
