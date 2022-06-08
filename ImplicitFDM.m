function price = ImplicitFDM(S0, K, r, T, sigma, div, Smax, dS, dt)
M = round(Smax / dS);
dS = Smax / M;
N = round(T / dt);
dt = T / N;
vetS = linspace(0, Smax, M+1);
% set up boundary conditions
B = max(K-vetS', 0);
% set up coefficients
A = zeros(M+1,M+1);
A(1,1) = 1;
for i=2:M
    A(i,i-1) = 0.5 * dt * (r - div - sigma^2 * (i-1)) * (i-1);
    A(i,i) = 1 + dt * (sigma^2 * (i-1)^2 + r);
    A(i,i+1) = 0.5 * dt * (- r + div - sigma^2 * (i-1)) * (i-1);
end
A(M+1,M+1) = 1;
F = inv(A) * B;
% solve backward in time
for j = N:-1:1
    B = F;
    F = inv(A) * B;
    for i = 2:M
        F(i) = max(F(i), K-vetS(i));
    end
end
% return price, possibly by linear interpolation outside the grid
price = F(M/2+1);
% price = ImplicitFD(50,50,0.05,5/12,0.4,0,100,2,5/1200)
end

