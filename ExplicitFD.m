function price = ExplicitFD(S0, K, r, T, sigma, div, Smax, dS, dt)
M = round(Smax/dS);
dS = Smax/M;
N = round(T/dt);
dt = T/N;
matval = zeros(M+1, N+1);
vetS = linspace(0, Smax, M+1);
veti = 0:M;
vetj = 0:N;
% set up boundary conditions
matval(:, N+1) = max(K-vetS, 0);
matval(1, :) = K * exp(- r * dt * (N - vetj));
matval(M+1, :) = 0;
% set up coefficients
a = 0.5 * dt * (sigma^2 * veti - r + div) .* veti;
b = 1 - dt * (sigma^2 * veti.^2 + r);
c = 0.5 * dt * (sigma^2 * veti + r - div) .* veti;
% solve backward in time
for j=N:-1:1
    for i=2:M
        hold = a(i)*matval(i-1,j+1)+b(i)*matval(i,j+1)+c(i)*matval(i+1,j+1);
        matval(i,j) = max(hold, K-vetS(i));
    end
end
% return price, possibly by linear interpolation outside the grid
price = interp1(vetS, matval(:,1), S0);
% price = ExplicitFD(50,50,0.05,5/12,0.4,0,100,2,5/1200)
end