function price = BinomialTreeCRR(S0, K, r, T, sigma, div, N)
% Binomial Tree Model of CRR
dt = T / N;
u = exp(sigma * sqrt(dt));
d = 1 / u;
p = (exp((r - div) * dt) - d) / (u - d);
discount = exp(- r * dt);
pu = discount * p;
pd = discount * (1 - p);
SVals = zeros(2 * N + 1, 1);
SVals(N + 1) = S0;
for i = 1:N
    SVals(N + 1 + i) = u * SVals(N + i);
    SVals(N + 1 - i) = d * SVals(N - i + 2);
end
PVals = zeros(2 * N + 1, 1);
for i = 1:2:2*N+1
    PVals(i) = max(K - SVals(i), 0);
end
for t = 1:N
    for i = (t+1):2:(2*N+1-t)
        hold = pu * PVals(i + 1) + pd * PVals(i - 1);
        PVals(i) = max(hold, K - SVals(i));
    end
end
price = PVals(N+1);
% price = BinomialTreeCRR(50,50,0.05,5/12,0.4,0,1000)
end

