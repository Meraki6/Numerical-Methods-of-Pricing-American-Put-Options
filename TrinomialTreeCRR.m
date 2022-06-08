function price = TrinomialTreeCRR(S0, K, r, T, sigma, div, N)
% Trinomial Tree Model
dt = T / N;
u = exp(sigma * sqrt(2 * dt));
d = 1 / u;
pu = ((exp((r-div)*dt/2)-exp(-sigma*sqrt(dt/2)))/(exp(sigma*sqrt(dt/2))-exp(-sigma*sqrt(dt/2))))^2;
pd = ((exp(sigma*sqrt(dt/2))-exp((r-div)*dt/2))/(exp(sigma*sqrt(dt/2))-exp(-sigma*sqrt(dt/2))))^2;
pm = 1 - pu - pd;
discount = exp(- r * dt);            
SVals{N + 1} = S0 * u^N * d.^(0:1:2*N+1);
PVals{N + 1} = max(K - SVals{N + 1}, 0);
for i = N:-1:1
    svals = zeros(1, 2 * i - 1);
    pvals = zeros(1, 2 * i - 1);
    for j = 1:(2 * i - 1)
        svals(j) = S0 * u^(i - 1) * d^(j - 1);   
        Early = discount * (pu * PVals{i + 1}(j) + pm * PVals{i + 1}(j + 1) + pd * PVals{i + 1}(j + 2));
        Internal = max(K - svals(j), 0);
        pvals(j) = max(Early, Internal);
    end
    SVals{i} = svals;
    PVals{i} = pvals;
end
price = PVals{1};
% price = TrinomialTreeCRR(50,50,0.05,5/12,0.4,0,1000)
end

