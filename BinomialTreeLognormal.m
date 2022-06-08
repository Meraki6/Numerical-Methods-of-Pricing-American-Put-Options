function price = BinomialTreeLognormal(S0, K, r, T, sigma, div, N)
% Binomial Tree Model of lognormal
dt = T / N;
u = exp((r - div - 0.5 * sigma ^ 2) * dt + sigma * sqrt(dt));
d = exp((r - div - 0.5 * sigma ^ 2) * dt - sigma * sqrt(dt));
p = 0.5;
discount = exp(- r * dt);            
SVals{N + 1} = S0 * u^N * (d/u).^(0:1:N);
PVals{N + 1} = max(K - SVals{N + 1}, 0);
for i = N:-1:1
    svals = zeros(1, i);
    pvals = zeros(1, i);
    for j = 1:i
        svals(j) = S0 * u^(i - 1) * (d/u)^(j - 1);   
        Early = discount * (p * PVals{i + 1}(j) + (1 - p) * PVals{i + 1}(j + 1));
        Internal = max(K - svals(j), 0);
        pvals(j) = max(Early, Internal);
    end
    SVals{i} = svals;
    PVals{i} = pvals;
end
price = PVals{1};
end

