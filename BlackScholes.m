function price = BlackScholes(S0, K, r, T, sigma, div)
% Black Scholes for European Put Option
d1 = (log(S0 / K) + (r - div + 0.5 * sigma^2) * T) / (sigma * sqrt(T));
d2 = d1 - sigma * sqrt(T);
price = K * exp(- r * T) * normcdf(- d2) - S0 * exp(- div * T) * normcdf(- d1);
end

