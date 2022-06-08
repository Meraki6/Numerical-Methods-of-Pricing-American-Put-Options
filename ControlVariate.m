function price = ControlVariate(S0, K, r, T, sigma, div, N, M, ctr)
EurBS = BlackScholes(S0, K, r, T, sigma, div);
EurMC = MCEuropean(S0, K, r, T, sigma, div, N, ctr);
AmeMC = 0;
for i=1:ctr
    dt = T/N;
    nudt = (r - div - 0.5 * sigma^2) * dt;
    sidt = sigma * sqrt(dt);
    randm = rand(round(M/2),N);
    randu = norminv(randm);
    RandMat = [nudt + sidt * randu; nudt - sidt * randu];
    AmeMC = AmeMC + LongstaffSchwartz(S0, K, r, T, N, M, RandMat);
end
AmeMC = AmeMC / ctr;
price = AmeMC + EurBS - EurMC;
% price = ControlVariate(50,50,0.05,5/12,0.4,0.02,5000,5000,500)
end

