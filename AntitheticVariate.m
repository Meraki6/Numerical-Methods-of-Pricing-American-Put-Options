function price = AntitheticVariate(S0, K, r, T, sigma, div, N, M, ctr)
price = 0;
for i=1:ctr
    dt = T/N;
    nudt = (r - div - 0.5 * sigma^2) * dt;
    sidt = sigma * sqrt(dt);
    randm = rand(round(M/2),N);
    randu = norminv(randm);
    randanti = norminv(1 - randm);
    RandMat = [nudt + sidt * randu; nudt - sidt * randu];
    RandMatAnti = [nudt + sidt * randanti; nudt - sidt * randanti];
    priceu = LongstaffSchwartz(S0, K, r, T, N, M, RandMat);
    priceuanti = LongstaffSchwartz(S0, K, r, T, N, M, RandMatAnti);
    price = price + (priceu + priceuanti) / 2;
end
price = price / ctr;
% price = AntitheticVariate(50,50,0.05,5/12,0.4,0,5000,5000,500)
end