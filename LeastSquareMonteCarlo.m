function Price = LeastSquareMonteCarlo(S0, K, r, T, sigma, N, M, ran)
dt = T/N;
t = 0:dt:T;
t = repmat(t',1,M);
R = exp((r-sigma^2/2)*dt + sigma*sqrt(dt)*ran);
S = cumprod([S0*ones(1,M); R]);
ExTime = (M+1)*ones(N,1); 
CF = zeros(size(S)); % Cash flow matrix
CF(end,:) = max(K-S(end,:),0); % Option only pays off if it is in the money
for ii = size(S)-1:-1:2
    Idx = find(S(ii,:) < K); % Find paths that are in the money at time ii
    X = S(ii,Idx)'; 
    X1 = X/S0;
    Y = CF(ii+1,Idx)' * exp(-r*dt); % Discounted cashflow
    R = [ ones(size(X1)) (1-X1) 1/2*(2-4*X1-X1.^2)];
    a = R\Y; % Linear regression step
    C = R*a; % Cash flows as predicted by the model
    Jdx = max(K-X,0) > C; % Immediate exercise better than predicted cashflow
    nIdx = setdiff((1:M),Idx(Jdx));
    CF(ii,Idx(Jdx)) = max(K-X(Jdx),0);
    ExTime(Idx(Jdx)) = ii;
    CF(ii,nIdx) = exp(-r*dt)*CF(ii+1,nIdx);
end
Price = mean(CF(2,:))*exp(-r*dt);
end

