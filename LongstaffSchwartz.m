function price = LongstaffSchwartz(S0,K,r,T,N,M,RandMat)
dt = T/N;
% discount rates over different time intervals
discountVet = exp(-r*dt*(1:N)');
a = zeros(3,1); % regression parameters
% generate sample paths
LogPaths = cumsum([log(S0)*ones(M,1) , RandMat] , 2);
SPaths = exp(LogPaths);
SPaths(:,1) = []; % get rid of starting prices
CashFlows = max(0, K - SPaths(:,N));
% first set exercise time at expiration for convenience
ExerciseTime = N*ones(M,1);
for step = N-1:-1:1
    InMoney = find(SPaths(:,step) < K);
    XData = SPaths(InMoney,step);
    RegrMat = [ones(length(XData),1), XData, XData.^2];
    YData = CashFlows(InMoney) .* discountVet(ExerciseTime(InMoney) - step);
    a = RegrMat \ YData;
    IntrinsicValue = K - XData;
    ContinuationValue = RegrMat * a;
    Exercise = find(IntrinsicValue > ContinuationValue);
    k = InMoney(Exercise);
    CashFlows(k) = IntrinsicValue(Exercise);
    ExerciseTime(k) = step;
end
price = mean(CashFlows.*discountVet(ExerciseTime));
end