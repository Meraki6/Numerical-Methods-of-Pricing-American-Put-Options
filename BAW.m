StartDate = 'Jan-1-2016';
EndDate = 'Jan-1-2017';
Basis = 1;
Compounding = -1;
Rates = 0.03;
RateSpec = intenvset('ValuationDate',StartDate,'StartDate',StartDate,'EndDate',EndDate, ...
'Rates',Rates,'Basis',Basis,'Compounding',Compounding);

DividendType = {'cash'};
DividendAmounts = ones(1,2) * 0.4;
ExDividendDates = {'Jul-1-2016','Dec-31-2016'};
AssetPrice = 50;
Volatility = 0.5;
StockSpec = stockspec(Volatility, AssetPrice, DividendType, DividendAmounts, ExDividendDates);

OptSpec = 'put';
Strike = 50;
Settle = 'Jan-1-2016';
Maturity = EndDate;
Price = optstockbybaw(RateSpec,StockSpec,Settle,Maturity,OptSpec,Strike)




