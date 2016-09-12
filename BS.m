function C = BS(S, K, volatility, r,q, time)

%% code to calculate the price of call option with BS formula
d1 = (log(S / K) + (r -q + volatility^2 / 2) * time) / volatility / sqrt(time);
d2=d1-volatility*sqrt(time);
C=S*normcdf(d1)*exp(-q*time) - K*normcdf(d2)*exp(-r*time);