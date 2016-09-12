clear all
clc;
rng(123);
load marketdata
tic;

%for numSim=1:6
numSim=1;
stDate=datestr(SettleDate(numSim,1));
enDate=datestr(MaturityDate(numSim,1));

% SettleDate='12-Feb-2016';
% MaturityDate='17-Jun-2016';
% NumDays=daysdif(SettleDate,MaturityDate,13);
M=daysdif(stDate,enDate,13);
% Stock price
S=72.17;

% Strike price
%K=75;
K=KK(numSim,1);
% Interest rate
r=0.007;

% Dividend yield
q=0.0093;
t=0; % now
% Maturity time
time=M/250;
%Cmarket=4.55;
Cmarket=CallMarket(numSim,1);
% Various volatilities
volatility=(0.05:0.0001:0.5)';
C=zeros(size(volatility,1),1);
%C=funcRoot(volatility, S, K, r, time);
for i=1:size(volatility,1)
    d(i,1) = (log(S / K) + (r -q + volatility(i) * volatility(i) / 2) * time) / volatility(i) / sqrt(time);
    d(i,2)=d(i,1)-volatility(i)*sqrt(time);
end
for i=1:size(volatility,1)
    C(i,1)=S*normcdf(d(i,1))*exp(-q*time) - K*normcdf(d(i,2))*exp(-r*time);
    
end

mm=zeros(size(C,1),1);
for j=1:size(volatility,1)
    mm(j,1)=(abs(C(j,1)-Cmarket));
end
[mmm, ind]=min(mm);
disp(volatility(ind));
vol(numSim,1)=volatility(ind);
toc;
%end % end loop in case someone wants to simulate all input data
%% confirmartion with Newton-Rahpson method
%vol2=volatilityFast(Cmarket, S, K, r, t, time, q) 
mean(vol)
%% confirmartion with Matlab function from financial toolbox
%volatilityAEM = blsimpv(S, K, r, time, Cmarket, 10, q);
toc;


            
            