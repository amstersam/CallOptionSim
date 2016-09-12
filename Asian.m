%% Calculation of the Asian Option
clear all
clc;
tic;
rng(123);
load marketdata


I=10000; %% 100000 for SLLN

dt=1/250;

% Interest rate
r=0.007;

% Dividend yield
q=0.0093;
%volatility
volatility=0.3410;
mu=r-q;
numValidations=1; %% picks the number of the option from the market data input file


    %% Initializations
    
    K=KK(numValidations,1);
    stDate=datestr(SettleDate(numValidations,1));
    enDate=datestr(MaturityDate(numValidations,1));
    
    M=daysdif(stDate,enDate,13);
    T=M/250;
%     S=zeros(I,M);
%     S(1:I,1)=72.17;
S0=72.17;
    Cmarket=CallMarket(numValidations,1);

S=zeros(I,M+1);
arithmetic_mean=zeros(I,1);
Payoff=zeros(I,1);


%MCmean_simulation=zeros(I,1);
sample=10;  % for CLT
MatrixForCLT=zeros(sample,1);

PayoffGeometric=zeros(I,1);
% for clt=1:sample %%

for j = 1:I %Simulations SLLN
    mTemp = (r - q - volatility^2/2)*dt;
    sTemp = volatility*sqrt(dt);
    Y= mTemp+sTemp*randn(1,M); 
    S(j,:) = cumsum([log(S0), Y],2);
    arithmetic_mean(j,1) = mean(exp(S(j,:))); 
    geometric_mean=exp(mean(S(j,:)));
    
    Payoff(j,1)= max([arithmetic_mean(j,1)-K;zeros(1)]);
    PayoffGeometric(j,1)= max([geometric_mean-K;zeros(1)]);

end
    AsianMC = mean(Payoff) 
    Geometric_asian=mean(PayoffGeometric);
     MCstd = std(Payoff)/sqrt(I);
%MatrixForCLT(clt,1)=MCmean;
%end
toc;
    muGeom = 1/2*(r-volatility^2/2)*(1+1/M);
    volGeom = sqrt((volatility^2)/3*(1+1/M)*(1+1/(2*M)));
    S0Geom = S0*exp(T*((volGeom^2)/2+muGeom-r));
    CallBS = BS(S0Geom,K,volGeom,r,q,T);
    C = cov(AsianMC,Geometric_asian);
    b = C(1,2)/C(1,1);
R = exp(-r*T);    
    CV = (R*AsianMC) -b*((R*Geometric_asian) - CallBS);
    CVmean = mean(CV);
    CVstd = std(CV)/sqrt(I);

lb = AsianMC -1.96* MCstd/2;
ub = AsianMC +1.96* MCstd/2;
lbCV = CVmean -1.96* CVstd/2;
ubCV = CVmean +1.96* CVstd/2;
Price_AM = AsianMC;
CI_AM = [lb ub];
Price = CVmean;
CI = [lbCV ubCV];
%% Not necessary for CLT
toExcel=zeros(I,1);
toExcelready=zeros(I,1);
toExcel(1,1)=Payoff(1,1);
for jjj=2:I
    
    toExcel(jjj,1)=toExcel(jjj-1,1)+Payoff(jjj,1);
    
end
for kkk=1:I
    toExcelready(kkk,1)=toExcel(kkk,1)/kkk;
end
