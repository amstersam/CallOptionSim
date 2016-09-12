%% Validation of the volatility for the European Option Call
clear all
clc;
tic;
rng(123);
load marketdata


I=20000; % Basic input 10000 - change to 20000 for ass. 2

dt=1/250;

% Interest rate
r=0.007;

% Dividend yield
q=0.0093;
%volatility
volatility=0.3410;
%drift
mu=r-q;


standdevpayoff=zeros(size(KK,1),1);
ExpPayoff=zeros(size(KK,1),1);
SE=zeros(size(KK,1),1);
lb=zeros(size(KK,1),1);
ub=zeros(size(KK,1),1);

numValidations=6; % select which option is to be simulated (1...6)

 
    %% Initializations
    
    K=KK(numValidations,1);
    stDate=datestr(SettleDate(numValidations,1));
    enDate=datestr(MaturityDate(numValidations,1));
% Maturity time    
    M=daysdif(stDate,enDate,13); % calculates the business days
    
    T=M/250;
    

sample=200; % sample size for CLT
MatrixForCLT=zeros(sample,1);
%for clt=1:sample  %% this loop was specifically used for the CLT 
% the computation time for the CLT verification is appr. 5h 
    S=zeros(I,M);
    S(1:I,1)=72.17;
    Cmarket=CallMarket(numValidations,1);
    
    for i=1:I
        for j=1:M
            S(i,1 + j) = S(i,j) * (1 + mu * dt + volatility * sqrt(dt) * random('Normal',0,1));
        end
    end


    Payoff=zeros(I,1);
    for k=1:I
        if S(k,M+1)>K
            Payoff(k,1)=(S(k,M+1)-K)*exp(-r*T);
        end
    end
    standdevpayoff(numValidations,1)=std(Payoff);
    ExpPayoff(numValidations,1)=mean(Payoff)
    sampleVar(numValidations,1)=var(Payoff);
    SE(numValidations,1)=sqrt(sampleVar(numValidations,1))/(I^0.5);
    

    
    lb(numValidations,1)=ExpPayoff(numValidations,1)-1.96*SE(numValidations,1);
    ub(numValidations,1)=ExpPayoff(numValidations,1)+1.96*SE(numValidations,1);
    CI=[lb ub];
    if (Cmarket<=ub(numValidations,1)) && (Cmarket>=lb(numValidations,1))
        disp('OK!');
    else 
        disp('Not OK!');
    end
    
%    MatrixForCLT(clt,1)=ExpPayoff(numValidations,1);
%    RE=SE/ExpPayoff;
    toc;
%end
%% Code for SSLN
%% Code specifically to create matrix for diagrams in excel 
% toExcel=zeros(I,1);
% toExcelready=zeros(I,1);
% toExcel(1,1)=Payoff(1,1);
% for jjj=2:I
%     
%     toExcel(jjj,1)=toExcel(jjj-1,1)+Payoff(jjj,1); 
%     
% end
% for kkk=1:I
%     toExcelready(kkk,1)=toExcel(kkk,1)/kkk;
% end
