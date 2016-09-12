%% Calculates the call by the Heston Model
clear all
clc;
load marketdata
rng(123);
tic;

%% Decide which case of data to run (numSim=1,...6)
numSim=1;
K=KK(numSim,1);
stDate=datestr(SettleDate(numSim,1));
enDate=datestr(MaturityDate(numSim,1));
M=daysdif(stDate,enDate,13);
T=M/250;

S0 = 72.17; 
V0 = 0.3410^2; 
r = 0.007;
q=0.0093;
mu=r-q;
k = 2;
theta =0.0625; 
ksi=0.2;


dt = 1/250; 
rho =-0.6;
numberOfSimulations = 10000;
VV=zeros(numberOfSimulations,1);
simPathSum = 0;
vPath=zeros(numberOfSimulations,1);
RootvPath=zeros(numberOfSimulations,1);
for i=1:numberOfSimulations
    
        v1 = V0;
        s1 = S0;

        for j=1:M
            z1 = randn;
            z2 = randn;
            
            
            s1 = s1 * exp((mu - 0.5 * v1) * dt + sqrt(v1) * z1 * sqrt(dt));
            v1 = v1 + k*(theta - v1)*dt + ksi*sqrt(v1)* (rho*z1 + sqrt(1- rho^2)*z2)*sqrt(dt);
            if v1<0
                disp('negative volatility');
                %v1=-v1;
                v1=0;
               
            else
                 
                %disp('ok!');
            end
        
        end
        vPath(i,1)=v1;
        RootvPath(i,1)=sqrt(v1);
        simPath=s1;
    
        simPathSum = simPathSum + exp((-T) * r) * max(simPath - K, 0);
   
end
toc;
HestonCall=simPathSum/numberOfSimulations

%% Adds data to plot in matrix and calculates convergence ready to be used by script plotresult.m
toExcel=zeros(numberOfSimulations,1);
toExcelready=zeros(numberOfSimulations,1);
toExcel(1,1)=RootvPath(1,1);
for jjj=2:numberOfSimulations
    
    toExcel(jjj,1)=toExcel(jjj-1,1)+RootvPath(jjj,1); 
    
end

 for kkk=1:numberOfSimulations
     toExcelready(kkk,1)=toExcel(kkk,1)/kkk;
 end