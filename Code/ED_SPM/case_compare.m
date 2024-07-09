% compare two pricing mechanism

%function result=case_compare(a,b,c,d,e,f,g)
a=0.9;
b=0.05;
c=1;
d=10;
e=0.2;
f=4;
g=0.5;
%% proposed opportunity pricing framework
[cost,energy_price,storage_price,reserve_price,SOC,Pdis,Pch,Pg,electricitycost]=MultiPeriod_ED_CC(a,b,c,d,e,f,g);

storagerevenue=(Pdis/0.95-Pch*0.95)*storage_price(2:end);
storagecost=sum(10*Pdis);
storageprofit=storagerevenue-storagecost;
generationcost=cost-storagecost;
generationrevenue=Pg*energy_price;
generationprofit=generationrevenue-generationcost;

%% generate price uncertainty
energy_price1=zeros(24,100);
ESenergycap=zeros(24,100);
parfor i=1:100
[energy_price1(:,i),ESenergycap(:,i)]=MultiPeriod_ED_scenario(a,b,c,d,e,f,g);
end
sigma=std(energy_price1');

%% calculate value function
Ts = 1; % time step

lambdaH=energy_price;

T = numel(lambdaH); % number of time steps

Dur = f; 
Pr = 1/(Dur); % normalized power rating wrt energy rating
P = Pr*Ts; % actual power rating taking time step size into account
eta = 0.95; % efficiency
mc = 10; % marginal discharge cost - degradation
ed = 1/(1500); % SoC sample granularity, 1500 segments/day
% ed = 1/20000; % SoC sample granularity for 1 month duration
ef = 0; % final SoC target level, use 0 if none
Ne = floor(1/ed)+1; % number of SOC samples
e0 = g;

vEnd = zeros(Ne,1);  % generate value function samples

vEnd(1:floor(ef*Ne)) = 1e2; % use 100 as the penalty for final discharge level


%%
tic
v = zeros(Ne, T+1); % initialize the value function series
% v(1,1) is the marginal value of 0% SoC at the beginning of day 1
% V(Ne, T) is the maringal value of 100% SoC at the beginning of the last operating day
v(:,end) = vEnd; % update final value function

% process index
es = (0:ed:1)';
Ne = numel(es);
% calculate soc after charge vC = (v_t(e+P*eta))
eC = es + P*eta; 
% round to the nearest sample 
iC = ceil(eC/ed)+1;
iC(iC > (Ne+1)) = Ne + 2;
iC(iC < 2) = 1;
% calculate soc after discharge vC = (v_t(e-P/eta))
eD = es - P/eta; 
% round to the nearest sample 
iD = floor(eD/ed)+1;
iD(iD > (Ne+1)) = Ne + 2;
iD(iD < 2) = 1;

for t = T:-1:1 % start from the last day and move backwards
    vi = v(:,t+1); % input value function from tomorrow
   % vo = CalcValueNoUnc(lambdaH(t), c, P, eta, vi, ed, iC, iD);
    vo = CalcValueNormal(lambdaH(t),sigma(t), mc, P, eta, vi, ed, iC, iD);
    v(:,t) = vo; % record the result 
end

tElasped = toc;

%% convert value function to 5 segments
vAvg = zeros(101,T+1);

NN = (Ne-1)/100;

vAvg(1,:) = v(1,:);

for i = 1:100
   vAvg(i+1,:) = mean(v((i-1)*NN + 1 + (1:NN),:)); 
end

%% perform the actual arbitrage
eS = zeros(T,1); % generate the SoC series
pS = eS; % generate the power series

SOC1 = g; % initial SoC
for t = 1:T % start from the first day and move forwards
    vv = v(:,t+1); % read the SoC value for this day
   [SOC1 , p] =  Arb_Value(lambdaH(t), vv, SOC1 , P, 1, eta, mc, size(v,1));
   eS(t) = SOC1; % record SoC
   pS(t) = p; % record Power
end

storageprofit1 = (sum(pS.*lambdaH) - sum(pS(pS>0)*mc))*ESenergycap(1);
storagerevenue1 = sum(pS.*lambdaH)*ESenergycap(1);
storagecost1=sum(pS(pS>0)*mc)*ESenergycap(1);
fprintf('Profit=%e, revenue=%e',storageprofit1, storagerevenue1)
solTimeOut = toc;

%% output
[generationcost1,electricitycost1]=MultiPeriod_ED_ESbid(a,b,c,d,e,f,pS);
cost1=generationcost1+storagecost1;
result=[storageprofit generationcost cost electricitycost;storageprofit1 generationcost1  cost1 electricitycost1];


