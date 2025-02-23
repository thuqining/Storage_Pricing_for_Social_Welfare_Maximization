function [system_cost1_av,system_cost2_av,storage_profit1_av,storage_profit2_av,consumer_payment1_av,consumer_payment2_av,price_std1,price_std2]=Cost_compare(sigma,c)
%% parameter
a=0.3; % wind ratio
b=0.05;% confidence level
%c=3; %sigma, original forecast error is too small
d=10; %day number
e=0.2; %storage capacity ratio (load capacity)
f=4; %storage duration
g=0.5; %initial SoC
h=10; %storage marginal cost

%sigma=50; % price sigma
scenarios=100; %scenarios number for MCS

%% DA price generation
 [energy_price]=DA_price(a,b,c,d,e,f,g,h);


%% calculate value function
Ts = 1; % time step
lambda=energy_price;
T = numel(lambda); % number of time steps
Dur = f; 
Pr = 1/(Dur); % normalized power rating wrt energy rating
P = Pr*Ts; % actual power rating taking time step size into account
eta = 0.95; % efficiency
mc = h; % marginal discharge cost - degradation
ed = 1/(1500); % SoC sample granularity, 1500 segments/day
ef = 0; % final SoC target level, use 0 if none
Ne = floor(1/ed)+1; % number of SOC samples
seg_num=5;
[VF,~] = value_fcn_calcu(lambda,seg_num,Ne, T, mc, P, eta, ed, ef, sigma);
% VF_av=mean(VF.k);

%% strategic bid for storage-1*24
storage_bid_d1=VF.k/eta+mc;
storage_bid_c1=VF.k*eta;


%SoC-dependent bid bounds
g0=0:1/(seg_num-1):1;
parfor i=1:seg_num
[bid_bound_c(i,1),bid_bound_d(i,1)]=MultiPeriod_ED_CC(a,b,c,d,e,f,g0(i),h);
end

%% Capped storage bid
storage_bid_d2=min(storage_bid_d1,bid_bound_d);

%storage_bid_c2=min(storage_bid_c1,bid_bound_c);
storage_bid_c2=storage_bid_c1;%no charge


%% single-period ED
% g1=g/seg_num*ones(scenarios,seg_num);
% g2=g/seg_num*ones(scenarios,seg_num);
g1=repmat([0.2 0.2 0.1 0 0],scenarios,1);
g2=repmat([0.2 0.2 0.1 0 0],scenarios,1);
system_cost1_av=zeros(scenarios,T);
storage_profit1_av=zeros(scenarios,T);
consumer_payment1_av=zeros(scenarios,T);
system_cost2_av=zeros(scenarios,T);
storage_profit2_av=zeros(scenarios,T);
consumer_payment2_av=zeros(scenarios,T);
price_1=zeros(scenarios,T);
price_2=zeros(scenarios,T);
% SOC_value1=zeros(scenarios,seg_num);
% SOC_value2=zeros(scenarios,seg_num);

for i=1:T
[system_cost1_av(:,i),storage_profit1_av(:,i),consumer_payment1_av(:,i),system_cost2_av(:,i),storage_profit2_av(:,i),consumer_payment2_av(:,i),price_1(:,i),price_2(:,i),SOC_value1,SOC_value2]=RT_price_single(a,scenarios,c,d,e,f,g1,g2,h,storage_bid_d1(:,i),storage_bid_c1(:,i),storage_bid_d2(:,i),storage_bid_c2(:,i),i);
g1=SOC_value1;
g2=SOC_value2;
end

system_cost1_av=mean(sum(system_cost1_av,2));
system_cost2_av=mean(sum(system_cost2_av,2));

storage_profit1_av=mean(sum(storage_profit1_av,2));
storage_profit2_av=mean(sum(storage_profit2_av,2));

consumer_payment1_av=mean(sum(consumer_payment1_av,2));
consumer_payment2_av=mean(sum(consumer_payment2_av,2));

price_std1=max(std(price_1,1));
price_std2=max(std(price_2,1));










