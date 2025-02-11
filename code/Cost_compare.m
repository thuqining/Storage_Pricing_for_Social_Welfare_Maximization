function [system_cost1_av,system_cost2_av,storage_profit1_av,storage_profit2_av,consumer_payment1_av,consumer_payment2_av]=Cost_compare(sigma,c)
%% parameter
a=0.3; % wind ratio
b=0.05;% confidence level
%c=1; %sigma, original forecast error is too small
d=10; %day number
e=0.2; %storage capacity ratio (load capacity)
f=4; %storage duration
g=1; %initial SoC
h=10; %storage marginal cost

%sigma=40; % price sigma
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
%% optimal bid for storage
storage_bid_d1=VF.k/eta+mc;
storage_bid_c1=VF.k*eta;

% [bid_bound_c,bid_bound_d]=MultiPeriod_ED_CC(a,b,c,d,e,f,0,h,l);

%Default bid and bid bounds
g0=0.2:0.2:1;
parfor i=1:seg_num
[bid_bound_c(i,1),bid_bound_d(i,1)]=MultiPeriod_ED_CC(a,b,c,d,e,f,g0(i),h);
end

%% Capped storage bid
storage_bid_d2=min(storage_bid_d1,bid_bound_d);
% storage_bid_c2=max(storage_bid_c1,bid_bound_c);
storage_bid_c2=min(storage_bid_c1,bid_bound_c);
%storage_bid_c2=storage_bid_c1;

[system_cost1_av,storage_profit1_av,consumer_payment1_av,system_cost2_av,storage_profit2_av,consumer_payment2_av]=RT_price(a,scenarios,c,d,e,f,g,h,storage_bid_d1,storage_bid_c1,storage_bid_d2,storage_bid_c2);





