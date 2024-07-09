%% risk-aversion case
a=0.3;
b=[0.001 0.005 0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35];
c=1;
d=10;
e=0.2;
f=4;
g=0.5;
X=size(b,2);
energy_price=zeros(24,X);
reserve_price=zeros(24,X);
storage_price=zeros(25,X);
SOC=zeros(25,X);
cost=zeros(1,X);
rho_g=zeros(76,24,X);
rho_e=zeros(24,X);
Pdis=zeros(24,X);
Pch=zeros(24,X);
Pg=zeros(24,X);

% Dch_value=zeros(24,365);
% Ddis_value=zeros(24,365);
% SOC_value1_sorted=zeros(23,365);
% storage_price_sorted=zeros(23,365);

parfor i=1:X
[cost(:,i),energy_price(:,i),storage_price(:,i),reserve_price(:,i),SOC(:,i),~]=MultiPeriod_ED_CC(a,b(i),c,d,e,f,g);
end
energy_price_average=mean(energy_price,1);
reserve_price_total=sum(reserve_price,1);
storage_price_average=mean(storage_price(2:end,:),1);
