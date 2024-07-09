
a=0.3;
b=0.05;
c=0.5;
d=10;
e=0.2;
f=4;
g=0.1:0.1:0.9;
X=size(g,2);
energy_price=zeros(24,X);
reserve_price=zeros(24,X);
storage_price=zeros(25,X);
SOC=zeros(25,X);
cost=zeros(1,X);
Dch_value=zeros(24,X);
Ddis_value=zeros(24,X);

parfor i=1:X
[cost(:,i),energy_price(:,i),storage_price(:,i),reserve_price(:,i),SOC(:,i),~]=MultiPeriod_ED_CC(a,b,c,d,e,f,g(i));
end

% T=24;
% load data.mat
% T=24; %time duration
% %% constants
% wind_ratio=a;% wind ratio
% netload_error=load_error-max(max(load_real))*(wind_real-wind_forecast)*wind_ratio;
% 
% netload_error=reshape(netload_error',96,365)';% netload error, [-, +]
% 
% % 初始化一个新的 365x24 的数组 B
% B = zeros(365, 24);
% 
% for i = 1:365
%     
%     reshapedRow = reshape(netload_error(i, :), 4, T);
%     B(i, :) = mean(reshapedRow, 1);
% end
% 
% netload_error=B;
% 
% epsilon=b; % probability level
% for i=1:T
% % [mu(i),sigma(i)] = normfit(netload_error(:,i),0.05); % μ, σ time different
% 
% mu(i)=mean(netload_error(:,i));
% sigma(i)=std(netload_error(:,i));
% 
% end
% 
% wind_forecast=reshape(wind_forecast',96,365)';
% netload=load_forecast-max(max(load_real))*wind_forecast*wind_ratio; % netload forcast
% C = zeros(365, T);
% 
% 
% for i = 1:365
%     reshapedRow = reshape(netload(i, :), 4, T);
%     C(i, :) = mean(reshapedRow, 1);
% end
% 
% netload=C;
% 
% 
% sigma=sigma*c;
% 
% 
% load_c=netload(d,:);
% 
% load_d1=zeros(1,T);
% load_d2=zeros(1,T);
% 
% 
% for i=1:T
% 
% 
% load_d1(:,i)= quantile(B(:,i), epsilon);
% 
% 
% end
% 
% for i=1:T
% 
% 
% load_d2(:,i)= quantile(B(:,i), 1-epsilon);
% 
% end
% 
% C1=0.95*0.95*load_d1./load_d2;
% C2=0.95*(load_d2/(0.95*0.95)-load_d1)./load_d2;
% C3=0.95./load_d2;
% C4=-20*mu*0.95./load_d2;
% 
% B1=load_d2./(load_d1*0.95*0.95);
% B2=(load_d1*0.95*0.95-load_d2)./(0.95*load_d1);
% B3=1./(load_d1*0.95);
% B4=20*(load_d2-0.95*0.95*load_d1-mu)./(load_d1*0.95);
% 
% charge_price=C1.*storage_price(2:end)'+C2.*energy_price'-C3.*reserve_price'+C4;
% discharge_price=B1.*storage_price(2:end)'+B2.*energy_price'+B3.*reserve_price'+B4;
% 
% 
