%@Create Time    : 2024/3/15
%@Author  : Ning Qi nq2176@columbia.edu
%@File    : MultiPeriod_ED_CC.m
%@Description : This is the main code for multi-period chance-constrained
%economic dispatch for system operator and can generate both dispatch and
%price for storage and generator
% storage: energy and reserve market

function [bid_bound_c,bid_bound_d]=MultiPeriod_ED_DT(a,b,c,d,e,f,g,h)
%% inputs
% a=0.3;
% b=0.05;
% c=1; %sigma, original forecast error is too small
% d=10; %day number
% e=0.2; %storage capacity ratio (load capacity)
% f=4; %storage duration
% g=0; %initial SoC
% h=10; %storage marginal cost
% 
MS=h;

load data.mat
%% case setting
T=24; %time duration
wind_ratio=a;% wind ratio
netload_error=load_error-max(max(load_real))*(wind_real-wind_forecast)*wind_ratio;
netload_error=reshape(netload_error',96,365)';% netload error, [-, +]
B = zeros(365, 24);
for i = 1:365
    reshapedRow = reshape(netload_error(i, :), 4, T);
    B(i, :) = mean(reshapedRow, 1);
end
netload_error=B;

epsilon=b; % probability level

for i=1:T
mu(i)=mean(netload_error(:,i));
sigma(i)=std(netload_error(:,i));
end

wind_forecast=reshape(wind_forecast',96,365)';
netload=load_forecast-max(max(load_real))*wind_forecast*wind_ratio; % netload forcast
C = zeros(365, T);
for i = 1:365
    reshapedRow = reshape(netload(i, :), 4, T);
    C(i, :) = mean(reshapedRow, 1);
end
netload=C;

sigma=sigma*c;% change the STD of forecast error for net load

load_c=netload(d,:);



load_d2=normrnd(mu, sigma);


if c==0
    load_d2=load_c*0;
end




num_e=1; %number of storage
num_g=size(data_generator,1); %number of generator
yita=0.95; % storage efficiency
% num_e=2; %number of storage
% yita=[0.8  0.95 ]'; % storage efficiency
SOC0=g; % initial SoC
ESpowercap=e*max(max(netload))/num_e; % storage power cap
ESenergycap=ESpowercap*f; % storage energy cap
G_min=data_generator(:,2)*0;% minimum capacity=0
G_max=data_generator(:,1)*1;
RU=data_generator(:,6);
RD=data_generator(:,6);
C_0=data_generator(:,8);
C_1=data_generator(:,9);
C_2=data_generator(:,10);

%% decision variables
Pch=sdpvar(num_e,T); % storage charge
Pdis=sdpvar(num_e,T); % storage discharge
Dch=binvar(num_e,T); % storage charge state
Ddis=binvar(num_e,T); % storage discharge state
SOC=sdpvar(num_e,T+1); % storage SOC
Pg=sdpvar(num_g,T); % generator power
rg=sdpvar(num_g,T); % generator power

%% constraints
cons=[];
% power balance-energy price dual
for j=1:T    
cons=[cons, sum(Pg(:,j),1)+sum(Pdis(:,j),1)-sum(Pch(:,j),1)>=load_c(:,j)+load_d2(:,j)];
end
% SoC recovery

cons=[cons,SOC(:,1)==SOC0];
for i=1:num_e
for j=2:T+1    
cons=[cons, SOC(i,j)==SOC(i,j-1)+((Pch(i,j-1)).*yita(i)+(-Pdis(i,j-1))./yita(i))/ESenergycap];
end
end


% storage power bounds
for j=1:T    
cons=[cons,Pch(:,j)>=0];
end

for j=1:T    
cons=[cons,Pdis(:,j)>=0];
end


for j=1:T    
cons=[cons,Pch(:,j)<=Dch(:,j)*ESpowercap];
end
% 
for j=1:T    
cons=[cons,Pdis(:,j)<=Ddis(:,j)*ESpowercap];
end


% complementary charging-discharging
for j=1:T    
cons=[cons,0<=Dch(:,j)+Ddis(:,j)<=1];
end

% generator power bounds
for j=1:T    
cons=[cons,G_min<=Pg(:,j)<=G_max-rg(:,j)];
end


% generator reserve
for j=1:T    
cons=[cons,0<=rg(:,j)];
end


for j=1:T    
cons=[cons,sum(rg(:,j),1)>=0.1*load_c(:,j)];
end



% generator RAMPUP
for j=2:T    
cons=[cons,Pg(:,j)-Pg(:,j-1)<=RU];
end

% generator RAMPDOWN
for j=2:T    
cons=[cons,-RD<=Pg(:,j)-Pg(:,j-1)];
end

% storage SoC bounds
for j=2:T+1    
cons=[cons, 0<=SOC(:,j)<=1];
end


%% objective
obj=0;

% % quadratic+marginal cost
for i=1:num_g
obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:))+C_2(i,1)*(Pg(i,:).*Pg(i,:));
end
obj=sum(obj,2);
for j=1:T
obj=obj+sum((Pdis(:,j))*MS,1);
end

options = sdpsettings('verbose', 1, 'solver', 'gurobi');
results= optimize(cons,obj,options);

Dch_value=value(Dch);
Ddis_value=value(Ddis);

%% using integer value to solve it again so that we can derive dual

cons=[];
% power balance-energy price dual
for j=1:T    
cons=[cons, sum(Pg(:,j),1)+sum(Pdis(:,j),1)-sum(Pch(:,j),1)>=load_c(:,j)+load_d2(:,j)];
end
% SoC recovery

cons=[cons,SOC(:,1)==SOC0];
for i=1:num_e
for j=2:T+1    
cons=[cons, SOC(i,j)==SOC(i,j-1)+((Pch(i,j-1)).*yita(i)+(-Pdis(i,j-1))./yita(i))/ESenergycap];
end
end


% storage power bounds
for j=1:T    
cons=[cons,Pch(:,j)>=0];
end

for j=1:T    
cons=[cons,Pdis(:,j)>=0];
end


for j=1:T    
cons=[cons,Pch(:,j)<=Dch_value(:,j)*ESpowercap];
end
% 
for j=1:T    
cons=[cons,Pdis(:,j)<=Ddis_value(:,j)*ESpowercap];
end


% generator power bounds
for j=1:T    
cons=[cons,G_min<=Pg(:,j)<=G_max-rg(:,j)];
end


% generator reserve
for j=1:T    
cons=[cons,0<=rg(:,j)];
end


for j=1:T    
cons=[cons,sum(rg(:,j),1)>=0.1*load_c(:,j)];
end



% generator RAMPUP
for j=2:T    
cons=[cons,Pg(:,j)-Pg(:,j-1)<=RU];
end

% generator RAMPDOWN
for j=2:T    
cons=[cons,-RD<=Pg(:,j)-Pg(:,j-1)];
end

% storage SoC bounds
for j=2:T+1    
cons=[cons, 0<=SOC(:,j)<=1];
end


%% objective
obj=0;

% % quadratic+marginal cost
for i=1:num_g
obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:))+C_2(i,1)*(Pg(i,:).*Pg(i,:));
end
obj=sum(obj,2);
for j=1:T
obj=obj+sum((Pdis(:,j))*MS,1);
end

options = sdpsettings('verbose', 1, 'solver', 'gurobi');
results= optimize(cons,obj,options);


Pdis_value=value(Pdis);
Pch_value=value(Pch);
SOC_value=value(SOC)';
Pg_value=value(Pg);
Pg_value=sum(Pg_value,1);
obj_value=value(obj);

energy_price=abs(dual(cons(1:T)));
storage_price=abs(dual(cons(T+2:2*T+1))/ESenergycap);



Dch_value=zeros(num_e,T);
Ddis_value=zeros(num_e,T);


for i=1:num_e
    for j=1:T
 if  Pdis_value(i,j)>=0.1
    Ddis_value(i,j)=1;
 end

  if  Pch_value(i,j)>=0.1
    Dch_value(i,j)=1;
  end
    end
end

if ~isnan(find(Ddis_value>0))
    bid_bound_d=max(storage_price(Ddis_value>0))/yita+MS;
    LMP_max=max(energy_price(Ddis_value>0));
    oc_max=max(storage_price(Ddis_value>0));
else
    bid_bound_d=max(energy_price);%uncleared
end


if ~isnan(find(Dch_value>0))
    bid_bound_c=min(storage_price(Dch_value>0))*yita;
    LMP_min=min(energy_price(Dch_value>0));
else
    bid_bound_c=min(energy_price);%uncleared
end



storage_bid_c=storage_price'*yita;
storage_bid_d=storage_price'/yita+MS;






% storage_bid_c1=storage_price1(2:end)*yita-MS;
% storage_bid_d1=storage_price1(2:end)/yita+MS;

%storage_profit=sum(storage_price(2:end).*(Pdis_value/yita-Pch_value*yita))-sum((Pch_value-rho_e1_value.*mu)*MS+(Pdis_value+rho_e_value.*mu)*MS)+sum((rho_e_value+rho_e1_value)*reserve_price);
%storage_profit=sum((Pdis_value/yita-Pch_value*yita)*energy_price)-sum((Pch_value-rho_e1_value.*mu)*MS+(Pdis_value+rho_e_value.*mu)*MS)+sum((rho_e_value+rho_e1_value)*reserve_price);
% 
% %storage_profit=sum(storage_price(2:end).*(Pdis_value/yita-Pch_value*yita))-sum((Pdis_value+rho_e_value.*mu)*MS)+sum((rho_e_value+rho_e1_value)*reserve_price);
% generator_profit=sum(Pg_value*energy_price)+sum(rho_g_value*reserve_price)-(obj_value-sum((Pch_value+rho_e_value.*mu)*MS+(Pdis_value-rho_e1_value.*mu)*MS));


% figure(1)
% set(gcf,'unit','centimeters','position',[0,0,8,4])
% 
% % 画出现有曲线
% plot(storage_bid_c,'g--','LineWidth',1) % Charge Bids 使用绿色虚线
% hold on
% plot(storage_bid_d,'b--','LineWidth',1) % Discharge Bids 使用蓝色虚线
% hold on
% plot(energy_price,'k.-','LineWidth',1)
% hold on
% 
% % 添加 SoC 曲线
% yyaxis right
% plot(0:24, SOC_value, 'r-.', 'LineWidth', 1) % SoC 使用红色点划线
% hold on
% 
% % 设置 x 轴刻度和标签
% set(gca,'xtick',[0:4:T])
% xlim([0 T])
% set(gca,'xticklabel',{'0','4','8','12','16','20','24'})
% set(gca,'FontName','Times New Roman','FontSize',8)
% 
% % 设置左侧 y 轴用于价格
% yyaxis left
% ylabel('\fontsize{8}\fontname{Times new roman}Price ($/MWh)')
% set(gca,'YColor','k') % 左侧 y 轴刻度颜色为黑色
% 
% % 设置右侧 y 轴用于 SoC
% yyaxis right
% ylabel('\fontsize{8}\fontname{Times new roman}SoC') % 第二个 y 轴标签
% set(gca,'YColor','k') % 右侧 y 轴刻度颜色为黑色
% set(gca,'ytick',0:0.2:1) % 右侧 y 轴刻度为 0:0.2:1
% 
% % 添加 x 轴标签
% xlabel('\fontsize{8}\fontname{Times new roman}Time (h)')
% 
% % 添加图例
% legend('Charge Bids','Discharge Bids','Risk-Aware LMP', 'SoC', 'Location', 'Best')
% 
% hold off

