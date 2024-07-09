%% initial test
% storage bid
function [obj_value,energy_cost]=MultiPeriod_ED_ESbid(a,b,c,d,e,f,g)
%% inputs
%generation cost function
%net load

load data.mat
T=24; %time duration
%% constants
wind_ratio=a;% wind ratio
netload_error=load_error-max(max(load_real))*(wind_real-wind_forecast)*wind_ratio;

netload_error=reshape(netload_error',96,365)';% netload error, [-, +]

% 初始化一个新的 365x24 的数组 B
B = zeros(365, 24);

for i = 1:365
    
    reshapedRow = reshape(netload_error(i, :), 4, T);
    B(i, :) = mean(reshapedRow, 1);
end

netload_error=B;

epsilon=b; % probability level
for i=1:T
% [mu(i),sigma(i)] = normfit(netload_error(:,i),0.05); % μ, σ time different

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

sigma=sigma*c;



load_c=netload(d,:);

load_d1=zeros(1,T);
load_d2=zeros(1,T);

% alpha_value=zeros(1,T);
% beta_value=zeros(1,T);
% gama_value=zeros(1,T);
% for i=1:T
% [alpha_value(i),beta_value(i),gama_value(i)]=fit_versatile(B(:,i));
% end

for i=1:T

load_d1(:,i)= norminv(epsilon, mu(i), sigma(i));
% if epsilon<1/6
% load_d1(:,i)= mu(i)-sqrt(2/(9*epsilon))*sigma(i);
% else 
%     load_d1(:,i)= mu(i)-sqrt(3)*(1-2*epsilon)*sigma(i);
% end
% 
%load_d1(:,i)= mu(i)-sqrt((1-epsilon)/(epsilon))*sigma(i);
%load_d1(:,i)= quantile(B(:,i), epsilon);
% load_d1(:,i)=gama_value(i)-imag(log((epsilon)^(-1/beta_value(i))-1)./alpha_value(i));

end

for i=1:T

load_d2(:,i)= norminv(1-epsilon, mu(i), sigma(i));

% if epsilon<1/6
% load_d2(:,i)= mu(i)+sqrt(2/(9*epsilon))*sigma(i);
% else 
%     load_d2(:,i)= mu(i)+sqrt(3)*(1-2*epsilon)*sigma(i);
% end

%load_d2(:,i)= mu(i)+sqrt((1-epsilon)/(epsilon))*sigma(i);
%load_d2(:,i)= quantile(B(:,i), 1-epsilon);
% load_d2(:,i)=gama_value(i)-imag(log((1-epsilon)^(-1/beta_value(i))-1)./alpha_value(i));
end


num_e=1; %number of storage
num_g=size(data_generator,1); %number of generator

ESpowercap=e*max(max(netload)); % storage power cap
ESenergycap=ESpowercap*f; % storage energy cap
G_min=data_generator(:,2)*0;
G_max=data_generator(:,1)*0.8;
%G_min=zeros(num_g,1);
% G_min=0;%% assume one generator
% G_max=sum(G_max);
% C_0=0*ones(num_g,1);
% C_1=43*ones(num_g,1);
% C_2=0.00292*ones(num_g,1)*0;
% C_0=10707*ones(num_g,1);
% C_1=153.16*ones(num_g,1);
% C_2=0.0346*ones(num_g,1);

C_0=data_generator(:,8);
C_1=data_generator(:,9);
C_2=data_generator(:,10);
% C_2=zeros(num_g,1);
% load generator_fit.mat
% C_0=bQuad(:,3);
% C_1=bQuad(:,2);
% C_2=bQuad(:,1);
Pg=sdpvar(num_g,T); % generator power
rho_g=sdpvar(num_g,T); % generator reserve ratio
%% constraints
cons=[];

for j=1:T    
cons=[cons, sum(Pg(:,j),1)+g(j)*ESenergycap==load_c(:,j)];
end

for j=1:T    
cons=[cons,sum(rho_g(:,j),1)==1];
end


for j=1:T   
cons=[cons,0<=rho_g(:,j)<=1];
end
% 
% 
% 
for j=1:T    
cons=[cons,G_min<=Pg(:,j)+rho_g(:,j)*load_d1(1,j)];
end

for j=1:T    
cons=[cons,Pg(:,j)+rho_g(:,j)*load_d2(1,j)<=G_max];
end


%% objective
obj=0;

% % quadratic
for i=1:num_g
obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:)+rho_g(i,:).*mu)+C_2(i,1)*(Pg(i,:).*Pg(i,:)+2*Pg(i,:).*rho_g(i,:).*mu+rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma));
end
% for i=1:num_g
% obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:))+C_2(i)*(Pg(i,:).*Pg(i,:));
% end
% cubic
% for i=1:num_g
% obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:)+rho_g(i,:).*mu)+C_2(i)*(Pg(i,:).*Pg(i,:)+2*Pg(i,:).*rho_g(i,:).*mu+rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma))+C_3(i)*(Pg(i,:).*Pg(i,:).*Pg(i,:)+3*Pg(i,:).*Pg(i,:).*rho_g(i,:).*mu+3*Pg(i,:).*rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma));
% end
% % quartic
% for i=1:num_g
% obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:)+rho_g(i,:).*mu)+C_2(i)*(Pg(i,:).*Pg(i,:)+2*Pg(i,:).*rho_g(i,:).*mu+rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma))+C_3(i)*(Pg(i,:).*Pg(i,:).*Pg(i,:)+3*Pg(i,:).*Pg(i,:).*rho_g(i,:).*mu+3*Pg(i,:).*rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma))+C_4(i)*(Pg(i,:).*Pg(i,:).*Pg(i,:).*Pg(i,:)+4*Pg(i,:).*Pg(i,:).*Pg(i,:).*rho_g(i,:).*mu+6*Pg(i,:).*Pg(i,:).*rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma)+rho_g(i,:).*rho_g(i,:).*rho_g(i,:).*rho_g(i,:).*(3*sigma.*sigma.*sigma.*sigma+6*sigma.*sigma.*mu.*mu+mu.*mu.*mu.*mu));
% end

obj=sum(obj,2);
options = sdpsettings('verbose', 1, 'solver', 'gurobi');
results= optimize(cons,obj,options);


obj_value=value(obj);
energy_price=abs(dual(cons(1:T)));
energy_cost=sum(load_c*energy_price);

% C1=yita*yita*load_d1./load_d2;
% C2=yita*(load_d2/(yita*yita)-load_d1)./load_d2;
% C3=yita./load_d2;
% C4=-20*mu*yita./load_d2;
% 
% B1=load_d2./(load_d1*yita*yita);
% B2=(load_d1*yita*yita-load_d2)./(yita*load_d1);
% B3=1./(load_d1*yita);
% B4=20*(load_d2-yita*yita*load_d1-mu)./(load_d1*yita);
% 
% charge_price=C1.*storage_price(2:end)'+C2.*energy_price'-C3.*reserve_price'+C4;
% discharge_price=B1.*storage_price(2:end)'+B2.*energy_price'-B3.*reserve_price'+B4;
% discharge_price(end)=storage_price(end-1);
% figure(1)
% set(gcf,'unit','centimeters','position',[0,0,5,3])
% plot(storage_price(1:end-1),'k-x','LineWidth',1,'MarkerSize',2)
% hold on
% % plot(charge_price,'r-.','LineWidth',1)
% % hold on
% plot(discharge_price,'b-.','LineWidth',1)
% hold on
% set(gca,'FontName','Times New Roman','FontSize',8)
% xlabel('\fontsize{8}\fontname{Times new roman}Time (h)')
% ylabel('\fontsize{8}\fontname{Times new roman}\theta_t ($/MWh)')

% SOC_ch=SOC_value(intersect(find(Dch_value==1),1:T));
% SOC_dis=SOC_value(intersect(find(Ddis_value==1),1:T));
% SOC_idle=SOC_value(intersect(intersect(find(Ddis_value==0),1:T),find(Dch_value==0)));
% 
% storage_price_dis=storage_price(intersect(find(Ddis_value==1),1:T));
% storage_price_ch=storage_price(intersect(find(Dch_value==1),1:T));
% storage_price_idle=storage_price(intersect(intersect(find(Ddis_value==0),1:T),find(Dch_value==0)));
% Dch_value=Dch_value';
% Ddis_value=Ddis_value';
% % 
% % [SOC_value1_sorted, idx] = sort(SOC_value1(1:T-1), 'descend');
% % storage_price_sorted = storage_price(idx+1);
% % 
% [SOC_ch_sorted, idx] = sort(SOC_ch, 'descend');
% storage_price_ch_sorted = storage_price_ch(idx);
% 
% 
% [SOC_dis_sorted, idx] = sort(SOC_dis, 'descend');
% storage_price_dis_sorted = storage_price_dis(idx);
% 
% [SOC_idle_sorted, idx] = sort(SOC_idle, 'descend');
% storage_price_idle_sorted = storage_price_idle(idx);
% 
% [SOC_value_sorted, idx] = sort(SOC_value(1:T), 'descend');
% storage_price_sorted = storage_price(idx);
% % % 
% figure(1)
% set(gcf,'unit','centimeters','position',[0,0,8,6])
% % plot(SOC_value_sorted,storage_price_sorted)
% plot(SOC_dis_sorted,storage_price_dis_sorted,'o')
% hold on
% plot(SOC_ch_sorted,storage_price_ch_sorted,'o')
% hold on
% plot(SOC_idle_sorted,storage_price_idle_sorted,'o')
% % % plot(SOC_ch,storage_price_ch)
% % % hold on
% % % plot(SOC_dis,storage_price_dis)
% % % hold on
% % % plot(SOC_idle,storage_price_idle)
% set(gca,'FontName','Times New Roman','FontSize',8)
% xlabel('\fontsize{8}\fontname{Times new roman}SoC')
% ylabel('\fontsize{8}\fontname{Times new roman}Opportunity Price ($/MWh)')
% 
% 
% chargeprice=mean(storage_price_ch_sorted)
% dischargeprice=mean(storage_price_dis_sorted )