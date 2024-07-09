% simulation days
a=0.3;
b=0.05;
c=1;
d=1:365;
e=0.2;
f=4;
g=0;
energy_price=zeros(24,365);
reserve_price=zeros(24,365);
storage_price=zeros(25,365);
SOC=zeros(25,365);
cost=zeros(1,365);
Dch_value=zeros(24,365);
Ddis_value=zeros(24,365);
Pch_value=zeros(24,365);
Pdis_value=zeros(24,365);
% SOC_value1_sorted=zeros(23,365);
% storage_price_sorted=zeros(23,365);

parfor i=1:365
[cost(:,i),energy_price(:,i),storage_price(:,i),reserve_price(:,i),SOC(:,i),Pdis_value(:,i),Pch_value(:,i),Dch_value(:,i),Ddis_value(:,i)]=MultiPeriod_ED_CC(a,b,c,d(i),e,f,g);
end

T=24;
energy_price_average=mean(energy_price,1);
reserve_price_total=sum(reserve_price,1);
storage_price_average=mean(storage_price(2:end,:),1);

% figure(1)
% set(gcf,'unit','centimeters','position',[0,0,8,6])
% for i=2
% % SOC_ch=SOC(intersect(find(Dch_value(:,i)==1)-1,1:T-1));
% % storage_price_ch=storage_price(intersect(find(Dch_value(:,i)==1),2:T));
% plot(SOC(i,:),storage_price(i,:),'o')
% hold on
% 
% end
% set(gca,'FontName','Times New Roman','FontSize',8)
% xlabel('\fontsize{8}\fontname{Times new roman}SoC')
% ylabel('\fontsize{8}\fontname{Times new roman}Opportunity Price ($/MWh)')

% % figure(2)
% % set(gcf,'unit','centimeters','position',[0,0,8,6])
% % for i=1:9
% % SOC_dis=SOC(intersect(find(Ddis_value(:,i)==1)-1,1:T-1));
% % storage_price_dis=storage_price(intersect(find(Ddis_value(:,i)==1),2:T));
% % plot(SOC_dis,storage_price_dis)
% % hold on
% % 
% % set(gca,'FontName','Times New Roman','FontSize',8)
% % xlabel('\fontsize{8}\fontname{Times new roman}SoC')
% % ylabel('\fontsize{8}\fontname{Times new roman}Opportunity Price ($/MWh)')
% % 
% % end
% % 
% % figure(3)
% % set(gcf,'unit','centimeters','position',[0,0,8,6])
% % for i=1:9
% % SOC_idle=SOC(intersect(intersect(find(Ddis_value(:,i)==0)-1,1:T-1),find(Dch_value(:,i)==0)-1));
% % storage_price_idle=storage_price(intersect(intersect(find(Ddis_value(:,i)==0),2:T),find(Dch_value(:,i)==0)));
% % plot(SOC_idle,storage_price_idle)
% % hold on
% % 
% % set(gca,'FontName','Times New Roman','FontSize',8)
% % xlabel('\fontsize{8}\fontname{Times new roman}SoC')
% % ylabel('\fontsize{8}\fontname{Times new roman}Opportunity Price ($/MWh)')
% % 
% % end
% 
% storage_price_average=mean(storage_price_sorted,1);
% energy_price_average=mean(energy_price,1);
% reserve_price_average=mean(reserve_price,1);
% 
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
% 
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
% B = zeros(365, T);
% 
% 
% for i = 1:365
%     reshapedRow = reshape(netload(i, :), 4, T);
%     B(i, :) = mean(reshapedRow, 1);
% end
% 
% netload=B;
% 
% netload_mean=mean(netload');
% 
% [netload_mean_sorted,idx] = sort(netload_mean, 'descend');
% storage_price_average_sorted1 = storage_price_average(idx);
% energy_price_average_sorted1 =energy_price_average(idx);
% reserve_price_average_sorted1=reserve_price_total(idx);
% cost_sorted1=cost(idx);
% 
% 
% figure(3)
% set(gcf,'unit','centimeters','position',[0,0,8,6])
% plot(netload_mean_sorted,storage_price_average_sorted1,'r-','LineWidth',1)
% set(gca,'FontName','Times New Roman','FontSize',8)
% xlabel('\fontsize{8}\fontname{Times new roman}Average Netload')
% ylabel('\fontsize{8}\fontname{Times new roman}Average Opportunity Price ($/MWh)')
% 
% 
% figure(4)
% set(gcf,'unit','centimeters','position',[0,0,8,6])
% plot(netload_mean_sorted,energy_price_average_sorted1,'r-','LineWidth',1)
% set(gca,'FontName','Times New Roman','FontSize',8)
% xlabel('\fontsize{8}\fontname{Times new roman}Average Netload')
% ylabel('\fontsize{8}\fontname{Times new roman}Average Electricity Price ($/MWh)')
% 
% 
% figure(5)
% set(gcf,'unit','centimeters','position',[0,0,8,6])
% plot(netload_mean_sorted,reserve_price_average_sorted1,'r-','LineWidth',1)
% set(gca,'FontName','Times New Roman','FontSize',8)
% xlabel('\fontsize{8}\fontname{Times new roman}Average Netload')
% ylabel('\fontsize{8}\fontname{Times new roman}Average Reserve Cost ($)')
% 
% figure(6)
% set(gcf,'unit','centimeters','position',[0,0,8,6])
% plot(netload_mean_sorted,cost_sorted1,'r-','LineWidth',1)
% set(gca,'FontName','Times New Roman','FontSize',8)
% xlabel('\fontsize{8}\fontname{Times new roman}Average Netload')
% ylabel('\fontsize{8}\fontname{Times new roman}System Cost ($)')
