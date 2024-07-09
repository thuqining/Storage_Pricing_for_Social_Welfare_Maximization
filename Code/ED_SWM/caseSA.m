%% YITA
a=[0.1 0.2 0.3 0.4 0.5];%wind ratio
b=0.05;%risk-aversion
c=1;%uncertainty
d=1;%load capacity
e=0.2;%power capacity
f=4;%duration
g=0.5;%initial SoC
m=0.95;%efficiency
n=20;%marginal cost
X=size(a,2);
energy_price=zeros(24,X);
reserve_price=zeros(24,X);
storage_price=zeros(25,X);
cost=zeros(1,X);
storageprofit=zeros(1,X);
energy_cost=zeros(1,X);
generationcost=zeros(1,X);

parfor i=1:X
[cost(:,i),energy_price(:,i),storage_price(:,i),reserve_price(:,i),storageprofit(:,i),energy_cost(:,i),generationcost(:,i)]=MultiPeriod_ED_CC_SA(a(i),b,c,d,e,f,g,m,n);
end
energy_price_average=mean(energy_price,1);
reserve_price_total=sum(reserve_price,1);
storage_price_average=mean(storage_price(2:end,:),1);
% reserve_price1=cost-sum(storage_price(2:end,:).*(Pdis/0.95-Pch*0.95),1)-sum(Pg.*energy_price,1);


% figure(1)
% set(gcf,'unit','centimeters','position',[0,0,8,4])
% for i=1
% % SOC_ch=SOC(intersect(find(Dch_value(:,i)==1)-1,1:T-1));
% % storage_price_ch=storage_price(intersect(find(Dch_value(:,i)==1),2:T));
% plot(SOC(i,:),storage_price(i,:))
% hold on
% 
% end
% set(gca,'FontName','Times New Roman','FontSize',8)
% xlabel('\fontsize{8}\fontname{Times new roman}SoC')
% ylabel('\fontsize{8}\fontname{Times new roman}Opportunity Price ($/MWh)')

