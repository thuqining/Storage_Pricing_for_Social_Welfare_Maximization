% wind generation case
a=0.1:0.1:0.9;
b=0.01;
c=1;
d=10;
e=0.2;
f=4;
g=0.5;
X=size(a,2);
energy_price=zeros(24,X);
reserve_price=zeros(24,X);
storage_price=zeros(24,X);
SOC=zeros(25,X);
cost=zeros(1,X);
Dch_value=zeros(24,X);
Ddis_value=zeros(24,X);
SOC_value1_sorted=zeros(23,X);
storage_price_sorted=zeros(23,X);

parfor i=1:X
[cost(:,i),energy_price(:,i),storage_price(:,i),reserve_price(:,i),SOC(:,i),~]=MultiPeriod_ED_CC(a(i),b,c,d,e,f,g);
end

T=24;

figure(1)
set(gcf,'unit','centimeters','position',[0,0,8,6])
for i=1:X
plot(SOC_value1_sorted(:,i),storage_price_sorted(:,i))
hold on
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}SoC')
ylabel('\fontsize{8}\fontname{Times new roman}Opportunity Price ($/MWh)')
end

% figure(2)
% set(gcf,'unit','centimeters','position',[0,0,8,6])
% for i=1:9
% SOC_dis=SOC(intersect(find(Ddis_value(:,i)==1)-1,1:T-1));
% storage_price_dis=storage_price(intersect(find(Ddis_value(:,i)==1),2:T));
% plot(SOC_dis,storage_price_dis)
% hold on
% 
% set(gca,'FontName','Times New Roman','FontSize',8)
% xlabel('\fontsize{8}\fontname{Times new roman}SoC')
% ylabel('\fontsize{8}\fontname{Times new roman}Opportunity Price ($/MWh)')
% 
% end
% 
% figure(3)
% set(gcf,'unit','centimeters','position',[0,0,8,6])
% for i=1:9
% SOC_idle=SOC(intersect(intersect(find(Ddis_value(:,i)==0)-1,1:T-1),find(Dch_value(:,i)==0)-1));
% storage_price_idle=storage_price(intersect(intersect(find(Ddis_value(:,i)==0),2:T),find(Dch_value(:,i)==0)));
% plot(SOC_idle,storage_price_idle)
% hold on
% 
% set(gca,'FontName','Times New Roman','FontSize',8)
% xlabel('\fontsize{8}\fontname{Times new roman}SoC')
% ylabel('\fontsize{8}\fontname{Times new roman}Opportunity Price ($/MWh)')
% 
% end

storage_price_average=mean(storage_price_sorted,1);
energy_price_average=mean(energy_price,1);
reserve_price_average=mean(reserve_price,1);


figure(2)
set(gcf,'unit','centimeters','position',[0,0,8,6])
plot(1:X,storage_price_average,'r-','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}\epsilon')
ylabel('\fontsize{8}\fontname{Times new roman}Average Opportunity Price ($/MWh)')
set(gca,'xtick',[1:1:X])
set(gca,'xticklabel',{'0.001','0.01','0.05','0.1','0.15','0.2','0.25','0.3','0.35'})


figure(3)
set(gcf,'unit','centimeters','position',[0,0,8,6])
plot(1:X,energy_price_average,'r-','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}\epsilon')
ylabel('\fontsize{8}\fontname{Times new roman}Average Electricity Price ($/MWh)')
set(gca,'xtick',[1:1:X])
set(gca,'xticklabel',{'0.001','0.01','0.05','0.1','0.15','0.2','0.25','0.3','0.35'})


figure(4)
set(gcf,'unit','centimeters','position',[0,0,8,6])
plot(1:X,reserve_price_average,'r-','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}\epsilon')
ylabel('\fontsize{8}\fontname{Times new roman}Average Reserve Cost ($)')
set(gca,'xtick',[1:1:X])
set(gca,'xticklabel',{'0.001','0.01','0.05','0.1','0.15','0.2','0.25','0.3','0.35'})


figure(5)
set(gcf,'unit','centimeters','position',[0,0,8,6])
plot(1:X,cost,'r-','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}\epsilon')
ylabel('\fontsize{8}\fontname{Times new roman}System Cost ($)')
set(gca,'xtick',[1:1:X])
set(gca,'xticklabel',{'0.001','0.01','0.05','0.1','0.15','0.2','0.25','0.3','0.35'})
