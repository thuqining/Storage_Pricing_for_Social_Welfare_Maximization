%% sigma
a=0.3;
b=0.05;
c=[0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.9 2.1 2.3 2.5];
d=10;
e=0.2;
f=4;
g=0.5;
X=size(c,2);
energy_price=zeros(24,X);
reserve_price=zeros(24,X);
storage_price=zeros(25,X);
SOC=zeros(25,X);
cost=zeros(1,X);
rho_g=zeros(1,24,X);
rho_e=zeros(24,X);
Pdis=zeros(24,X);
Pch=zeros(24,X);
Pg=zeros(24,X);


parfor i=1:X
[cost(:,i),energy_price(:,i),storage_price(:,i),reserve_price(:,i),SOC(:,i),~]=MultiPeriod_ED_CC(a,b,c(i),d,e,f,g);
end

% storage_price(find(abs(SOC)<1e-6))=0;
energy_price_average=mean(energy_price,1);
reserve_price_total=sum(reserve_price,1);
storage_price_average=mean(storage_price(2:end,:),1);

T=24;

figure(1)
set(gcf,'unit','centimeters','position',[0,0,8,4])
plot(c,storage_price_average,'LineWidth',1)
hold on
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Scale Factor of \sigma')
ylabel('\fontsize{8}\fontname{Times new roman}\theta_{t} ($/MWh)')
xlim([0.5 2.5])
set(gca,'xticklabel',{'0.5','0.9','1.3','1.7','2.1','2.5'})

% insetPosition = [0.2, 0.6, 0.25, 0.25]; % [left, bottom, width, height]
% 
% % 创建插图
% axes('Position', insetPosition);
% box on; % 添加边框
% hold on;
% plot(g,storage_price(2,:),'LineWidth',1)
% set(gca,'xticklabel',{'0.1','0.2','0.3','0.4'})
% set(gca,'FontName','Times New Roman','FontSize',6)
% % 设置插图的x轴范围为[0.05, 0.15]
% xlim([0.1, 0.4]);

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

% storage_price_average=mean(storage_price_sorted,1);
% 
% figure(2)
% set(gcf,'unit','centimeters','position',[0,0,8,6])
% plot(c,storage_price_average,'r-','LineWidth',1)
% set(gca,'FontName','Times New Roman','FontSize',8)
% xlabel('\fontsize{8}\fontname{Times new roman}Scale Factor of \sigma')
% ylabel('\fontsize{8}\fontname{Times new roman}Average Opportunity Price ($/MWh)')


