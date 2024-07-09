% storage power
a=0.3;
b=0.05;
c=1;
d=10;
e=[0.2 0.6 1];
f=8;
g=0.5;
X=size(e,2);
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

parfor i=1:X
[cost(:,i),energy_price(:,i),storage_price(:,i),reserve_price(:,i),SOC(:,i),~]=MultiPeriod_ED_CC(a,b,c,d,e(i),f,g);
end
% parfor i=1:X
% [cost(:,i),energy_price(:,i),reserve_price(:,i),~]=MultiPeriod_ED_CC_withoutES(a,b,c,d,e(i),f,g);
% end

energy_price_average=mean(energy_price,1);
reserve_price_total=sum(reserve_price,1);
storage_price_average=mean(storage_price(2:end,:),1);

T=24;

figure(1)
set(gcf,'unit','centimeters','position',[0,0,5,3])
% for i=1:X
plot(e,storage_price_average,'.-','LineWidth',1)
hold on
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Storage Capacity (%)')
ylabel('\fontsize{8}\fontname{Times new roman}\theta_{t}($/MWh)')
% end

figure(2)
set(gcf,'unit','centimeters','position',[0,0,5,3])
% for i=1:X
plot(e,energy_price_average,'.-','LineWidth',1)
hold on
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Storage Capacity (%)')
ylabel('\fontsize{8}\fontname{Times new roman}\lambda_{t}($/MWh)')
% end
figure(3)
set(gcf,'unit','centimeters','position',[0,0,5,3])
% for i=1:X
plot(e,reserve_price_total,'.-','LineWidth',1)
hold on
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Storage Capacity (%)')
ylabel('\fontsize{8}\fontname{Times new roman}\pi_{t} ($)')
% end
figure(4)
set(gcf,'unit','centimeters','position',[0,0,5,3])
% for i=1:X
plot(e,cost,'.-','LineWidth',1)
hold on
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Storage Capacity (%)')
ylabel('\fontsize{8}\fontname{Times new roman}Cost ($)')
% end


