%% plot figure1
load result1.mat
X=size(epsilon,2);
figure(1)
set(gcf,'unit','centimeters','position',[0,0,5,3])
plot(energy_price','LineWidth',1)
set(gca,'xtick',[1:1:X])
set(gca,'xticklabel',{'0.005','0.01','0.05','0.1','0.15','0.2','0.25','0.3'})
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}\epsilon')
ylabel('\fontsize{8}\fontname{Times new roman}\lambda ($/MWh)')

figure(2)
set(gcf,'unit','centimeters','position',[0,0,5,3])
plot(storage_price','LineWidth',1)
set(gca,'xtick',[1:1:X])
set(gca,'xticklabel',{'0.005','0.01','0.05','0.1','0.15','0.2','0.25','0.3'})
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}\epsilon')
ylabel('\fontsize{8}\fontname{Times new roman}\theta ($/MWh)')

figure(3)
set(gcf,'unit','centimeters','position',[0,0,5,3])
plot(reserve_price','LineWidth',1)
set(gca,'xtick',[1:1:X])
set(gca,'xticklabel',{'0.005','0.01','0.05','0.1','0.15','0.2','0.25','0.3'})
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}\epsilon')
ylabel('\fontsize{8}\fontname{Times new roman}\pi ($)')
% 定义插图的位置和大小
insetPosition = [0.2, 0.6, 0.25, 0.25]; % [left, bottom, width, height]

% 创建插图
axes('Position', insetPosition);
box on; % 添加边框
hold on;

% 在插图中绘制曲线
plot(reserve_price(1:end-1,:)','LineWidth',1)
set(gca,'xticklabel',{'0.01','0.05','0.05','0.1'})
set(gca,'FontName','Times New Roman','FontSize',6)
% 设置插图的x轴范围为[0.05, 0.15]
xlim([1, 4]);

% 显示图形
hold off;

figure(4)
set(gcf,'unit','centimeters','position',[0,0,5,3])
plot(cost','LineWidth',1)
set(gca,'xtick',[1:1:X])
set(gca,'xticklabel',{'0.005','0.01','0.05','0.1','0.15','0.2','0.25','0.3'})
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}\epsilon')
ylabel('\fontsize{8}\fontname{Times new roman}Cost ($)')
% 定义插图的位置和大小
insetPosition = [0.2, 0.6, 0.25, 0.25]; % [left, bottom, width, height]

% 创建插图
axes('Position', insetPosition);
box on; % 添加边框
hold on;

% 在插图中绘制曲线
plot(cost(1:end-1,:)','LineWidth',1)
set(gca,'xticklabel',{'0.01','0.05','0.05','0.1'})
set(gca,'FontName','Times New Roman','FontSize',6)
% 设置插图的x轴范围为[0.05, 0.15]
xlim([1, 4]);

% 显示图形
hold off;