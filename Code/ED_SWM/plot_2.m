storage_price=[26.30285511	25.82626085	24.3352351	23.33666709	22.55169436	22.00626923	21.61824809	21.24584309	20.86528509	20.422772	6.085627548
27.01857444	26.75488895	26.48547054	26.18565374	25.64349746	25.08455863	24.60495577	24.10145055	23.64336651	22.91468373	5.153441325
27.97474114	27.37436075	27.05630544	26.80289966	26.57240263	26.34780991	26.06919262	25.66410904	25.20223316	13.00197612	3.818864837
29.50462948	27.74755994	27.47927507	27.22778401	26.98453223	26.77190096	26.55827509	26.344895	26.11611939	8.793104587	2.835327946
37.74819015	28.33091239	27.84808295	27.49140468	27.24145841	27.01466233	26.79181399	26.51665203	21.81043789	3.912949185	2.324691711
39.8931152	29.34103974	28.60228912	27.82768248	27.19516093	26.81603782	26.51531231	26.35984649	9.508677868	3.303102095	2.270774505
];

g=0:0.1:1;
figure(1)
set(gcf,'unit','centimeters','position',[0,0,8,4])
plot(g,storage_price,'LineWidth',1)
hold on
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}SoC')
ylabel('\fontsize{8}\fontname{Times new roman}\theta_{t} ($/MWh)')
insetPosition = [0.2, 0.6, 0.25, 0.25]; % [left, bottom, width, height]

% 创建插图
axes('Position', insetPosition);
box on; % 添加边框
hold on;
plot(g,storage_price,'LineWidth',1)
set(gca,'xticklabel',{'0.1','0.2','0.3','0.4'})
set(gca,'FontName','Times New Roman','FontSize',6)
% 设置插图的x轴范围为[0.05, 0.15]
xlim([0.1, 0.4]);

insetPosition = [0.2, 0.6, 0.25, 0.25]; % [left, bottom, width, height]

% 创建插图
axes('Position', insetPosition);
box on; % 添加边框
hold on;
plot(g,storage_price,'LineWidth',1)
set(gca,'xticklabel',{'0.6','0.7','0.8','0.9'})
set(gca,'FontName','Times New Roman','FontSize',6)
% 设置插图的x轴范围为[0.05, 0.15]
xlim([0.6, 0.9]);
