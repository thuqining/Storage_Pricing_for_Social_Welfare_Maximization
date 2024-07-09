a=[155805.5091	204534.9104	8905814.265	9185946.985	8986265.125	9283581.939	14127000.63	15160147.12
148401.7001	211682.0732	8918452.9	9199974.315	9006235.321	9297856.55	14014072.71	15813716.47
142488.5444	203387.4009	8937623.213	9211966.48	9031018.92	9309848.715	13818305.85	16109186.22
279311.1546	376780.6903	8966622.489	9202349.224	9064504.724	9300231.459	15961793.5	17254209
279565.3816	377030.6036	9043288.242	9325646.385	9141170.477	9423528.621	15811717.06	18428247.76
275804.6038	369916.126	9126937.849	9500929.435	9224820.084	9598825.408	15840186.28	19180398.44
292678.5414	390292.8836	9218666.057	9617722.469	9316548.292	9715604.704	15968903.97	19331219.15
501421.1122	598461.7386	9354419.89	9784241.966	9452302.125	9882124.201	19222323.02	19382490.75
490488.8273	587655.9707	9558664.667	9953821.703	9656546.902	10051703.94	19239576.99	19655084.3
501557.7814	589434.3209	9783068.317	10202375.54	9880950.552	10300271.51	20171024.41	20720712.53
427659.4966	525451.9173	10035786.01	10358006.62	10133668.24	10455888.85	21138252.93	20599961.78
];
c=0.5:0.2:2.5;
figure(1)
set(gcf,'unit','centimeters','position',[0,0,5,3])
plot(c,a(:,1),'LineWidth',1)
hold on
plot(c,a(:,2),'LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Scale Factor of \sigma')
ylabel('\fontsize{8}\fontname{Times new roman}Storage Profit ($)')
xlim([0.5 2.5])
set(gca,'xticklabel',{'0.5','0.9','1.3','1.7','2.1','2.5'})


figure(2)
set(gcf,'unit','centimeters','position',[0,0,5,3])
plot(c,a(:,3),'LineWidth',1)
hold on
plot(c,a(:,4),'LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Scale Factor of  \sigma')
ylabel('\fontsize{8}\fontname{Times new roman}Generation Cost ($)')
xlim([0.5 2.5])
set(gca,'xticklabel',{'0.5','0.9','1.3','1.7','2.1','2.5'})


figure(3)
set(gcf,'unit','centimeters','position',[0,0,5,3])
plot(c,a(:,5),'LineWidth',1)
hold on
plot(c,a(:,6),'LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Scale Factor of  \sigma')
ylabel('\fontsize{8}\fontname{Times new roman}System Cost ($)')
xlim([0.5 2.5])
set(gca,'xticklabel',{'0.5','0.9','1.3','1.7','2.1','2.5'})

figure(4)
set(gcf,'unit','centimeters','position',[0,0,5,3])
plot(c,a(:,7),'LineWidth',1)
hold on
plot(c,a(:,8),'LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Scale Factor of  \sigma')
ylabel('\fontsize{8}\fontname{Times new roman}Electricity Payment ($)')
xlim([0.5 2.5])
set(gca,'xticklabel',{'0.5','0.9','1.3','1.7','2.1','2.5'})
