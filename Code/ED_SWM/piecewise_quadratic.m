%@Create Time    : 2024/4/15
%@Author  : Ning Qi nq2716@columbia.edu
%@File    : piecewise_quadratic.m
%@Description : This is the code for fitting generator function
load data.mat

data = [data_generator(:,1),data_generator(:,8),data_generator(:,9),data_generator(:,10)];

% 按照二次系数从大到小排序
sortedData = sortrows(data, [4, 3, 2, -1]);


% 初始化累积容量和累积成本数组
all_P_shifted = [];
all_Cost = [];
cumulative_capacity = 0;

% 用于绘图的设置
figure;
hold on; % 允许在同一图形上绘制多个曲线
grid on;
xlabel('Cumulative Capacity (MW)');
ylabel('Cumulative Cost ($)');

% 遍历每个发电机，绘制其成本曲线
for i = 1:size(sortedData, 1)
    capacity = sortedData(i, 1);
    const = sortedData(i, 2);
    linear = sortedData(i, 3);
    quadratic = sortedData(i, 4);
    
    % 生成当前机组的功率范围
    P = linspace(0, capacity, 100);
    % 计算成本
    Cost = quadratic * P.^2 + linear * P + const;
    
    % 累积容量并调整横坐标
    P_shifted = P + cumulative_capacity;
    
    % 累积成本计算：当前机组的成本加上之前的累积成本
    if i > 1
        Cost = Cost + all_Cost(end);
    end
    
    % 绘制当前机组的成本曲线
    plot(P_shifted, Cost, 'LineWidth', 2);
    
    % 更新数据集合
    all_P_shifted = [all_P_shifted, P_shifted];
    all_Cost = [all_Cost, Cost];
    
    % 更新累积容量
    cumulative_capacity = cumulative_capacity + capacity;
end

% 图例和标题
title('Piecewise Quadratic Cost Curves');
%legend(arrayfun(@(i) sprintf('Generator %d', i), 1:size(data, 1), 'UniformOutput', false), 'Location', 'northwest');
hold off;

% 保存数据
save('generator_data.mat', 'all_P_shifted', 'all_Cost');

% 加载数据（如果之后需要）
% load('generator_data.mat', 'all_P_shifted', 'all_Cost');

% 定义二次、三次和四次模型函数
quadModel = @(b, x) b(1) * x.^2 + b(2) * x + b(3);
cubicModel = @(b, x) b(1) * x.^3 + b(2) * x.^2 + b(3) * x + b(4);
quarticModel = @(b, x) b(1) * x.^4 + b(2) * x.^3 + b(3) * x.^2 + b(4) * x + b(5);

% 初始猜测
initialGuessQuad = [0.0001, 0.1, 200];
initialGuessCubic = [0, 0.0001, 0.1, 200];
initialGuessQuartic = [0, 0, 0.0001, 0.1, 200];

% 设置非负约束
lb = zeros(size(initialGuessQuartic));  % 最小值设为0

% 使用 lsqcurvefit 进行拟合
options = optimset('Display', 'off');
bQuad = lsqcurvefit(quadModel, initialGuessQuad, all_P_shifted, all_Cost, lb(1:3), [], options);
bCubic = lsqcurvefit(cubicModel, initialGuessCubic, all_P_shifted, all_Cost, lb(1:4), [], options);
bQuartic = lsqcurvefit(quarticModel, initialGuessQuartic, all_P_shifted, all_Cost, lb, [], options);

% 计算拟合曲线的值
P_fit = linspace(min(all_P_shifted), max(all_P_shifted), 1000);
fitQuad = quadModel(bQuad, P_fit);
fitCubic = cubicModel(bCubic, P_fit);
fitQuartic = quarticModel(bQuartic, P_fit);

% 绘制原始数据点和拟合曲线
figure;
plot(all_P_shifted, all_Cost, 'ko', 'MarkerFaceColor', 'k');  % 原始数据点
hold on;
plot(P_fit, fitQuad, 'r-', 'LineWidth', 2);  % 二次拟合
plot(P_fit, fitCubic, 'b--', 'LineWidth', 2);  % 三次拟合
plot(P_fit, fitQuartic, 'g-.', 'LineWidth', 2);  % 四次拟合
xlabel('Cumulative Capacity (MW)');
ylabel('Cumulative Cost ($)');
title('Constrained Polynomial Fitting of Cumulative Cost Curve');
legend('Data', 'Quadratic Fit', 'Cubic Fit', 'Quartic Fit', 'Location', 'best');
grid on;
hold off;









