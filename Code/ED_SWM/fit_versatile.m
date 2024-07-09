
function [alpha,beta,gama]=fit_versatile(data)
[F, x_values] = ecdf(data);

modelFun = @(b, x) (1 + exp(-max(0.001, b(1)) * (x - b(2)))).^(-max(0.001, b(3)));

% 选择合理的初始猜测
initialGuess = [0.1, median(x_values), 0.1];  % 假设 b(1) 和 b(3) 不应过小

% 使用 nlinfit 拟合模型
options = statset('nlinfit');
options.RobustWgtFun = 'bisquare';  % 使用鲁棒性权重函数
[beta_hat, ~, ~, ~, ~] = nlinfit(x_values, F, modelFun, initialGuess, options);
alpha=beta_hat(1);
beta=beta_hat(2);
gama=beta_hat(3);
% 
%     % 打印拟合参数
%     disp('Fitted parameters:');
%     disp(['Alpha = ', num2str(beta_hat(1))]);
%     disp(['Gamma = ', num2str(beta_hat(2))]);
%     disp(['Beta = ', num2str(beta_hat(3))]);
% 
%     % 计算拟合模型的 CDF
%     fittedCDF = modelFun(beta_hat, x_values);
% 
%     % 绘制经验 CDF 和拟合 CDF
%     figure;
%     hold on;
%     plot(x_values, F, 'o', 'DisplayName', 'Empirical CDF');
%     plot(x_values, fittedCDF, '-', 'DisplayName', 'Fitted CDF');
%     xlabel('Data values');
%     ylabel('Cumulative probability');
%     title('Empirical vs. Fitted CDF');
%     legend show;
%     grid on;
%     hold off;