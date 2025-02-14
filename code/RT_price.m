%@Create Time    : 2024/3/15
%@Author  : Ning Qi nq2176@columbia.edu
%@File    : MultiPeriod_ED_CC.m
%@Description : This is the main code for multi-period chance-constrained
%economic dispatch for system operator and can generate both dispatch and
%price for storage and generator
% storage: energy and reserve market

function [system_cost1_av,storage_profit1_av,consumer_payment1_av,system_cost2_av,storage_profit2_av,consumer_payment2_av]=RT_price(a,b,c,d,e,f,g,h,bd1,bc1,bd2,bc2)

%% inputs
% a=0.3;
% b=100;
% c=3; %sigma, original forecast error is too small
% d=10; %day number
% e=0.2; %storage capacity ratio (load capacity)
% f=4; %storage duration
% g=0; %initial SoC
% h=10; %storage marginal cost
% l=0; %mu=0,shifting load forecast error
% bd1=60*ones(5,24);
% bc1=0*ones(5,24);
% % bd2=30*ones(1,24);
% % bc2=0*ones(1,24);
% % bd3=40*ones(1,24);
% % bc3=0*ones(1,24);

MS=h;

load data.mat
%% case setting
T=24; %time duration
wind_ratio=a;% wind ratio
netload_error=load_error-max(max(load_real))*(wind_real-wind_forecast)*wind_ratio;
netload_error=reshape(netload_error',96,365)';% netload error, [-, +]
B = zeros(365, 24);
for i = 1:365
    reshapedRow = reshape(netload_error(i, :), 4, T);
    B(i, :) = mean(reshapedRow, 1);
end
netload_error=B;

epsilon=b; % probability level

for i=1:T
mu(i)=mean(netload_error(:,i));
sigma(i)=std(netload_error(:,i));
end

wind_forecast=reshape(wind_forecast',96,365)';
netload=load_forecast-max(max(load_real))*wind_forecast*wind_ratio; % netload forcast
C = zeros(365, T);
for i = 1:365
    reshapedRow = reshape(netload(i, :), 4, T);
    C(i, :) = mean(reshapedRow, 1);
end
netload=C;

sigma=sigma*c;% change the STD of forecast error for net load

load_c=netload(d,:);

load_scenario=zeros(b,T);

rng(42);

for i=1:b

load_scenario(i,:)=normrnd(mu,sigma);

end

load_scenario=load_scenario+ones(b,1)*load_c;


num_seg=size(bd1,1); %number of segment
num_g=size(data_generator,1); %number of generator
yita=0.95; % storage efficiency
% num_e=2; %number of storage
% yita=[0.8  0.95 ]'; % storage efficiency
SOC0=g./num_seg; % initial SoC
ESpowercap=e*max(max(netload)); % storage power cap
ESenergycap=ESpowercap*f; % storage energy cap
G_min=data_generator(:,2)*0;% minimum capacity=0
G_max=data_generator(:,1)*1;
RU=data_generator(:,6);
RD=data_generator(:,6);
C_0=data_generator(:,8);
C_1=data_generator(:,9);
C_2=data_generator(:,10);
delta_E=1./num_seg;

system_cost1 = zeros(b, 1);
storage_profit1= zeros(b, 1);
consumer_payment1=zeros(b, 1);


system_cost2 = zeros(b, 1);
storage_profit2= zeros(b, 1);
consumer_payment2=zeros(b, 1);


% system_cost3 = zeros(b, 1);
% storage_profit3= zeros(b, 1);

%% DP bids
parfor s = 1:b
    %% decision variables
    Pch = sdpvar(num_seg, T); % storage charge
    Pdis = sdpvar(num_seg, T); % storage discharge
    Dch = binvar(1, T); % storage charge state
    Dseg=binvar(num_seg-1, T); % storage segment
    SOC = sdpvar(num_seg, T); % storage SOC
    Pg = sdpvar(num_g, T); % generator power
    rg=sdpvar(num_g, T); % generator power
    lu=sdpvar(1, T); % generator power
    ld=sdpvar(1, T); % generator power


    %% constraints
    cons = [];
    % power balance-energy price dual
    for j = 1:T    
        cons = [cons, sum(Pg(:, j), 1) + sum(Pdis(:, j), 1) - sum(Pch(:, j), 1)+ld(:,j) == load_scenario(s, j)+lu(:,j)];
    end
    cons = [cons,lu>=0];
    cons = [cons,ld>=0];

   % storage SoC bounds
    for i = 1:num_seg
        for j = 1:T
            if j==1
            cons = [cons, SOC(i, j) == SOC0 + ((Pch(i, j)) .* yita + (-Pdis(i, j)) ./ yita) / ESenergycap];
            else
            cons = [cons, SOC(i, j) == SOC(i, j-1) + ((Pch(i, j)) .* yita + (-Pdis(i, j)) ./ yita) / ESenergycap];
            end
        end
    end
  
    % storage power bounds
    for j = 1:T    
        cons = [cons, Pch(:, j) >= 0];
        cons = [cons, Pdis(:, j) >= 0];
        cons = [cons, Pch(:, j) <= ESpowercap];
        cons = [cons, Pdis(:, j) <= ESpowercap];
        cons = [cons, sum(Pch(:, j)./ESpowercap,1) <= Dch(1, j)];
        cons = [cons, sum(Pdis(:, j)./ESpowercap,1) <=1-Dch(1, j)];
    end


    for j=1:T
    
        for i=2:num_seg
    
        cons = [cons, SOC(i, j) <=delta_E*Dseg(i-1,j)];    
    
        end
    
        for i=1:num_seg-1
    
        cons = [cons, SOC(i, j) >=delta_E*Dseg(i,j)];    
    
        end
    
       cons = [cons, SOC(num_seg, j)>=0];
       cons = [cons, SOC(1, j)<=delta_E];
       cons=[cons, sum(SOC(:,j),1)<=1];
    
    end


       for j=1:T

           cons = [cons, SOC(:, j)>=0];

       end



    % generator power bounds
for j=1:T    
cons=[cons,G_min<=Pg(:,j)<=G_max-rg(:,j)];
end


% generator reserve
for j=1:T    
cons=[cons,0<=rg(:,j)];
end


for j=1:T    
cons=[cons,sum(rg(:,j),1)>=0.1*load_c(:,j)];
end



    % generator RAMPUP and RAMPDOWN
    for j = 2:T    
        cons = [cons, Pg(:, j) - Pg(:, j-1) <= RU];
        cons = [cons, -RD <= Pg(:, j) - Pg(:, j-1)];
    end


    %% objective
    obj = 0;
    % quadratic + marginal cost
    for i = 1:num_g
        obj = obj + C_0(i, 1) + C_1(i, 1) * Pg(i, :) + C_2(i, 1) * (Pg(i, :) .* Pg(i, :));
    end
    obj = sum(obj, 2);


    for j = 1:T
        obj = obj + sum(bd1(:,j).*Pdis(:,j)-bc1(:,j).*Pch(:,j),1);
    end

    for j=1:T
        obj = obj + 2000*ld(:,j);
    end


    % optimization
    options = sdpsettings('verbose', 1, 'solver', 'gurobi');
    results = optimize(cons, obj, options);

    Dch_value=value(Dch);
    Dseg_value=value(Dseg);

   %% dual constraints
    cons = [];
    % power balance-energy price dual
    for j = 1:T    
        cons = [cons, sum(Pg(:, j), 1) + sum(Pdis(:, j), 1) - sum(Pch(:, j), 1)+ld(:,j) == load_scenario(s, j)+lu(:,j)];
    end
    cons = [cons,lu>=0];
    cons = [cons,ld>=0];


    for i = 1:num_seg
        for j = 1:T
            if j==1
            cons = [cons, SOC(i, j) == SOC0 + ((Pch(i, j)) .* yita + (-Pdis(i, j)) ./ yita) / ESenergycap];
            else
            cons = [cons, SOC(i, j) == SOC(i, j-1) + ((Pch(i, j)) .* yita + (-Pdis(i, j)) ./ yita) / ESenergycap];
            end
        end
    end
  
    % storage power bounds
    for j = 1:T    
        cons = [cons, Pch(:, j) >= 0];
        cons = [cons, Pdis(:, j) >= 0];
        cons = [cons, Pch(:, j) <= ESpowercap];
        cons = [cons, Pdis(:, j) <= ESpowercap];
        cons = [cons, sum(Pch(:, j)./ESpowercap,1) <= Dch_value(1, j)];
        cons = [cons, sum(Pdis(:, j)./ESpowercap,1) <=1-Dch_value(1, j)];
    end


    for j=1:T
    
        for i=2:num_seg
    
        cons = [cons, SOC(i, j) <=delta_E*Dseg_value(i-1,j)];    
    
        end
    
        for i=1:num_seg-1
    
        cons = [cons, SOC(i, j) >=delta_E*Dseg_value(i,j)];    
    
        end
    
       cons = [cons, SOC(num_seg, j)>=0];
       cons = [cons, SOC(1, j)<=delta_E];
       cons=[cons, sum(SOC(:,j),1)<=1];
    
    end


       for j=1:T

           cons = [cons, SOC(:, j)>=0];

       end


    % generator power bounds
for j=1:T    
cons=[cons,G_min<=Pg(:,j)<=G_max-rg(:,j)];
end


% generator reserve
for j=1:T    
cons=[cons,0<=rg(:,j)];
end


for j=1:T    
cons=[cons,sum(rg(:,j),1)>=0.1*load_c(:,j)];
end




    % generator RAMPUP and RAMPDOWN
    for j = 2:T    
        cons = [cons, Pg(:, j) - Pg(:, j-1) <= RU];
        cons = [cons, -RD <= Pg(:, j) - Pg(:, j-1)];
    end


    %% objective
    obj = 0;
    % quadratic generation
    for i = 1:num_g
        obj = obj + C_0(i, 1) + C_1(i, 1) * Pg(i, :) + C_2(i, 1) * (Pg(i, :) .* Pg(i, :));
    end
    obj = sum(obj, 2);


    % storage bids
    for j = 1:T
        obj = obj + sum(bd1(:,j).*Pdis(:,j)-bc1(:,j).*Pch(:,j),1);
    end

    for j=1:T
        obj = obj + 2000*ld(:,j);
    end


    % optimization
    options = sdpsettings('verbose', 1, 'solver', 'gurobi');
    results = optimize(cons, obj, options);


    energy_price=abs(dual(cons(1:T)));

    % store result in a local variable
    temp_cost1 = 0;
    temp_cost2=0;
    temp_cost3=0;
    

    for i = 1:num_g
        temp_cost1 = temp_cost1 + C_0(i, 1) + C_1(i, 1) * value(Pg(i, :)) + C_2(i, 1) * (value(Pg(i, :)) .* value(Pg(i, :)));
    end

    temp_cost1  = sum(temp_cost1 , 2);
    

    temp_cost1 =temp_cost1 +sum(MS*(value(sum(Pdis,1))));

    temp_cost2 = sum(energy_price'.*sum((value(Pdis-Pch)),1)-MS*(value(sum(Pdis,1))));
    temp_cost3=sum(energy_price'.*load_scenario(s,:));


       
    system_cost1(s, 1) = temp_cost1;
    storage_profit1(s, 1) = temp_cost2;
    consumer_payment1(s, 1)=temp_cost3;


end

% %% capped

parfor s = 1:b
    %% decision variables
    Pch = sdpvar(num_seg, T); % storage charge
    Pdis = sdpvar(num_seg, T); % storage discharge
    Dch = binvar(1, T); % storage charge state
    Dseg=binvar(num_seg-1, T); % storage segment
    SOC = sdpvar(num_seg, T); % storage SOC
    Pg = sdpvar(num_g, T); % generator power
    rg=sdpvar(num_g, T); % generator power
    lu=sdpvar(1, T); % generator power
    ld=sdpvar(1, T); % generator power


    %% constraints
    cons = [];
    % power balance-energy price dual
    for j = 1:T    
        cons = [cons, sum(Pg(:, j), 1) + sum(Pdis(:, j), 1) - sum(Pch(:, j), 1)+ld(:,j) == load_scenario(s, j)+lu(:,j)];
    end
    cons = [cons,lu>=0];
    cons = [cons,ld>=0];

   % storage SoC bounds
    for i = 1:num_seg
        for j = 1:T
            if j==1
            cons = [cons, SOC(i, j) == SOC0 + ((Pch(i, j)) .* yita + (-Pdis(i, j)) ./ yita) / ESenergycap];
            else
            cons = [cons, SOC(i, j) == SOC(i, j-1) + ((Pch(i, j)) .* yita + (-Pdis(i, j)) ./ yita) / ESenergycap];
            end
        end
    end
  
    % storage power bounds
    for j = 1:T    
        cons = [cons, Pch(:, j) >= 0];
        cons = [cons, Pdis(:, j) >= 0];
        cons = [cons, Pch(:, j) <= ESpowercap];
        cons = [cons, Pdis(:, j) <= ESpowercap];
        cons = [cons, sum(Pch(:, j)./ESpowercap,1) <= Dch(1, j)];
        cons = [cons, sum(Pdis(:, j)./ESpowercap,1) <=1-Dch(1, j)];
    end


    for j=1:T
    
        for i=2:num_seg
    
        cons = [cons, SOC(i, j) <=delta_E*Dseg(i-1,j)];    
    
        end
    
        for i=1:num_seg-1
    
        cons = [cons, SOC(i, j) >=delta_E*Dseg(i,j)];    
    
        end
    
       cons = [cons, SOC(num_seg, j)>=0];
       cons = [cons, SOC(1, j)<=delta_E];
       cons=[cons, sum(SOC(:,j),1)<=1];
    
    end


       for j=1:T

           cons = [cons, SOC(:, j)>=0];

       end


    % generator power bounds

for j=1:T    
cons=[cons,G_min<=Pg(:,j)<=G_max-rg(:,j)];
end


% generator reserve
for j=1:T    
cons=[cons,0<=rg(:,j)];
end


for j=1:T    
cons=[cons,sum(rg(:,j),1)>=0.1*load_c(:,j)];
end




    % generator RAMPUP and RAMPDOWN
    for j = 2:T    
        cons = [cons, Pg(:, j) - Pg(:, j-1) <= RU];
        cons = [cons, -RD <= Pg(:, j) - Pg(:, j-1)];
    end


    %% objective
    obj = 0;
    % quadratic + marginal cost
    for i = 1:num_g
        obj = obj + C_0(i, 1) + C_1(i, 1) * Pg(i, :) + C_2(i, 1) * (Pg(i, :) .* Pg(i, :));
    end
    obj = sum(obj, 2);

    for j = 1:T
        obj = obj + sum(bd2(:,j).*Pdis(:,j)-bc2(:,j).*Pch(:,j),1);
    end

    for j=1:T
        obj = obj + 2000*ld(:,j);
    end


    % optimization
    options = sdpsettings('verbose', 1, 'solver', 'gurobi');
    results = optimize(cons, obj, options);

    Dch_value=value(Dch);
    Dseg_value=value(Dseg);

   %% dual constraints
    cons = [];
    % power balance-energy price dual
    for j = 1:T    
        cons = [cons, sum(Pg(:, j), 1) + sum(Pdis(:, j), 1) - sum(Pch(:, j), 1)+ld(:,j) == load_scenario(s, j)+lu(:,j)];
    end
    cons = [cons,lu>=0];
    cons = [cons,ld>=0];



    for i = 1:num_seg
        for j = 1:T
            if j==1
            cons = [cons, SOC(i, j) == SOC0 + ((Pch(i, j)) .* yita + (-Pdis(i, j)) ./ yita) / ESenergycap];
            else
            cons = [cons, SOC(i, j) == SOC(i, j-1) + ((Pch(i, j)) .* yita + (-Pdis(i, j)) ./ yita) / ESenergycap];
            end
        end
    end
  
    % storage power bounds
    for j = 1:T    
        cons = [cons, Pch(:, j) >= 0];
        cons = [cons, Pdis(:, j) >= 0];
        cons = [cons, Pch(:, j) <= ESpowercap];
        cons = [cons, Pdis(:, j) <= ESpowercap];
        cons = [cons, sum(Pch(:, j)./ESpowercap,1) <= Dch_value(1, j)];
        cons = [cons, sum(Pdis(:, j)./ESpowercap,1) <=1-Dch_value(1, j)];
    end


    for j=1:T
    
        for i=2:num_seg
    
        cons = [cons, SOC(i, j) <=delta_E*Dseg_value(i-1,j)];    
    
        end
    
        for i=1:num_seg-1
    
        cons = [cons, SOC(i, j) >=delta_E*Dseg_value(i,j)];    
    
        end
    
       cons = [cons, SOC(num_seg, j)>=0];
       cons = [cons, SOC(1, j)<=delta_E];
       cons=[cons, sum(SOC(:,j),1)<=1];
    
    end


       for j=1:T

           cons = [cons, SOC(:, j)>=0];

       end




for j=1:T    
cons=[cons,G_min<=Pg(:,j)<=G_max-rg(:,j)];
end


% generator reserve
for j=1:T    
cons=[cons,0<=rg(:,j)];
end


for j=1:T    
cons=[cons,sum(rg(:,j),1)>=0.1*load_c(:,j)];
end


    % generator RAMPUP and RAMPDOWN
    for j = 2:T    
        cons = [cons, Pg(:, j) - Pg(:, j-1) <= RU];
        cons = [cons, -RD <= Pg(:, j) - Pg(:, j-1)];
    end


    %% objective
    obj = 0;
    % quadratic generation
    for i = 1:num_g
        obj = obj + C_0(i, 1) + C_1(i, 1) * Pg(i, :) + C_2(i, 1) * (Pg(i, :) .* Pg(i, :));
    end
    obj = sum(obj, 2);


    % storage bids
    for j = 1:T
        obj = obj + sum(bd2(:,j).*Pdis(:,j)-bc2(:,j).*Pch(:,j),1);
    end

    for j=1:T
        obj = obj + 2000*ld(:,j);
    end

    % optimization
    options = sdpsettings('verbose', 1, 'solver', 'gurobi');
    results = optimize(cons, obj, options);


    energy_price=abs(dual(cons(1:T)));

    % store result in a local variable
    temp_cost1 = 0;
    temp_cost2=0;
    

    for i = 1:num_g
        temp_cost1 = temp_cost1 + C_0(i, 1) + C_1(i, 1) * value(Pg(i, :)) + C_2(i, 1) * (value(Pg(i, :)) .* value(Pg(i, :)));
    end

    temp_cost1  = sum(temp_cost1 , 2);

    temp_cost1 =temp_cost1 +sum(MS*(value(sum(Pdis,1))));

    temp_cost2 = sum(energy_price'.*sum((value(Pdis-Pch)),1)-MS*(value(sum(Pdis,1))));

    temp_cost3=sum(energy_price'.*load_scenario(s,:));


       
    system_cost2(s, 1) = temp_cost1;
    storage_profit2(s, 1) = temp_cost2;
    consumer_payment2(s, 1)=temp_cost3;

end

% [a,b]=max(system_cost1-system_cost2);
% 
% 
% % data = [system_cost1-system_cost2];    % 示例数据
% % sortedData = sort(data);  % 对数据进行排序
% % n = length(sortedData);
% % b = ceil(0.95 * n);     % 计算位置，下标从 1 开始
% % q10 = sortedData(b);    % 得到 25% 分位数对应的值
% % b=find(data==q10);
% 
% 
% system_cost1_av=system_cost1(b);
% system_cost2_av=system_cost2(b);
% 
% storage_profit1_av=storage_profit1(b);
% storage_profit2_av=storage_profit2(b);
% 
% consumer_payment1_av=consumer_payment1(b);
% consumer_payment2_av=consumer_payment2(b);

system_cost1_av=mean(system_cost1(system_cost1~=0));
system_cost2_av=mean(system_cost2(system_cost2~=0));

storage_profit1_av=mean(storage_profit1(storage_profit1~=0));
storage_profit2_av=mean(storage_profit2(storage_profit2~=0));

consumer_payment1_av=mean(consumer_payment1(consumer_payment1~=0));
consumer_payment2_av=mean(consumer_payment2(consumer_payment2~=0));





