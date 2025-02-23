%@Create Time    : 2024/3/15
%@Author  : Ning Qi nq2176@columbia.edu
%@File    : MultiPeriod_ED_CC.m
%@Description : This is the main code for multi-period chance-constrained
%economic dispatch for system operator and can generate both dispatch and
%price for storage and generator
% storage: energy and reserve market

function [system_cost1,storage_profit1,consumer_payment1,system_cost2,storage_profit2,consumer_payment2,energy_price1,energy_price2,SOC_value1,SOC_value2]=RT_price_single(a,b,c,d,e,f,g1,g2,h,bd1,bc1,bd2,bc2,clear_time)

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



mu=mean(netload_error(:,clear_time));
sigma=std(netload_error(:,clear_time));


wind_forecast=reshape(wind_forecast',96,365)';
netload=load_forecast-max(max(load_real))*wind_forecast*wind_ratio; % netload forcast
C = zeros(365, T);
for i = 1:365
    reshapedRow = reshape(netload(i, :), 4, T);
    C(i, :) = mean(reshapedRow, 1);
end
netload=C;


sigma=sigma*c;% change the STD of forecast error for net load

load_c=netload(d,clear_time);

load_scenario=zeros(b,1);

rng(30);
stream = RandStream('mt19937ar','Seed',30); % 明确指定流

for i=1:b

load_scenario(i,:)=normrnd(mu,sigma);

end

load_scenario=load_scenario+ones(b,1)*load_c;


num_seg=size(bd1,1); %number of segment
num_g=size(data_generator,1); %number of generator
yita=0.95; % storage efficiency
% num_e=2; %number of storage
% yita=[0.8  0.95 ]'; % storage efficiency
SOC0_1=g1; % initial SoC
SOC0_2=g2; % initial SoC
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
energy_price1=zeros(b, 1);
SOC_value1=zeros(b, num_seg);


system_cost2 = zeros(b, 1);
storage_profit2= zeros(b, 1);
consumer_payment2=zeros(b, 1);
energy_price2=zeros(b, 1);
SOC_value2=zeros(b, num_seg);

% system_cost3 = zeros(b, 1);
% storage_profit3= zeros(b, 1);

%% DP bids
parfor s = 1:b
    %% decision variables
    Pch = sdpvar(num_seg, 1); % storage charge
    Pdis = sdpvar(num_seg, 1); % storage discharge
    Dch = binvar(1, 1); % storage charge state
    Dseg=binvar(num_seg-1, 1); % storage segment
    SOC = sdpvar(num_seg, 1); % storage SOC
    Pg = sdpvar(num_g, 1); % generator power
    %rg = sdpvar(num_g,1); % generator reserve
    lu=sdpvar(1, 1); 
    ld=sdpvar(1, 1); 


    %% constraints
    cons = [];
    % power balance-energy price dual
    cons = [cons, sum(Pg(:, 1), 1) + sum(Pdis(:, 1), 1) - sum(Pch(:, 1), 1) == load_scenario(s, 1)+lu(:,1)-ld(:,1)];

    cons=[cons,lu>=0];
    cons=[cons,ld>=0];

   % storage SoC bounds
    for i = 1:num_seg
            cons = [cons, SOC(i, 1) == SOC0_1(s,i) + ((Pch(i, 1)) .* yita + (-Pdis(i, 1)) ./ yita) / ESenergycap];
    end
  
    % storage power bounds
        cons = [cons, Pch(:, 1) >= 0];
        cons = [cons, Pdis(:, 1) >= 0];
        cons = [cons, Pch(:, 1) <= ESpowercap];
        cons = [cons, Pdis(:, 1) <= ESpowercap];

        binary_constraints = [];

        binary_constraints = [binary_constraints, sum(Pch(:, 1)./ESpowercap,1) <= Dch(1, 1)];
        binary_constraints = [binary_constraints, sum(Pdis(:, 1)./ESpowercap,1) <=1-Dch(1, 1)];
        for i=2:num_seg
        binary_constraints = [binary_constraints, SOC(i, 1) <=delta_E*Dseg(i-1,1)];    
        end
        for i=1:num_seg-1
        binary_constraints = [binary_constraints, SOC(i, 1) >=delta_E*Dseg(i,1)];    
        end
        

% SoC dependent bid

        cons = [cons, SOC(num_seg, 1)>=0];
        cons = [cons, SOC(1, 1)<=delta_E];
        cons=[cons, sum(SOC(:,1),1)<=1];
        cons = [cons, SOC(:, 1)>=0];

    % generator power bounds
  

   
    cons = [cons, G_min <= Pg(:,1) <= G_max];  
    % cons = [cons, 0 <= rg(:,1)];    
    % cons = [cons, sum(rg(:,1),1) >= 0.1 * load_c(:,1)];



        cons = [cons, binary_constraints];


    %% objective
    obj = 0;
    % quadratic + marginal cost
    for i = 1:num_g
        obj = obj + C_0(i, 1) + C_1(i, 1) * Pg(i, :) + C_2(i, 1) * (Pg(i, :) .* Pg(i, :));
        %obj = obj + C_1(i, 1) * Pg(i, :) + C_2(i, 1) * (Pg(i, :) .* Pg(i, :));
    end

    obj = obj + sum(bd1(:,1).*Pdis(:,1)-bc1(:,1).*Pch(:,1),1);

    obj=obj+2000*sum(ld);% loss of load


    % optimization
    options = sdpsettings('verbose', 1, 'solver', 'gurobi');
    results = optimize(cons, obj, options);

    Dch_value=value(Dch);
    Dseg_value=value(Dseg);


    % remove binary_constraints
cons(end - length(binary_constraints) + 1:end) = [];

% **替换二进制变量相关的约束**
new_binary_constraints = [];


% **加入新的约束**
cons = [cons, new_binary_constraints];

new_binary_constraints= [new_binary_constraints, sum(Pch(:, 1)./ESpowercap,1) <= Dch_value(1, 1)];
new_binary_constraints= [new_binary_constraints, sum(Pdis(:, 1)./ESpowercap,1) <=1-Dch_value(1, 1)];
for i=2:num_seg
new_binary_constraints = [new_binary_constraints, SOC(i, 1) <=delta_E*Dseg_value(i-1,1)];    
end
for i=1:num_seg-1
new_binary_constraints = [new_binary_constraints, SOC(i, 1) >=delta_E*Dseg_value(i,1)];    
end

cons = [cons, new_binary_constraints];

%% **求解 LP（计算对偶变量）**
  
 results = optimize(cons, obj, options);

 energy_price=abs(dual(cons(1)));

    % store result in a local variable
    temp_cost1 = 0;
    temp_cost2=0;
    temp_cost3=0;
    

    for i = 1:num_g
        %temp_cost1 = temp_cost1 + C_1(i, 1) * value(Pg(i, :)) + C_2(i, 1) * (value(Pg(i, :)) .* value(Pg(i, :)));
        temp_cost1 = temp_cost1 + C_0(i, 1) + C_1(i, 1) * value(Pg(i, :)) + C_2(i, 1) * (value(Pg(i, :)) .* value(Pg(i, :)));
    end

    temp_cost1  = sum(temp_cost1 , 2);
    

    temp_cost1 =temp_cost1 +sum(MS*(value(sum(Pdis,1))));

    %temp_cost2 = sum(energy_price'.*sum((value(Pdis-Pch)),1));
    temp_cost2 = sum(energy_price'.*sum((value(Pdis-Pch)),1)-MS*(value(sum(Pdis,1))));
    temp_cost3=sum(energy_price'.*load_scenario(s,:));


       
    system_cost1(s, 1) = temp_cost1;
    storage_profit1(s, 1) = temp_cost2;
    consumer_payment1(s, 1)=temp_cost3;
    energy_price1(s,:)=energy_price';
    SOC_value1(s,:)=value(SOC);


end

% %% capped

parfor s = 1:b
    %% decision variables
    Pch = sdpvar(num_seg, 1); % storage charge
    Pdis = sdpvar(num_seg, 1); % storage discharge
    Dch = binvar(1, 1); % storage charge state
    Dseg=binvar(num_seg-1, 1); % storage segment
    SOC = sdpvar(num_seg, 1); % storage SOC
    Pg = sdpvar(num_g, 1); % generator power
    %rg = sdpvar(num_g,1); % generator reserve
    lu=sdpvar(1, 1); 
    ld=sdpvar(1, 1); 


    %% constraints
    cons = [];
    % power balance-energy price dual
    cons = [cons, sum(Pg(:, 1), 1) + sum(Pdis(:, 1), 1) - sum(Pch(:, 1), 1) == load_scenario(s, 1)+lu(:,1)-ld(:,1)];

    cons=[cons,lu>=0];
    cons=[cons,ld>=0];

   % storage SoC bounds
    for i = 1:num_seg
            cons = [cons, SOC(i, 1) == SOC0_2(s,i) + ((Pch(i, 1)) .* yita + (-Pdis(i, 1)) ./ yita) / ESenergycap];
    end
  
    % storage power bounds
        cons = [cons, Pch(:, 1) >= 0];
        cons = [cons, Pdis(:, 1) >= 0];
        cons = [cons, Pch(:, 1) <= ESpowercap];
        cons = [cons, Pdis(:, 1) <= ESpowercap];

        binary_constraints = [];

        binary_constraints = [binary_constraints, sum(Pch(:, 1)./ESpowercap,1) <= Dch(1, 1)];
        binary_constraints = [binary_constraints, sum(Pdis(:, 1)./ESpowercap,1) <=1-Dch(1, 1)];
        for i=2:num_seg
        binary_constraints = [binary_constraints, SOC(i, 1) <=delta_E*Dseg(i-1,1)];    
        end
        for i=1:num_seg-1
        binary_constraints = [binary_constraints, SOC(i, 1) >=delta_E*Dseg(i,1)];    
        end
        

% SoC dependent bid

        cons = [cons, SOC(num_seg, 1)>=0];
        cons = [cons, SOC(1, 1)<=delta_E];
        cons=[cons, sum(SOC(:,1),1)<=1];
        cons = [cons, SOC(:, 1)>=0];

    % generator power bounds
  
    cons = [cons, G_min <= Pg(:,1) <= G_max];  
    % cons = [cons, 0 <= rg(:,1)];    
    % cons = [cons, sum(rg(:,1),1) >= 0.1 * load_c(:,1)];



        cons = [cons, binary_constraints];


    %% objective
    obj = 0;
    % quadratic + marginal cost
    for i = 1:num_g
        obj = obj + C_0(i, 1) + C_1(i, 1) * Pg(i, :) + C_2(i, 1) * (Pg(i, :) .* Pg(i, :));
        %obj = obj + C_1(i, 1) * Pg(i, :) + C_2(i, 1) * (Pg(i, :) .* Pg(i, :));
    end

    obj = obj + sum(bd2(:,1).*Pdis(:,1)-bc2(:,1).*Pch(:,1),1);

    obj=obj+2000*sum(ld);% loss of load


    % optimization
    options = sdpsettings('verbose', 1, 'solver', 'gurobi');
    results = optimize(cons, obj, options);

    Dch_value=value(Dch);
    Dseg_value=value(Dseg);


    % remove binary_constraints
cons(end - length(binary_constraints) + 1:end) = [];

% **替换二进制变量相关的约束**
new_binary_constraints = [];


% **加入新的约束**
cons = [cons, new_binary_constraints];

new_binary_constraints= [new_binary_constraints, sum(Pch(:, 1)./ESpowercap,1) <= Dch_value(1, 1)];
new_binary_constraints= [new_binary_constraints, sum(Pdis(:, 1)./ESpowercap,1) <=1-Dch_value(1, 1)];
for i=2:num_seg
new_binary_constraints = [new_binary_constraints, SOC(i, 1) <=delta_E*Dseg_value(i-1,1)];    
end
for i=1:num_seg-1
new_binary_constraints = [new_binary_constraints, SOC(i, 1) >=delta_E*Dseg_value(i,1)];    
end

cons = [cons, new_binary_constraints];

%% **求解 LP（计算对偶变量）**
  
 results = optimize(cons, obj, options);

 energy_price=abs(dual(cons(1)));

    % store result in a local variable
    temp_cost1 = 0;
    temp_cost2=0;
    temp_cost3=0;
    

    for i = 1:num_g
        %temp_cost1 = temp_cost1 + C_1(i, 1) * value(Pg(i, :)) + C_2(i, 1) * (value(Pg(i, :)) .* value(Pg(i, :)));
        temp_cost1 = temp_cost1 + C_0(i, 1) + C_1(i, 1) * value(Pg(i, :)) + C_2(i, 1) * (value(Pg(i, :)) .* value(Pg(i, :)));
    end

    temp_cost1  = sum(temp_cost1 , 2);
    

    temp_cost1 =temp_cost1 +sum(MS*(value(sum(Pdis,1))));
     %temp_cost2 = sum(energy_price'.*sum((value(Pdis-Pch)),1));
    temp_cost2 = sum(energy_price'.*sum((value(Pdis-Pch)),1)-MS*(value(sum(Pdis,1))));
    temp_cost3=sum(energy_price'.*load_scenario(s,:));


       
    system_cost2(s, 1) = temp_cost1;
    storage_profit2(s, 1) = temp_cost2;
    consumer_payment2(s, 1)=temp_cost3;
    energy_price2(s,:)=energy_price';
    SOC_value2(s,:)=value(SOC);


end


% system_cost1_av=mean(system_cost1);
% system_cost2_av=mean(system_cost2);
% 
% storage_profit1_av=mean(storage_profit1);
% storage_profit2_av=mean(storage_profit2);
% 
% consumer_payment1_av=mean(consumer_payment1);
% consumer_payment2_av=mean(consumer_payment2);
% 
% price_std1=std(energy_price1);
% price_std2=std(energy_price2);









