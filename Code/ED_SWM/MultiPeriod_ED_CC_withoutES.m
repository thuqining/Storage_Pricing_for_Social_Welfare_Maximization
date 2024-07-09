%% initial test
%@Create Time    : 2024/3/15
%@Author  : Ning Qi nq2176@columbia.edu
%@File    : MultiPeriod_ED_CC_withoutES.m
%@Description : This is the main code for multi-period chance-constrained
%economic dispatch for system operator and can generate both dispatch and
%price for storage and generator
% storage: no ES

function [obj_value,energy_price,reserve_price,rho_g_value,Pg_value]=MultiPeriod_ED_CC_withoutES(a,b,c,d,e,f,g)
%% inputs
%generation cost function
%net load

load data.mat
T=24; %time duration
%% constants
wind_ratio=a;% wind ratio
netload_error=load_error-max(max(load_real))*(wind_real-wind_forecast)*wind_ratio;

netload_error=reshape(netload_error',96,365)';% netload error, [-, +]

% 初始化一个新的 365x24 的数组 B
B = zeros(365, 24);

for i = 1:365
    
    reshapedRow = reshape(netload_error(i, :), 4, T);
    B(i, :) = mean(reshapedRow, 1);
end

netload_error=B;

epsilon=b; % probability level
for i=1:T
% [mu(i),sigma(i)] = normfit(netload_error(:,i),0.05); % μ, σ time different

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

sigma=sigma*c;



load_c=netload(d,:);

load_d1=zeros(1,T);
load_d2=zeros(1,T);

% alpha_value=zeros(1,T);
% beta_value=zeros(1,T);
% gama_value=zeros(1,T);
% for i=1:T
% [alpha_value(i),beta_value(i),gama_value(i)]=fit_versatile(B(:,i));
% end

for i=1:T

load_d1(:,i)= norminv(epsilon, mu(i), sigma(i));
% if epsilon<1/6
% load_d1(:,i)= mu(i)-sqrt(2/(9*epsilon))*sigma(i);
% else 
%     load_d1(:,i)= mu(i)-sqrt(3)*(1-2*epsilon)*sigma(i);
% end
% 
%load_d1(:,i)= mu(i)-sqrt((1-epsilon)/(epsilon))*sigma(i);
%load_d1(:,i)= quantile(B(:,i), epsilon);
% load_d1(:,i)=gama_value(i)-imag(log((epsilon)^(-1/beta_value(i))-1)./alpha_value(i));

end

for i=1:T

load_d2(:,i)= norminv(1-epsilon, mu(i), sigma(i));

% if epsilon<1/6
% load_d2(:,i)= mu(i)+sqrt(2/(9*epsilon))*sigma(i);
% else 
%     load_d2(:,i)= mu(i)+sqrt(3)*(1-2*epsilon)*sigma(i);
% end

%load_d2(:,i)= mu(i)+sqrt((1-epsilon)/(epsilon))*sigma(i);
%load_d2(:,i)= quantile(B(:,i), 1-epsilon);
% load_d2(:,i)=gama_value(i)-imag(log((1-epsilon)^(-1/beta_value(i))-1)./alpha_value(i));
end


num_e=1; %number of storage
num_g=size(data_generator,1); %number of generator
%num_g=1; %number of generator
yita=0.95; % storage efficiency
SOC0=g; % initial SoC
ESpowercap=e*max(max(netload)); % storage power cap
ESenergycap=ESpowercap*f; % storage energy cap
G_min=data_generator(:,2)*0;
G_max=data_generator(:,1)*1;
%G_min=zeros(num_g,1);
% G_min=0;%% assume one generator
% G_max=sum(G_max);
% C_0=0*ones(num_g,1);
% C_1=43*ones(num_g,1);
% C_2=0.00292*ones(num_g,1)*0;
% C_0=10707*ones(num_g,1);
% C_1=153.16*ones(num_g,1);
% C_2=0.0346*ones(num_g,1);

C_0=data_generator(:,8);
C_1=data_generator(:,9);
C_2=data_generator(:,10);
% C_2=zeros(num_g,1);
% load generator_fit.mat
% C_0=bQuad(:,3);
% C_1=bQuad(:,2);
% C_2=bQuad(:,1);
%% decision variables

Pg=sdpvar(num_g,T); % generator power
rho_g=sdpvar(num_g,T); % generator reserve ratio
%% constraints
cons=[];

for j=1:T    
cons=[cons, sum(Pg(:,j),1)==load_c(:,j)];
end

for j=1:T    
cons=[cons,sum(rho_g(:,j),1)==1];
end

%cons=[cons,SOC(:,1)==SOC(:,T+1)];

for j=1:T   
cons=[cons,0<=rho_g(:,j)<=1];
end

for j=1:T    
cons=[cons,G_min<=Pg(:,j)+rho_g(:,j)*load_d1(1,j)];
end

for j=1:T    
cons=[cons,Pg(:,j)+rho_g(:,j)*load_d2(1,j)<=G_max];
end



%% objective
obj=0;

% % quadratic
for i=1:num_g
obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:)+rho_g(i,:).*mu)+C_2(i,1)*(Pg(i,:).*Pg(i,:)+2*Pg(i,:).*rho_g(i,:).*mu+rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma));
end

% for i=1:num_g
% obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:))+C_2(i)*(Pg(i,:).*Pg(i,:));
% end
% cubic
% for i=1:num_g
% obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:)+rho_g(i,:).*mu)+C_2(i)*(Pg(i,:).*Pg(i,:)+2*Pg(i,:).*rho_g(i,:).*mu+rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma))+C_3(i)*(Pg(i,:).*Pg(i,:).*Pg(i,:)+3*Pg(i,:).*Pg(i,:).*rho_g(i,:).*mu+3*Pg(i,:).*rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma));
% end
% % quartic
% for i=1:num_g
% obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:)+rho_g(i,:).*mu)+C_2(i)*(Pg(i,:).*Pg(i,:)+2*Pg(i,:).*rho_g(i,:).*mu+rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma))+C_3(i)*(Pg(i,:).*Pg(i,:).*Pg(i,:)+3*Pg(i,:).*Pg(i,:).*rho_g(i,:).*mu+3*Pg(i,:).*rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma))+C_4(i)*(Pg(i,:).*Pg(i,:).*Pg(i,:).*Pg(i,:)+4*Pg(i,:).*Pg(i,:).*Pg(i,:).*rho_g(i,:).*mu+6*Pg(i,:).*Pg(i,:).*rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma)+rho_g(i,:).*rho_g(i,:).*rho_g(i,:).*rho_g(i,:).*(3*sigma.*sigma.*sigma.*sigma+6*sigma.*sigma.*mu.*mu+mu.*mu.*mu.*mu));
% end

obj=sum(obj,2);
options = sdpsettings('verbose', 1, 'solver', 'gurobi');
results= optimize(cons,obj,options);

Pg_value=value(Pg);
Pg_value=sum(Pg_value,1);
rho_g_value=value(rho_g);

obj_value=value(obj);
% 
energy_price=abs(dual(cons(1:T)));
reserve_price=dual(cons(T+1:2*T));
