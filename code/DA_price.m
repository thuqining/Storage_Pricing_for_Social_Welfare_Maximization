%@Create Time    : 2024/3/15
%@Author  : Ning Qi nq2176@columbia.edu
%@File    : MultiPeriod_ED_CC.m
%@Description : This is the main code for multi-period chance-constrained
%economic dispatch for system operator and can generate both dispatch and
%price for storage and generator
% storage: energy and reserve market

function [energy_price]=DA_price(a,b,c,d,e,f,g,h)

%% inputs
% a=0.3;
% b=0.25;
% c=1; %sigma, original forecast error is too small
% d=10; %day number
% e=0.2; %storage capacity ratio (load capacity)
% f=4; %storage duration
% g=1; %initial SoC
% h=10; %storage marginal cost
% l=0; %mu=0,shifting load forecast error

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


num_e=1; %number of storage
num_g=size(data_generator,1); %number of generator
yita=0.95; % storage efficiency
% num_e=2; %number of storage
% yita=[0.8  0.95 ]'; % storage efficiency
SOC0=g; % initial SoC
ESpowercap=e*max(max(netload))/num_e; % storage power cap
ESenergycap=ESpowercap*f; % storage energy cap
G_min=data_generator(:,2)*0;% minimum capacity=0
G_max=data_generator(:,1)*1;
RU=data_generator(:,6);
RD=data_generator(:,6);
C_0=data_generator(:,8);
C_1=data_generator(:,9);
C_2=data_generator(:,10);

%% decision variables
%% decision variables
Pch=sdpvar(num_e,T); % storage charge
Pdis=sdpvar(num_e,T); % storage discharge
Dch=binvar(num_e,T); % storage charge state
Ddis=binvar(num_e,T); % storage discharge state
SOC=sdpvar(num_e,T+1); % storage SOC
Pg=sdpvar(num_g,T); % generator power
rg=sdpvar(num_g,T); % generator power

%% constraints
cons=[];
% power balance-energy price dual
for j=1:T    
cons=[cons, sum(Pg(:,j),1)+sum(Pdis(:,j),1)-sum(Pch(:,j),1)==load_c(:,j)];
end
% SoC recovery

cons=[cons,SOC(:,1)==SOC0];
for i=1:num_e
for j=2:T+1    
cons=[cons, SOC(i,j)==SOC(i,j-1)+((Pch(i,j-1)).*yita(i)+(-Pdis(i,j-1))./yita(i))/ESenergycap];
end
end


% storage power bounds
for j=1:T    
cons=[cons,Pch(:,j)>=0];
end

for j=1:T    
cons=[cons,Pdis(:,j)>=0];
end


for j=1:T    
cons=[cons,Pch(:,j)<=Dch(:,j)*ESpowercap];
end
% 
for j=1:T    
cons=[cons,Pdis(:,j)<=Ddis(:,j)*ESpowercap];
end


% complementary charging-discharging
for j=1:T    
cons=[cons,0<=Dch(:,j)+Ddis(:,j)<=1];
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



% generator RAMPUP
for j=2:T    
cons=[cons,Pg(:,j)-Pg(:,j-1)<=RU];
end

% generator RAMPDOWN
for j=2:T    
cons=[cons,-RD<=Pg(:,j)-Pg(:,j-1)];
end

% storage SoC bounds
for j=2:T+1    
cons=[cons, 0<=SOC(:,j)<=1];
end


%% objective
obj=0;

% % quadratic+marginal cost
for i=1:num_g
obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:))+C_2(i,1)*(Pg(i,:).*Pg(i,:));
end
obj=sum(obj,2);
for j=1:T
obj=obj+sum((Pdis(:,j))*MS,1);
end

options = sdpsettings('verbose', 1, 'solver', 'gurobi');
results= optimize(cons,obj,options);

Dch_value=value(Dch);
Ddis_value=value(Ddis);

%% using integer value to solve it again so that we can derive dual

cons=[];
% power balance-energy price dual
for j=1:T    
cons=[cons, sum(Pg(:,j),1)+sum(Pdis(:,j),1)-sum(Pch(:,j),1)==load_c(:,j)];
end
% SoC recovery

cons=[cons,SOC(:,1)==SOC0];
for i=1:num_e
for j=2:T+1    
cons=[cons, SOC(i,j)==SOC(i,j-1)+((Pch(i,j-1)).*yita(i)+(-Pdis(i,j-1))./yita(i))/ESenergycap];
end
end


% storage power bounds
for j=1:T    
cons=[cons,Pch(:,j)>=0];
end

for j=1:T    
cons=[cons,Pdis(:,j)>=0];
end


for j=1:T    
cons=[cons,Pch(:,j)<=Dch_value(:,j)*ESpowercap];
end
% 
for j=1:T    
cons=[cons,Pdis(:,j)<=Ddis_value(:,j)*ESpowercap];
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



% generator RAMPUP
for j=2:T    
cons=[cons,Pg(:,j)-Pg(:,j-1)<=RU];
end

% generator RAMPDOWN
for j=2:T    
cons=[cons,-RD<=Pg(:,j)-Pg(:,j-1)];
end

% storage SoC bounds
for j=2:T+1    
cons=[cons, 0<=SOC(:,j)<=1];
end


%% objective
obj=0;

% % quadratic+marginal cost
for i=1:num_g
obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:))+C_2(i,1)*(Pg(i,:).*Pg(i,:));
end
obj=sum(obj,2);
for j=1:T
obj=obj+sum((Pdis(:,j))*MS,1);
end

options = sdpsettings('verbose', 1, 'solver', 'gurobi');
results= optimize(cons,obj,options);


energy_price=abs(dual(cons(1:T)));
