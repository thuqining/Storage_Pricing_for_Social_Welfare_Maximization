%@Create Time    : 2024/3/15
%@Author  : Ning Qi nq2176@columbia.edu
%@File    : MultiPeriod_ED_CC.m
%@Description : This is the main code for multi-period chance-constrained
%economic dispatch for system operator and can generate both dispatch and
%price for storage and generator
% storage: energy and reserve market

function [obj_value,energy_price,storage_price,reserve_price,SOC_value,Pdis_value,Pch_value,Dch_value,Ddis_value]=MultiPeriod_ED_CC(a,b,c,d,e,f,g)
%% inputs
%generation cost function
%historical net load
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
load_d1=zeros(1,T);% lower quantile
load_d2=zeros(1,T);% higher quantile

% learning versatile distribution
% alpha_value=zeros(1,T);
% beta_value=zeros(1,T);
% gama_value=zeros(1,T);
% for i=1:T
% [alpha_value(i),beta_value(i),gama_value(i)]=fit_versatile(B(:,i));
% end

for i=1:T
%gaussion distribution
load_d1(:,i)= norminv(epsilon, mu(i), sigma(i));
%robust approximation1
% if epsilon<1/6
% load_d1(:,i)= mu(i)-sqrt(2/(9*epsilon))*sigma(i);
% else 
%     load_d1(:,i)= mu(i)-sqrt(3)*(1-2*epsilon)*sigma(i);
% end
%robust approximation2 
%load_d1(:,i)= mu(i)-sqrt((1-epsilon)/(epsilon))*sigma(i);
%Empirical
%load_d1(:,i)= quantile(B(:,i), epsilon);%
%versatile
% load_d1(:,i)=gama_value(i)-imag(log((epsilon)^(-1/beta_value(i))-1)./alpha_value(i));
end

for i=1:T
%gaussion distribution    
load_d2(:,i)= norminv(1-epsilon, mu(i), sigma(i));
%robust approximation1
% if epsilon<1/6
% load_d2(:,i)= mu(i)+sqrt(2/(9*epsilon))*sigma(i);
% else 
%     load_d2(:,i)= mu(i)+sqrt(3)*(1-2*epsilon)*sigma(i);
% end
%robust approximation2 
%load_d2(:,i)= mu(i)+sqrt((1-epsilon)/(epsilon))*sigma(i);
%Empirical
%load_d2(:,i)= quantile(B(:,i), 1-epsilon);
%versatile
% load_d2(:,i)=gama_value(i)-imag(log((1-epsilon)^(-1/beta_value(i))-1)./alpha_value(i));
end


num_e=1; %number of storage
num_g=size(data_generator,1); %number of generator
%num_g=1; %number of generator
yita=0.95; % storage efficiency
SOC0=g; % initial SoC
ESpowercap=e*max(max(netload)); % storage power cap
ESenergycap=ESpowercap*f; % storage energy cap
G_min=data_generator(:,2)*0;% minimum capacity=0
G_max=data_generator(:,1)*1;
%G_min=0;%% assume one generator
% G_max=sum(G_max);
C_0=data_generator(:,8);
C_1=data_generator(:,9);
C_2=data_generator(:,10);
%C_2=zeros(num_g,1);
% load generator_fit.mat
% C_0=bQuad(:,3);
% C_1=bQuad(:,2);
% C_2=bQuad(:,1);
%% decision variables
Pch=sdpvar(num_e,T); % storage charge
Pdis=sdpvar(num_e,T); % storage discharge
Dch=binvar(num_e,T); % storage charge state
Ddis=binvar(num_e,T); % storage discharge state
SOC=sdpvar(num_e,T+1); % storage SOC
Pg=sdpvar(num_g,T); % generator power
rho_g=sdpvar(num_g,T); % generator reserve ratio
rho_e=sdpvar(num_e,T); % storage reserve ratio
%% constraints
cons=[];
% power balance-energy price dual
for j=1:T    
cons=[cons, sum(Pg(:,j),1)+sum(Pdis(:,j),1)-sum(Pch(:,j),1)==load_c(:,j)];
end
% SoC recovery
cons=[cons,SOC(:,1)==SOC0];
% SoC-P-opportunity price dual
for j=2:T+1    
cons=[cons, SOC(:,j)==SOC(:,j-1)+(-Pdis(:,j-1)/yita+Pch(:,j-1)*yita)/ESenergycap];
end
% reserve-reserve price dual
for j=1:T    
cons=[cons,sum(rho_g(:,j),1)+sum(rho_e(:,j),1)==1];
end

% generator reserve ratio
for j=1:T   
cons=[cons,0<=rho_g(:,j)<=1];
end

% storage reserve ratio
for j=1:T   
cons=[cons,0<=rho_e(:,j)<=1];
end

% storage power bounds
for j=1:T    
cons=[cons,0<=Pch(:,j)<=Dch(:,j)*ESpowercap,0<=Pdis(:,j)<=Ddis(:,j)*ESpowercap];
end

% storage power bounds
for j=1:T    
cons=[cons,-ESpowercap<=-Pch(:,j)+rho_e(:,j)*load_d1(1,j)];
end
% storage power bounds 
for j=1:T    
cons=[cons,Pdis(:,j)+rho_e(:,j)*load_d2(1,j)<=ESpowercap];
end

% complementary charging-discharging
for j=1:T    
cons=[cons,0<=Dch(:,j)+Ddis(:,j)<=1];
end

% generator power bounds
for j=1:T    
cons=[cons,G_min<=Pg(:,j)+rho_g(:,j)*load_d1(1,j)];
end

% generator power bounds
for j=1:T    
cons=[cons,Pg(:,j)+rho_g(:,j)*load_d2(1,j)<=G_max];
end

% storage SoC bounds
for j=1:T    
cons=[cons,SOC(:,j)<=1-(Pch(:,j)-rho_e(:,j)*load_d1(1,j))*yita/ESenergycap];
end

% storage SoC bounds
for j=1:T    
cons=[cons,(Pdis(:,j)+rho_e(:,j)*load_d2(1,j))/(yita*ESenergycap)<=SOC(:,j)];
end

%% objective
obj=0;

% % quadratic+marginal cost
for i=1:num_g
obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:)+rho_g(i,:).*mu)+C_2(i,1)*(Pg(i,:).*Pg(i,:)+2*Pg(i,:).*rho_g(i,:).*mu+rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma));
end
obj=obj+(Pdis+rho_e.*mu)*20;
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

% Dch_value=zeros(1,T);
% Ddis_value=zeros(1,T);
% Pdis_value=value(Pdis);
% Pch_value=value(Pch);
% Dch_value(Pch_value>=0.01)=1;
% Ddis_value(Pdis_value>=0.01)=1;
Dch_value=value(Dch);
Ddis_value=value(Ddis);
%% put integer variables into MIQP again
% C_0=bCubic(:,4);
% C_1=bCubic(:,3);
% C_2=bCubic(:,2);
% C_3=bCubic(:,1);

%% using integer value to solve it again so that we can generate dual
cons=[];

for j=1:T    
cons=[cons, sum(Pg(:,j),1)+sum(Pdis(:,j),1)-sum(Pch(:,j),1)==load_c(:,j)];
end

cons=[cons,SOC(:,1)==SOC0];

for j=2:T+1    
cons=[cons, SOC(:,j)-SOC(:,j-1)==(-Pdis(:,j-1)/yita+Pch(:,j-1)*yita)/ESenergycap];
end

for j=1:T    
cons=[cons,sum(rho_g(:,j),1)+sum(rho_e(:,j),1)==1];
end


for j=1:T   
cons=[cons,0<=rho_g(:,j)<=1];
end

for j=1:T   
cons=[cons,0<=rho_e(:,j)<=1];
end


for j=1:T    
cons=[cons,0<=Pch(:,j)<=Dch_value(:,j)*ESpowercap,0<=Pdis(:,j)<=Ddis_value(:,j)*ESpowercap];
end


for j=1:T    
cons=[cons,-ESpowercap<=-Pch(:,j)+rho_e(:,j)*load_d1(1,j)];
end
% 
for j=1:T    
cons=[cons,Pdis(:,j)+rho_e(:,j)*load_d2(1,j)<=ESpowercap];
end


for j=1:T    
cons=[cons,G_min<=Pg(:,j)+rho_g(:,j)*load_d1(1,j)];
end

for j=1:T    
cons=[cons,Pg(:,j)+rho_g(:,j)*load_d2(1,j)<=G_max];
end

for j=1:T    
cons=[cons,SOC(:,j)<=1-(Pch(:,j)-rho_e(:,j)*load_d1(1,j))*yita/ESenergycap];
end
% 
for j=1:T    
cons=[cons,(Pdis(:,j)+rho_e(:,j)*load_d2(1,j))/(yita*ESenergycap)<=SOC(:,j)];
end

%% objective
obj=0;

% quadratic
for i=1:num_g
obj=obj+C_0(i,1)+C_1(i,1)*(Pg(i,:)+rho_g(i,:).*mu)+C_2(i,1)*(Pg(i,:).*Pg(i,:)+2*Pg(i,:).*rho_g(i,:).*mu+rho_g(i,:).*rho_g(i,:).*(mu.*mu+sigma.*sigma));
end
obj=obj+(Pdis+rho_e.*mu)*20;
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

% options = sdpsettings('verbose', 1, 'solver', 'ipopt');
% options.ipopt.tol = 0.001;%
options = sdpsettings('verbose', 1, 'solver', 'gurobi');
results= optimize(cons,obj,options );

Pdis_value=value(Pdis);
Pch_value=value(Pch);
SOC_value=value(SOC)';
Pg_value=value(Pg);
Pg_value=sum(Pg_value,1);
rho_e_value=value(rho_e)';
rho_g_value=value(rho_g);

obj_value=value(obj);
% 
energy_price=abs(dual(cons(1:T)));
storage_price=abs(dual(cons(T+1:2*T+1))/ESenergycap);
reserve_price=dual(cons(2*T+2:3*T+1));
reserve_upbound_e=abs(dual(cons(end-2*T+1:end-T)));
reserve_downbound_e=abs(dual(cons(end-T+1:end)));
reserve_upbound_g=abs(dual(cons(end-3*T+1:end-2*T)));
reserve_downbound_g=abs(dual(cons(end-4*T+1:end-3*T)));
reserve_upbound_p=abs(dual(cons(end-5*T+1:end-4*T)));
reserve_downbound_p=abs(dual(cons(end-6*T+1:end-5*T)));
