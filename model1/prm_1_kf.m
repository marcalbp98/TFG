% Marçal Bravo (2019), TFG - Determination of a population rate model with
% Kalman filtering in a slow oscilations.
% ESCI - UPF/UPC/UB. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Setting the parameters
%
%   theta = threshold d'activació
%   k = pendent velocitat de canvi
%   Ta = escala temporal d'adaptació
%   w = Força connectivitat
%   lambda = efecte d'adaptació
%   r = rate
%   a = adaptació
dt=0.1;
r=0.9;
a=0.4;
w=0.4;
lambda=1;
I=0;
theta=150;
k=100;
Ta=500;
dr = -r + activationfunction(w*r-lambda*a+I,theta,k);
da=(-a + r)/Ta;  % Slow variable Ta

%Process noise covariance matrix
sigmar= 0.003;
sigman= 0.003;
Q  = [sigmar^2 0; 0 sigman^2];
%Initial Points
m0 = [dr;da;];
P0 = 0.001*eye(2);

QL = chol(Q,'lower');
steps=500;
%%
% Simulating data
%
T = [];
X = [];
Y = [];
t=0;
x= m0;
randn('state',33);
for k=1:steps
    x = [x(1)+(-r + activationfunction((w*r-lambda*a+I),theta,k))*dt;
         x(2)+((-a + r)/Ta)*dt];
    noise = QL*randn(size(QL,1),1);
    y= x;
    x = x + noise;
    t = t + dt;
    X = [X x];
    Y = [Y y];
    T = [T t];
end
figure, plot(T,X(1,:),'.',T,Y(1,:),'-');
title('Rate simulation')
legend('Measurements','Trajectory');
xlabel('{\it x}_1');
ylabel('{\it x}_2');

figure, plot(T,X(2,:),'.',T,Y(2,:),'-');
title('Adaptation simulation')
legend('Measurements','Trajectory');
xlabel('{\it x}_1');
ylabel('{\it x}_2');
rmse_raw = sqrt(mean(sum((Y(2,:) - X(2,:)).^2,1)))

%%
% Filter KF
%
%%
% Kalman filter
%
Vest = m0; %value estimate
Eest = P0; %error in estimate ( error covariance matrix)
kf_m = zeros(size(m,1),size(X,2));
kf_m(:,1) = m0;
Em= Q; % error measurement is constant
for k=2:size(X,2)
    Vm= X(:,k); %value measurement
    KG = Eest/(Eest+Em); % Kalman Gain
    kf_m(:,k) = Vest+KG*[Vm-Vest];
    Eest= [eye(2)-KG]*Eest; %update error of estimate
    Vest = [Vest(1)+(-r + activationfunction((w*r-lambda*a+I),theta,k))*dt;Vest(2)+((-a + r)/Ta)*dt];
end

figure, plot(T,Y(1,:),'k-',T,X(1,:),'r.',T,kf_m(1,:),'b--');
title('KF estimate');
legend('True','Measurements','Estimate');

rmse_kf = sqrt(mean((X(1,:)-kf_m(1,:)).^2))


figure, plot(T,Y(2,:),'k-',T,X(2,:),'r.',T,kf_m(2,:),'b--');
title('KF estimate');
legend('True','Measurements','Estimate');

rmse_ckf = sqrt(mean((X(2,:)-MM(2,:)).^2))