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
%   Lambda = efecte d'adaptació
%   r = rate
%   a = adaptació
dt=0.1;
r=0.9;
a=0.4;
w=0.4;
Lambda=1;
I=0;
Theta=150;
k=100;
Ta=500;
dr = -r + activationfunction(w*r-Lambda*a+I,Theta,k);
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
    x = [x(1)+(-r + activationfunction((w*r-Lambda*a+I),Theta,k))*dt;
         x(2)+((-a + r)/Ta)*dt];
    noise = QL*randn(size(QL,1),1);
    y= x;
    x = x + noise;
    t = t + dt;
    X = [X x];
    Y = [Y y];
    T = [T t];
end
% figure, plot(T,X(1,:),'.',T,Y(1,:),'-');
% title('Rate simulation')
% legend('Measurements','Trajectory');
% xlabel('{\it x}_1');
% ylabel('{\it x}_2');
% 
% figure, plot(T,X(2,:),'.',T,Y(2,:),'-');
% title('Adaptation simulation')
% legend('Measurements','Trajectory');
% xlabel('{\it x}_1');
% ylabel('{\it x}_2');

rmse_raw_r = sqrt(mean(sum((Y(2,:) - X(1,:)).^2,1)))
rmse_raw_a = sqrt(mean(sum((Y(2,:) - X(2,:)).^2,1)))

%%
% Kalman filter
%
H=eye(2);
m = m0;
P = P0;
R=Q*0.001;
A= [-1+(exp((-r*w + Theta + a*Lambda)/k)*w)/(((1 + exp((-r*w +Theta + a*Lambda)/k))^2)*k) -((exp((-r*w +Theta + a*Lambda)/k)*Lambda)/(((1 + exp((-r*w +Theta+ a*Lambda)/k))^2)*k));1/Ta (-1)];%jacobia
kf_m = zeros(size(m,1),size(Y,2));
kf_P = zeros(size(P,1),size(P,2),size(Y,2));
for k=1:size(Y,2)
    m = [m(1)+(-r + activationfunction((w*r-Lambda*a+I),Theta,k))*dt;
         m(2)+((-a + r)/Ta)*dt];
    P = A*P*A' + Q;

    S = H*P*H' + R;
    K = P*H'/S;
    m = m + K*(X(:,k) - H*m);
    P = P - K*S*K';

    kf_m(:,k) = m;
    kf_P(:,:,k) = P;
end

figure, plot(T,Y(1,:),'k-',T,X(1,:),'r.',T,kf_m(1,:),'b--');
title('KF estimate for firing rate');
legend('True','Measurements','Estimate');

rmse_kf_r = sqrt(mean((X(1,:)-kf_m(1,:)).^2))


figure, plot(T,Y(2,:),'k-',T,X(2,:),'r.',T,kf_m(2,:),'b--');
title('KF estimate for adaptation rate');
legend('True','Measurements','Estimate');

rmse_kf_a = sqrt(mean((X(2,:)-kf_m(2,:)).^2))