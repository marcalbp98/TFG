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
dt=0.01;
r=0.5;
a=0.4;
w=0.4;
Lambda=0.5;
I=0;

Theta=0.2;
k=30;
Ta=400;
%dr = (-r + activationfunction(w*r-Lambda*a+I,Theta,k))*dt;
%da=((-a + r)/Ta)*dt;  % Slow variable Ta

%Process noise covariance matrix
sigmar= 0.003;
sigman= 0.003;
Q  = [sigmar^2 0; 0 sigman^2];

%Measurements noise
R= 0.00001;
%Initial Points
m0 = [r;a;];
P0 = 0.001*eye(2);

QL = chol(Q,'lower');
steps=5000;
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
    x = [x(1)+(-x(1) + activationfunction((w*x(1)-Lambda*x(2)+I),Theta,k))*dt;
         x(2)+((-x(2) + x(1))/Ta)*dt];
    x = x + QL*randn(size(QL,1),1);
    y = x + sqrt(R)*randn;
    t = t + dt;
    X = [X x];
    Y = [Y y];
    T = [T t];
end
% figure, plot(T,Y(1,:),'.',T,X(1,:),'-');
% title('Rate simulation')
% legend('Measurements','Trajectory');
% xlabel('{\it x}_1');
% ylabel('{\it x}_2');
% 
% figure, plot(T,X(2,:),'k-',T,Y(2,:),'.');
% title('Adaptation simulation')
% legend('Measurements','Trajectory');
% xlabel('{\it x}_1');
% ylabel('{\it x}_2');

rmse_raw_r = sqrt(mean(sum((Y(1,:) - X(1,:)).^2,1)))
rmse_raw_a = sqrt(mean(sum((Y(2,:) - X(2,:)).^2,1)))

%%
% Kalman filter
%
H=eye(2);
m = m0;
P = P0;
R=Q;
kf_m = zeros(size(m,1),size(Y,2));
kf_P = zeros(size(P,1),size(P,2),size(Y,2));
for k=1:size(Y,2)
    epsilon= exp((-I - m(1)*w +Theta+ m(2)*Lambda)/k);
    A= [-1+((epsilon)*w)/(((1 + epsilon)^2)*k) -(epsilon*Lambda)/(((1 + epsilon)^2)*k);1/Ta (-1)/Ta];%jacobia
    m = [m(1)+(-m(1) + activationfunction((w*m(1)-Lambda*m(2)+I),Theta,k))*dt;
         m(2)+((-m(2) + m(1))/Ta)*dt];
    P = A*P*A' + Q;

    S = H*P*H' + R;
    K = P*H'/S;
    m = m + K*(Y(:,k) - H*m);
    P = P - K*S*K';

    kf_m(:,k) = m;
    kf_P(:,:,k) = P;
end

figure, plot(T,X(1,:),'k-',T,Y(1,:),'r.',T,kf_m(1,:),'b--');
title('KF estimate for firing rate');
legend('True','Measurements','Estimate');
xlabel('Time{\it (s)}');
ylabel('Voltage{\it (V)}');

rmse_kf_r = sqrt(mean((X(1,:)-kf_m(1,:)).^2))


figure, plot(T,X(2,:),'k-',T,Y(2,:),'r.',T,kf_m(2,:),'b--');
title('KF estimate for adaptation rate');
legend('True','Measurements','Estimate');
xlabel('Time{\it (s)}');
ylabel('Voltage{\it (V)}');

rmse_kf_a = sqrt(mean((X(2,:)-kf_m(2,:)).^2))

%clearvars -except rmse_raw_r rmse_raw_a