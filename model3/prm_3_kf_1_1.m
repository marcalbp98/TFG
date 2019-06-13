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
%excitatory population:
r_e=0.9;
a_e=0.05;
s_e=0.2;
Lambda_e=0.9;
I_e=0;

%inhibitory population
r_i=0.2;
a_i=0.05;
s_i=0.2;
Lambda_i=0.9;
I_i=0;

wii=0.2;
wei=0.4;
wie=0.4;
wee=0.2;
Theta=0.3;
k=0.50;
Delta=0.3;

Te=3;
Ti=5;
Ta=400;
Ts=1000;

% dr_e = -r_e + activationfunction(wee*r_e*s_e-wie*r_i*s_i-Lambda_e*a_e+Ie)
% dr_i = -r_i + activationfunction(wei*r_e*s_e-wii*r_i*s_i-Lambda_i*a_i+Ii)

%Process noise covariance matrix
sigma_r_e= 0.003;
sigma_r_i= 0.003;
sigma_a_e= 0.003;
sigma_a_i= 0.003;
sigma_s_e= 0.003;
sigma_s_i= 0.003;

Q  = [sigma_r_e^2 0 0 0 0 0; 0 sigma_r_i^2 0 0 0 0;0 0 sigma_a_e^2 0 0 0; 0 0 0 sigma_a_i^2 0 0;0 0 0 0 sigma_s_e^2 0; 0 0 0 0 0 sigma_s_i^2];

%Measurements noise
R= Q;
%Initial Points
m0 = [r_e;r_i;a_e;a_i;s_e;s_i];
P0 = 0.001*eye(6);

QL = chol(Q,'lower');
RL=  chol(R,'lower');

steps=50000;
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
    x = [x(1)+((-x(1) + activationfunction(wee*x(1)*x(5)-wie*x(2)*x(6)-Lambda_e*x(3)+I_e,Theta,k))/Te)*dt;
         x(2)+((-x(2) + activationfunction(wei*x(1)*x(5)-wii*x(2)*x(6)-Lambda_i*x(4)+I_i,Theta,k))/Ti)*dt;
         x(3)+((-x(3) + x(1))/Ta)*dt;
         x(4)+((-x(4) + x(2))/Ta)*dt;
         x(5)+((-x(5)-1-Delta*x(1)*x(5))/Ts)*dt;
         x(6)+((-x(6)-1-Delta*x(2)*x(6))/Ts)*dt;
];
    x = x + QL*randn(size(QL,1),1);
    y = x + RL*randn(size(RL,1),1);
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
%%
% Kalman filter
%
H=eye(6);
m = m0;
P = P0;
kf_m = zeros(size(m,1),size(Y,2));
kf_P = zeros(size(P,1),size(P,2),size(Y,2));
for k=1:size(Y,2)
    epsilon_e= exp((-wee*m(1)*m(5)+wie*m(2)*m(6)+Lambda_e*m(3)-I_e+Theta)/k);
    epsilon_i= exp((-wei*m(1)*m(5)+wii*m(2)*m(6)+Lambda_i*m(4)-I_i+Theta)/k);
    A= [(-1+(epsilon_e*m(5)*wee)/(((1+epsilon_e)^2)*k))/Te (-(epsilon_e*wie*m(6))/(((1 + epsilon_e)^2)*k))/Te ((epsilon_e*Lambda_e)/(((1 + epsilon_e)^2)*k))/Te 0 ((epsilon_e*wee*m(1))/(((1 + epsilon_e)^2)*k))/Te (-(epsilon_e*wie*m(2))/(((1 + epsilon_e)^2)*k))/Te;
        ((epsilon_e*m(5)*wei)/(((1+epsilon_e)^2)*k))/Ti -1+(-(epsilon_e*wii*m(6))/(((1 + epsilon_e)^2)*k))/Ti 0 ((epsilon_e*Lambda_i)/(((1 + epsilon_e)^2)*k))/Ti ((epsilon_e*wei*m(1))/(((1 + epsilon_e)^2)*k))/Ti (-(epsilon_e*wii*m(2))/(((1 + epsilon_e)^2)*k))/Ti;
        1/Ta 0 -1 0 0 0;
        0 1/Ta 0 -1 0 0;
        -Delta*m(5)/Ts 0 0 0 -1-Delta*m(1)/Ts 0;
        0 -Delta*m(6)/Ts 0 0 0 -1-Delta*m(2)/Ts;];%jacobia
    m = [m(1)+((-m(1) + activationfunction(wee*m(1)*m(5)-wie*m(2)*m(6)-Lambda_e*m(3)+I_e,Theta,k))/Te)*dt;
         m(2)+((-m(2) + activationfunction(wei*m(1)*m(5)-wii*m(2)*m(6)-Lambda_i*m(4)+I_i,Theta,k))/Ti)*dt;
         m(3)+((-m(3) + m(1))/Ta)*dt;
         m(4)+((-m(4) + m(2))/Ta)*dt;
         m(5)+((-m(5)-1-Delta*m(1)*m(5))/Ts)*dt;
         m(6)+((-m(6)-1-Delta*m(2)*m(6))/Ts)*dt;];
    P = A*P*A' + Q;

    S = H*P*H' + R;
    K = P*H'/S;
    m = m + K*(Y(:,k) - H*m);
    P = P - K*S*K';

    kf_m(:,k) = m;
    kf_P(:,:,k) = P;
end

figure, plot(T,X(1,:),'k-',T,Y(1,:),'r.',T,kf_m(1,:),'b--');
title('KF estimate for firing rate of exitatory neurons');
legend('True','Measurements','Estimate');
xlabel('Time{\it (s)}');
ylabel('Voltage{\it (V)}');

rmse_kf_r = sqrt(mean((X(1,:)-kf_m(1,:)).^2))


figure, plot(T,X(2,:),'k-',T,Y(2,:),'r.',T,kf_m(2,:),'b--');
title('KF estimate for firing rate of inhibitory neurons');
legend('True','Measurements','Estimate');
xlabel('Time{\it (s)}');
ylabel('Voltage{\it (V)}');

rmse_kf_a = sqrt(mean((X(2,:)-kf_m(2,:)).^2))

clearvars -except rmse_raw_r rmse_raw_a