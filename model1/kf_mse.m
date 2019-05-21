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
Lambda=0.5;
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

%Measurements noise
R= 0.00001;
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
for k=1:steps
    x = [x(1)+(-r + activationfunction((w*r-Lambda*a+I),Theta,k))*dt;
         x(2)+((-a + r)/Ta)*dt];
    x = x + QL*randn(size(QL,1),1);
    y = x + sqrt(R)*randn;
    t = t + dt;
    X = [X x];
    Y = [Y y];
    T = [T t];
end

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
    m = m + K*(Y(:,k) - H*m);
    P = P - K*S*K';

    kf_m(:,k) = m;
    kf_P(:,:,k) = P;
end

rmse_kf_r = sqrt(mean((kf_m(1,:)-X(1,:)).^2));
rmse_kf_a = sqrt(mean((kf_m(2,:)-X(2,:)).^2));

%clearvars -except rmse_raw_r rmse_raw_a