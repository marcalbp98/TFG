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
Ta=5000;
dr = -r + activationfunction(w*r-lambda*a+I,theta,k);
da=(-a + r)/Ta;  % Slow variable

%drtest=[];
%datest=[];
%for I = 1:100
%    drtest(I,:) = -r + activationfunction((w*r-lambda*a+I),theta,k);
%    datest(I,:)=-a + r;
%end
%plot(drtest,'-')
sigmar= 0.003;
sigman= 0.003;
Q  = [sigmar^2 0; 0 sigman^2];
 %Process noise covariance matrix
R = 0.1;
m0 = [dr;da;];

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
% Kalman filtering
%
m = m0;
P = eye(2);
kf_m = zeros(size(m,1),size(X,2));
kf_P = zeros(size(P,1),size(P,2),size(X,2));
for k=1:size(X,2)
    m = ((-a + r)/Ta)*m;
    P = ((-a + r)/Ta)*P*((-a + r)/Ta)' + Q;

    S = ((-a + r)/Ta)*P*((-a + r)/Ta)' + R;
    K = P*((-a + r)/Ta)'/S;

    m = m + K*(X(2,k) - ((-a + r)/Ta)*m);
    P = P - K*S*K'; %best estimate

    kf_m(:,k) = m;
    kf_P(:,:,k) = P;
end

rmse_raw = sqrt(mean(sum((Y(2,:) - X(2,:)).^2,1)))
rmse_kf = sqrt(mean(sum((kf_m(2,:) - X(2,:)).^2,1)))

figure, plot(T,Y(2,:),'-',T,X(2,:),'.',...
    T,kf_m(2,:),'r-');
title('Kalman filter')
legend('True Trajectory','Measurements','Filter Estimate');
xlabel('{\it x}_1');
ylabel('{\it x}_2');


