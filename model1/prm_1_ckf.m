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
R = 0.003;
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
% Filter CKF
%

m = m0;
P = P0;

n = size(m,1);

XI = sqrt(n) * [eye(n) -eye(n)];
W  = ones(1,2*n)/(2*n);

%
% Do the filtering
%  
MM = zeros(size(m,1),length(X));
PP = zeros(size(P,1),size(P,2),length(X));
for k=1:length(X)

    % Form the sigma points for dynamic model
    SX = repmat(m,1,size(XI,2)) + chol(P,'lower') * XI;

    % Propagate through the dynamic model
    HX = [SX(1,:)+(-r + activationfunction((w*r-lambda*a+I),theta,k))*dt; SX(2,:)+((-a + r)/Ta)*dt];

    % Compute the predicted mean and covariance
    m = zeros(size(m));
    P = zeros(size(P));
    for i=1:size(HX,2)
        m = m + W(i) * HX(:,i);
    end
    for i=1:size(HX,2)
        P = P + W(i) * (HX(:,i) - m) * (HX(:,i) - m)';
    end
    P = P + Q;

    % Form sigma points for measurement step and
    % propagate throught the measurement model
    SX = repmat(m,1,size(XI,2)) + chol(P,'lower') * XI;
    HY = sin(SX(1,:)); 

    % Compute the updated mean and covariance
    mu = zeros(size(HY,1),1);
    S  = zeros(size(HY,1),size(HY,1));
    C  = zeros(size(SX,1),size(HY,1));
    for i=1:size(SX,2)
        mu = mu + W(i) * HY(:,i);
    end
    for i=1:size(SX,2)
        S = S + W(i) * (HY(:,i) - mu) * (HY(:,i) - mu)';
        C = C + W(i) * (SX(:,i) - m) * (HY(:,i) - mu)';
    end
    S = S + R;

    % Compute the gain and updated mean and covariance  
    K = C/S;
    m = m + K*(X(k) - mu);
    P = P - K*S*K';

    MM(:,k) = m;
    PP(:,:,k) = P;
end  

figure, plot(T,Y(1,:),'k-',T,X(1,:),'r.',T,MM(1,:),'b--');
title('CKF estimate');
legend('True','Measurements','Estimate');

rmse_ckf = sqrt(mean((X(1,:)-MM(1,:)).^2))


figure, plot(T,Y(2,:),'k-',T,X(2,:),'r.',T,MM(2,:),'b--');
title('CKF estimate');
legend('True','Measurements','Estimate');

rmse_ckf = sqrt(mean((X(2,:)-MM(2,:)).^2))