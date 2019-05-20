%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimate pendulum state with CKF and CRTS as in Examples 6.2
% and 10.2 of the book
%
% Simo Sarkka (2013), Bayesian Filtering and Smoothing,
% Cambridge University Press. 
%
% Last updated: $Date: 2013/08/26 12:58:41 $.
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Simulate data
%
    pendulum_sim;

%%
% Filter
%

    m = m0;
    P = P0;
    
    n = size(m,1);

    XI = sqrt(n) * [eye(n) -eye(n)];
    W  = ones(1,2*n)/(2*n);
    
    %
    % Do the filtering
    %  
    MM = zeros(size(m,1),length(Y));
    PP = zeros(size(P,1),size(P,2),length(Y));
    for k=1:length(Y)

        % Form the sigma points for dynamic model
        SX = repmat(m,1,size(XI,2)) + chol(P,'lower') * XI;

        % Propagate through the dynamic model
        HX = [SX(1,:)+SX(2,:)*DT; SX(2,:)-g*sin(SX(1,:))*DT];
  
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
        m = m + K*(Y(k) - mu);
        P = P - K*S*K';
    
        MM(:,k) = m;
        PP(:,:,k) = P;
    end  
  
    h = plot(T,Y,'k.',T,X(1,:),'r-',T,MM(1,:),'b--');
    set(h,'Linewidth',5);
    title('CKF estimate');
    legend('Measurements','True','Estimate');
  
    rmse_ckf = sqrt(mean((X(1,:)-MM(1,:)).^2))

%%
% Smoother
%

    ms = m;
    Ps = P;
    MS = zeros(size(m,1),length(Y));
    PS = zeros(size(P,1),size(P,2),length(Y));
    MMS(:,k) = m;
    PPS(:,:,k) = P;
    for k=size(MM,2)-1:-1:1
        m = MM(:,k);
        P = PP(:,:,k);
    
        % Form the sigma points for dynamic model
        SX = repmat(m,1,size(XI,2)) + chol(P,'lower') * XI;

        % Propagate through the dynamic model
        HX = [SX(1,:)+SX(2,:)*DT; SX(2,:)-g*sin(SX(1,:))*DT];
 
        % Compute the predicted mean and covariance
        % and the cross-covariance D.
        mp = zeros(size(m));
        Pp = zeros(size(P));
        D  = zeros(size(P));
        for i=1:size(HX,2)
            mp = mp + W(i) * HX(:,i);
        end
        for i=1:size(HX,2)
            Pp = Pp + W(i) * (HX(:,i) - mp) * (HX(:,i) - mp)';
            D = D + W(i) * (SX(:,i) - m) * (HX(:,i) - mp)';
        end
        Pp = Pp + Q;

        Ck = D/Pp;
        ms = m + Ck*(ms - mp);
        Ps = P + Ck*(Ps - Pp)*Ck';
        MMS(:,k) = ms;
        PPS(:,:,k) = Ps;
    end
  
    h = plot(T,Y,'k.',T,X(1,:),'r-',...
             T,MM(1,:),'b--',T,MMS(1,:),'g--');
    set(h,'Linewidth',5);
    title('CKF and CRTS estimates');
    legend('Measurements','True','CKF','CRTS');
    
    rmse_crts = sqrt(mean((X(1,:)-MMS(1,:)).^2))
    
%%
% Plot the filtering result
%

    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MM(1,:),'r');
    legend('True Angle','Measurements','CKF Estimate');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
    
%%
% Plot the smoothing result
%
    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MMS(1,:),'r');
    legend('True Angle','Measurements','CRTSS Estimate');
    xlabel('Time{\it t}');
ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
