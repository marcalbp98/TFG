% I. System
%
% x+ = F_x * x + F_u * u + F_n * n
% y = H * x + v
%
% x : state vector   - P : cova. matrix
% u : control vector
% n : perturbation vector   - Q : cova. matrix
% y : measurement vector
% v : measurement noise   - R : cova. matrix
%
% F_x : transition matrix
% F_u : control matrix
% F_n : pert. matrix
% H : measurement matrix
%
% II. Initialization
%
% Define F_x, F_u, and H
%
% Precise x, P, Q, R (initial state)
%
% III. Temporal loop
%
% IIIa. Prediction of mean(x) and P at the arrival of u
%
%       x+ = F_x * x + F_u * u
%       P+ = F_x * P * F_x' + F_n * Q * F_n'
%
% IIIb. correction of mean(x) and P at the arrival of y
%
%       e = H * x       - expectation
%       E = H * P * H'
%
%       z = y - e       - innovation
%       Z = R + E
%
%       K = P * H' * Z^-1   -Kalman gain
%
%       x+ = x + K * z      - Best estimate we can have in our system,
%                             considering uncertanty in perturbation and in sensors
%       P+ = P - K * H * P // P - K * Z * K' #second equation is simetric