function [ x1 ] =activationfunction(x,theta,k)
x1 =1/(1+exp(-(x-theta)/k));
end