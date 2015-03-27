function [ lambda ] = freqtoLambda( f )
% converts Hz to m
c = 299792458;

lambda =  c ./ f;

end