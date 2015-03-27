function [ lambda ] = eVtoLambda( eV )
% converts electron-volts to lambda (um)
% AD - june 2014
hc = 1.23984193; 

lambda =  hc / eV;

end