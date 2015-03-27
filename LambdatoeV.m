function [ eV ] = LambdatoeV( lambda )
% converts lambda (um) to energy (eV)
% AD - june 2014
hc = 1.23984193; 

eV =  hc ./ lambda;

end