function [ omega ] = eV2freq( eV )
% converts electron-volts to frequency
% AD - june 2014
hbar = 6.58211928e-16; %Planck's constant

omega =  eV / hbar;

end

