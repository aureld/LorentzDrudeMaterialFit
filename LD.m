function varargout = LD(lambda,filename)

%***********************************************************************
% Physical constants
%***********************************************************************
twopic = 1.883651567308853e+09; % twopic=2*pi*c where c is speed of light
omegalight = twopic*(lambda.^(-1)); % angular frequency of light (rad/s)
invsqrt2 = 0.707106781186547;  % 1/sqrt(2)
ehbar = 1.51926751447914e+015; % e/hbar where hbar=h/(2*pi) and e=1.6e-19

data2 = load(filename);
if (mod(length(data2),3) ==0) % if epsInf is absent, we set it to 1.0
    data = [1.0 ;data2];
else
    data = data2;
end
eps_inf = data(1);

%drude
f(1) = data(2);
omegap = data(3);
Gamma(1) = data(4);
omega(1) = 0.0;

%Lorentz
for i=1:(length(data)/3-1)
    f(i+1) = data(i*3+2);
    omega(i+1) = data(i*3+3);
    Gamma(i+1) = data(i*3+4);
end
order = length(omega);

 epsilon_D = ones(size(lambda)) - ((f(1)*omegap^2) *...
    (omegalight.^2 + 1i*Gamma(1)*omegalight).^(-1));

epsilon_L = zeros(size(lambda));
    % Lorentzian contributions
    for k = 2:order
        epsilon_L = epsilon_L + (f(k)*omegap^2)*...
            (((omega(k)^2)*ones(size(lambda)) - omegalight.^2) -...
            1i*Gamma(k)*omegalight).^(-1);
    end

    % Drude and Lorentz contributions combined
    epsilon = epsilon_D + epsilon_L;

switch nargout
    case 1 % one output variable assigned
        varargout{1} = epsilon;

    case 2 % two output variables assigned

        % Real part of dielectric constant
        varargout{1} = real(epsilon);

        % Imaginary part of dielectric constant
        varargout{2} = imag(epsilon);
end
