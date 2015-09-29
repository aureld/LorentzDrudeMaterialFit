 
path = '.';
filename = 'Silver_Johnson';
extension = '.txt';
data = load(fullfile(path,[filename extension]));

%lambda = freqtoLambda(data(:,1))*1e9;
lambda = data(:,1);
n = data(:,2);
k = data(:,3);
epsR = n.^2 - k.^2;
epsI = 2*n.*k;
[epsRcalc,epsIcalc] = LD(lambda*1e-6,fullfile(path,[filename '_coefs.txt']));

ncalc = sqrt( (sqrt( epsRcalc.^2 + epsIcalc.^2 ) + epsRcalc) /2 );
kcalc = sqrt( (sqrt( epsRcalc.^2 + epsIcalc.^2 ) - epsRcalc) /2 );

error = norm(epsRcalc - epsR).^2 + norm(epsIcalc - epsI).^2 / (norm(epsR).^2 + norm(epsI).^2) ;



fig = figure(1);
subplot(2,2,1);
plot(lambda, epsR,lambda,epsRcalc)
xlabel('wavelength (nm)');
ylabel('\Re(\epsilon)');
legend('data', 'LD')
subplot(2,2,2) 
plot(lambda, epsI,lambda,epsIcalc)
xlabel('wavelength (nm)');
ylabel('\Im(\epsilon)');
legend('data', 'LD')

subplot(2,2,3)
plot(lambda, n,lambda,ncalc);
xlabel('wavelength (nm)');
ylabel('n');
legend('data', 'LD')
subplot(2,2,4)
plot(lambda, k,lambda,kcalc);
xlabel('wavelength (nm)');
ylabel('k');
legend('data', 'LD')
suplabel(sprintf('avg. error: %3.3f', error));

print(fig, fullfile(path,[filename '_figs.jpg']), '-djpeg');


fig = figure(2);
data = load(fullfile(path,[filename '_fitness.csv']));
plot(data(:,1),data(:,2));
xlabel('generation #');
ylabel('fitness score');
title('solution progression');



