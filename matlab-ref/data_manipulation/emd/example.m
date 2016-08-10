%Reproduces in some extent the experiment shown in page 7 of paper:
% Yannis Kopsinis and Steve McLaughlin, “Improved EMD Using
% Doubly-Iterative Sifting and High Order Spline Interpolation,” EURASIP
% Journal on Advances in Signal Processing, vol. 2008, Article ID 128293.
%
% The time-frequency toolbox (freeware) is needed for the signal
% construction . You can download from: http://tftb.nongnu.org/

N = 30000;% # of data samples
t = 1:N;
p = N/4;% period of the 2 sinusoidal FM's

% sinusoidal FM 1
f1=1/25;
a1=1;
[xd(1,:),IF(1,:)]=fmconst(N,f1/(N/5000));
xd(1,:)=real(xd(1,:));
xd(1,:) = a1.*xd(1,:);
AA(1,:)= a1*ones(1,N);

% sinusoidal FM 2
f2=1/55;
a2_min=0.5;
a2_max=15;
a2=a2_min:(a2_max-a2_min)/(N-1):a2_max;
[xd(2,:),IF(2,:)]=fmconst(N,f2/(N/5000));
xd(2,:)=real(xd(2,:));
xd(2,:) = a2.*xd(2,:);
AA(2,:)= a2;
AM = AA;
x = xd(1,:)+xd(2,:);
figure;
subplot(3,1,1);
plot(linspace(0.5,15,N),xd(2,:),'r');hold on;plot(linspace(0.5,15,N),xd(1,:));
set(gca,'ytick',[])


options = emdoptimset('stop','Nsifts','N',500,'imfs',1);
IMF = emd(x,options);
subplot(3,1,2)
plot(IMF)

options = emdoptimset(@di,'stop','Nsifts','N',500,'imfs',1,'in_N',20,'in_Nskip',50);
IMF = emd(x,options);
subplot(3,1,3)
plot(IMF)

