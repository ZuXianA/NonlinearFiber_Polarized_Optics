clear all;close all;clc;

x = linspace(0.7,2,1000);   % 波长范围 0.7~2 um

% figure(1);
% subplot(211),plot(x,Sellmier(x)),title('Sellmier'),grid on;
% xlabel('\lambda  /  \mum'),ylabel('index of refraction');
% subplot(212),plot(x,Newton(x)),title('Newton'),grid on;
% xlabel('\lambda  /  \mum'),ylabel('index of refraction');

figure(2);
subplot(221),plot(x,Sellmier(x)),title('Sellmier'),grid on;
xlabel('\lambda  /  \mum'),ylabel('index of refraction');
subplot(223),plot(x,Newton(x)),title('Newton'),grid on;
xlabel('\lambda  /  \mum'),ylabel('index of refraction');

[D,p1] = dispersionD(x);
subplot(222),plot(x(3:end),D,'b',x(3:end),zeros(length(x(3:end))),'r--');
title('$D=-\frac{\lambda}{c}\frac{d^2n}{d\lambda^2}$','Interpreter','latex');
xlabel('\lambda  /  \mum'),ylabel('色散参量 D');grid on;
text(0.8,-5e-17,'正常色散区',FontSize=13);
text(1.6,3e-17,'反常色散区',FontSize=13);

[deltaTao,v,p] = DeltaTao(x,x(p1(1)+2));
subplot(224),plot(x(3:end),deltaTao,[x(p+3),x(p+3)],[0,18e-17],'r--');
title("\Delta\tau = D*\Delta\lambda");grid on;
xlabel('\lambda  /  \mum'),ylabel('时延差 \Delta\tau');

fprintf('D0 = %g\n', x(p1(1)+2));

function n = Sellmier(wavelength)       % 构造 Sellmeier 公式
B1 = 0.6961663;
B2 = 0.4079426;
B3 = 0.8974794;
lambda1 = 0.0684043;
lambda2 = 0.1162414;
lambda3 = 9.896161;
n = B1*wavelength.^2 ./ (wavelength.^2-lambda1^2) + B2*wavelength.^2 ./ (wavelength.^2-lambda2^2) + B3*wavelength.^2 ./ (wavelength.^2-lambda3^2);
n = sqrt(n+1);
end

function n = Newton(wavelength)         % 构造 Newton 公式
c1 = 1.45084;
c2 = -0.00343;
c3 = 0.00292;
n = c1 + c2*wavelength.^2 + c3*wavelength.^(-2);
end

function [D,p] = dispersionD(wavelength)     % 色散参量 D 的波长表达式，图标题写清楚了具体形式
c = 299792458;   % 光速                      % p 是 D=0 所在位置的下标
nn = Sellmier(wavelength);
d2n = diff(nn,2);
D = -wavelength(3:end) .* d2n / c;
p = find(D>=0);
end

function [deltaTao,v,p] = DeltaTao(wavelength,D0)    % 归一化时延差用色散参量表示，图标题写清楚了具体形式
DeltaLamda = wavelength(3:end) - D0;
deltaTao = dispersionD(wavelength) .* DeltaLamda;
[v,p] = min(deltaTao>0);
end
