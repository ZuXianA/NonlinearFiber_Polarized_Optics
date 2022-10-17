clear all;close all;clc;

c = 299792458;
D = 0.092 / 4 * (1550-1312.^4./1550^3) * 1e-12 / 1e-9 / 1e3;
beta2 = -1550e-9^2 / 2 / pi / c*D;

num_of_L = 11;
Ls = linspace(0,100,num_of_L);          % 光纤长度

N = 2^8;
twin = 50e-12;
dt = twin/N;
df = 1/twin;
fwin = 1/dt;
t = linspace(-twin/2,twin/2-dt,N);
f = linspace(-fwin/2,fwin/2-df,N);

FWHM = 1e-12;                       % full width at half maxima
morder = 1;                         % mth order super-Gaussian pulse; m=1 for Gaussian pulse
Nc = 1;                             % 决定是否有啁啾

a = zeros(length(Ls),N);               % 时域初始化
A = zeros(length(Ls),N);               % 频域初始化
a0 = zeros(length(Ls),N);               % 时域初始化

for row = 1:length(Ls)          % 计算每一个Nc对应的数值
    a(row,:) = exp(-log(2) / 2*abs(2.0*t/FWHM).^(2*morder)).*exp(-1i*(2*log(2)/FWHM^2*sqrt(Nc^2-1))*t.^2);
    A(row,:) = fftshift(ifft(fftshift(a(row,:))));    % 傅里叶变换(选取的基不同，所以是ifft)
    A(row,:) =  A(row,:).*exp(1i*0.5*beta2*Ls(row)*(2*pi*f).^2);
    a0(row,:) = fftshift(fft(fftshift(A(row,:))));
end

% 绘制三维线图 
u = a0;
[X,Y] = meshgrid(1:size(u,1),1:size(u,2));
Z = u(1:num_of_L,:);
for row = 1:length(Ls)  
    Z(row,:) = Z(row,:).*conj(Z(row,:));   % 二维坐标系下的纵坐标(模值)
end
for j = 1:length(Ls)
    Y(:,j) = t*1e12;      % 二维坐标系下的横坐标(波长,um)
end
for j = 1:num_of_L
    X(:,j) = Ls(j);      % 三维坐标系下的横坐标(光纤长度,m)
end

plot3(X,Y,Z),grid on
if Nc == 1
    title('Gaussian pulse (Chirp-free, Nc=1)',FontSize=18);
else
    title('Chirped Gaussian pulse (Nc=2)',FontSize=18);
end
xlabel('Fiber length / m'),ylabel('time'),zlabel('magnitude');
