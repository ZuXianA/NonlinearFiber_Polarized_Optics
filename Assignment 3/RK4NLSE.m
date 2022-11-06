clear all; close all; clc;

global f;                   % 定义全局变量，后面有函数调用
global beta2;
global gamma;

c = 299792458;                       % 光速 m/s
D = 17 * 1e-12/1e-9/1e3;        % 17 ps/nm/km
beta2=-1550e-9^2/2/pi/c*D;
gamma=0.09;                            % 非线性系数，不要乱改
 % 离散化网络
N=2^10;  twin=100e-12;  dt=twin/N;  df=1/twin;  fwin=1/dt;
t=linspace(-twin/2,twin/2-dt,N)';  f=linspace(-fwin/2,fwin/2-df,N)';
 
P0 = 1;                 % 初始脉冲峰值功率，不要乱改
T0=1e-12;            % 脉冲时域宽度
LD=T0^2/abs(beta2);                         fprintf('Dispersion length is %gm\n',LD);
LNL = 1 / gamma / P0;                       fprintf('Nonlinear length is %gm\n',LNL);                 
LP = pi * T0^2 / 2 / abs(beta2);           % 计算孤子周期--http://www.aikelabs.com/wiki/443.htm  
% fprintf('Soliton evolution period is %gm\n',LP);

a0 =  sech(t/T0);            % 输入脉冲，不要乱改，需要配合 γ 和 P0
A0 = fftshift(ifft(fftshift(a0)));


L=35;                           % 18 35 43 74
M=250;  dL=L/M;           % M不能取太小，否则计算不准确
 
B=A0;
for k=1:M
    z=(k-1)*dL;               % 四阶 RK 法则
    K1=askdBdz(z,B);  K2=askdBdz(z+0.5*dL,B+0.5*K1*dL);
    K3=askdBdz(z+0.5*dL,B+0.5*K2*dL);  K4=askdBdz(z+dL,B+K3*dL);
    B=B+1/6*(K1+2*K2+2*K3+K4)*dL;
end
A=B.*exp(1i*0.5*beta2*L*(2*pi*f).^2);  
a=fftshift(fft(fftshift(A)));

[left,right,P1] = findFWHM(a);        % 计算经历一段光纤传播后FWHM
t0 = (right-left) * dt / 1.76;               % 脉冲宽度
LD=t0^2/abs(beta2);                         fprintf('Dispersion length is %gm\n',LD);
LNL = 1 / gamma / P1;                       fprintf('Nonlinear length is %gm\n',LNL);   

subplot(311),hold on;legend;axis([-10 10 0 inf]);
plot(t*1e12,a0.*conj(a0),'r','DisplayName','Input Pulse Shape');
plot(t*1e12,a.*conj(a),'k-.','DisplayName','Output Pulse Shape')

subplot(312),hold on;legend;axis([1.54 1.56 0 inf]);
plot(c./(f+c/1550e-9)*1e6,A0.*conj(A0),'r','DisplayName','Input Pulse Spectrum');
plot(c./(f+c/1550e-9)*1e6,A.*conj(A),'k-.','DisplayName','Output Pulse Spectrum');

subplot(313),plot(t(2:end)*1e12,-diff(phase(a)) / 2 / pi);      % 画瞬时频率      

function  [low,high,v] = findFWHM(sig)
    sig = sig.*conj(sig);
    [v,p] = max(sig);
    pos = find(sig>=0.5*v);
    low = pos(1); high = pos(end); 
end
function res=askdBdz(z,B)
    global f;
    global beta2;
    global gamma;
    a=fftshift(fft(fftshift(B.*exp(1i*0.5*beta2*z*(2*pi*f).^2))));
    res=1i*gamma*fftshift(ifft(fftshift(a.*conj(a).*a))).*exp(-1i*0.5*beta2*z*(2*pi*f).^2);
end


