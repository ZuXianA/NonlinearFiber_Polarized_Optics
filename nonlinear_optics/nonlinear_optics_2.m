clear all;close all;clc;

% KDP晶体色散公式（3.2.2）
% 负单轴晶体第一类相位匹配角公式（3.2.5）

num = 2;                                       % 第二问
if num == 1
    lambda_p = 347.2 / 1000;        % nm->um
    lambda_s = 694.4 / 1000;
    lambda_i = 694.4 / 1000;
else
    lambda_p = linspace(250,450,901) / 1000;
    lambda_s = 2*lambda_p;
    lambda_i = 2*lambda_p;
end

% 计算折射率
[no_w_2,ne_w_2] = KDP(lambda_i);
[no_2w_2,ne_2w_2] = KDP(lambda_p);

% 计算匹配角
inside = (ne_2w_2./no_w_2) .* (no_2w_2-no_w_2) ./ (no_2w_2-ne_2w_2);
theta = asin(sqrt(inside)) *180 / pi;

% 如果是第一问，请把后面的代码注释掉
d_theta = abs(theta-50.5480);
linex = linspace(300,450,91);       % 画 y=0.8 的虚线
liney = 8 * ones(size(linex));
plot(lambda_p*1e3,d_theta,'k.-',linex,liney,'r--'),grid on;
xlabel('\lambda / nm'),ylabel('\Delta\theta / °'),ylim([0,8.5]);

function [no_2, ne_2] = KDP(wavelength) 
no_2 = 2.259276+0.78056./(77.26408*wavelength.^2-1)+(0.032513*wavelength.^2)./(0.0025*wavelength.^2-1);
ne_2 = 2.132668+0.703319./(81.42631*wavelength.^2-1)+(0.00807*wavelength.^2)./(0.0025*wavelength.^2-1);
end
