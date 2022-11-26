import numpy as np
import matplotlib.pyplot as plt

# KDP晶体色散公式（3.2.2）
# 负单轴晶体第一类相位匹配角公式（3.2.5）

def KDP(wavelength):
    no_2 = 2.259276+0.78056/(77.26408*wavelength**2-1)+(0.032513*wavelength**2)/(0.0025*wavelength**2-1)
    ne_2 = 2.132668+0.703319/(81.42631*wavelength**2-1)+(0.00807*wavelength**2)/(0.0025*wavelength**2-1)
    return no_2,ne_2

def main():
    num = 2                           # 第一问
    if num == 1:
        lambda_p = 347.2 / 1000       # nm->um
        lambda_s = 694.4 / 1000
        lambda_i = 694.4 / 1000
    else:
        lambda_p = np.linspace(250,450,901) / 1000
        lambda_s = 2*lambda_p
        lambda_i = 2*lambda_p

    # 计算折射率
    [no_w_2,ne_w_2] = KDP(lambda_i)
    [no_2w_2,ne_2w_2] = KDP(lambda_p)
    # 计算匹配角
    inside = (ne_2w_2/no_w_2) * (no_2w_2-no_w_2) / (no_2w_2-ne_2w_2)
    theta = np.arcsin(np.sqrt(inside)) *180 / np.pi

    d_theta = np.abs(theta-50.5480)
    linex = np.linspace(300,450,101);       # 画 y=8 的虚线
    liney = 8 * np.ones(np.size(linex))

    # 如果是第一问，请把后面的代码注释掉
    plt.plot(lambda_p*1000,d_theta,'k-.'),plt.plot(linex,liney,'r--')
    plt.xlim([300,450]),plt.ylim([0,8.5]),plt.grid()
    plt.xlabel(r'$\lambda$ / nm'),plt.ylabel(r'$\Delta$$\theta$ / °')
    plt.show()

if __name__ == "__main__":
    main()