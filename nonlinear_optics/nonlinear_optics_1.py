import numpy as np
import scipy.misc as sm
import matplotlib.pyplot as plt

c = 299792458                   # 光速 m/s

def Sellmeier(wavelength):      # 构造 Sellmeier 公式
    B1 = 0.6961663
    B2 = 0.4079426
    B3 = 0.8974794
    lambda1 = 0.0684043
    lambda2 = 0.1162414
    lambda3 = 9.896161
    n = B1*wavelength**2 / (wavelength**2-lambda1**2) + B2*wavelength**2 / (wavelength**2-lambda2**2) + B3*wavelength**2 / (wavelength**2-lambda3**2)
    return np.sqrt(n+1)

def Newton(wavelength):         # 构造 Newton 公式
    c1 = 1.45084
    c2 = -0.00343
    c3 = 0.00292
    n = c1 + c2*wavelength**2 + c3*wavelength**(-2)
    return n

def dispersionD(wavelength):    # 色散参量 D 的波长表达式，图标题写清楚了具体形式
    d2n = sm.derivative(func=Sellmeier,x0=wavelength,dx=1e-6,n=2)
    D = -wavelength * d2n / c
    p = np.where(D>=0)[0][0]         # p返回的是零点所在位置
    return D,p

def DeltaTao(wavelength,D,D0):       # 归一化时延差用色散参量表示，图标题写清楚了具体形式
    DeltaLambda = wavelength - D0
    deltaTao = D * DeltaLambda
    return deltaTao

def main():
    x = np.linspace(start=0.7,stop=2,num=1000)    # 波长的范围 1~2 um
    nnS = Sellmeier(wavelength=x)
    nnN = Newton(wavelength=x)

    line0 = np.zeros(len(x))                    # 画一条 D=0 的直线
    D,p1 = dispersionD(wavelength=x)
    
    '''下面是画图区域'''
    plt.subplot(221),plt.plot(x,nnS),plt.title('Sellmeier'),plt.grid()
    plt.xlabel(r'$\lambda$ / $\mu$m',loc='right'),plt.ylabel('index of refraction')

    plt.subplot(223),plt.plot(x,nnN),plt.title('Newton'),plt.grid()
    plt.xlabel(r'$\lambda$ / $\mu$m',loc='right'),plt.ylabel('index of refraction') 

    plt.subplot(222),plt.plot(x,D),plt.plot(x,line0,'r--'),plt.grid()
    plt.title(r'$D=-\frac{\lambda}{c}*\frac{d^2n}{d\lambda^2}$')
    plt.xlabel(r'$\lambda$ / $\mu$m',loc='right'),plt.ylabel('Dispersion parameter D') 
    plt.text(1,-1e-11,'Normal dispersion',fontsize=8)
    plt.text(1.6,9e-12,'Anomalous dispersion',fontsize=8)

    plt.subplot(224),plt.plot(x,DeltaTao(wavelength=x,D=D,D0=x[p1]))
    plt.title(r'$\Delta\tau = D*\Delta\lambda$'),plt.grid()
    plt.xlabel(r'$\lambda$ / $\mu$m',loc='right'),plt.ylabel(r'Normalized delay difference $\Delta\tau$') 
    plt.show()

if __name__ == '__main__':
    main()
