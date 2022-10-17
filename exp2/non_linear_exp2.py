import numpy as np
import matplotlib.pyplot as plt

c = 299792458                       # 光速
D = 0.092 / 4 * (1550-1312**4/1550**3) * 1e-12 / 1e-9 / 1e3
beta2 = -1550e-9**2 / 2 / np.pi / c*D

start_length = 0
stop_length = 100
num_of_length = 11
L = np.linspace(start=start_length,stop=stop_length,num=num_of_length,dtype=np.int0)      # 光纤长度

'''定义离散网络'''
N = 2**8
twin = 50e-12
dt = twin/N
df = 1/twin
fwin = 1/dt
t = np.linspace(start=-twin/2,stop=twin/2-dt,num=N)
f = np.linspace(start=-fwin/2,stop=fwin/2-df,num=N)

'''define the input Gaussian/super-Gaussian pulse'''
FWHM = 1e-12                       # full width at half maxima
morder = 1                         # mth order super-Gaussian pulse; m=1 for Gaussian pulse
Nc = 2 

a = np.zeros([len(L),int(N)],dtype=complex)
A = np.zeros([len(L),int(N)],dtype=complex)
a0 = np.zeros([len(L),int(N)],dtype=complex)

for row in range(len(L)):
    a[row,:] = np.exp(-np.log(2) / 2*np.abs(2*t/FWHM)**(2*morder)) * np.exp(-1j*2*np.log(2) / FWHM**2 * np.sqrt(Nc**2-1)* t**2)
    A[row,:] = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(a[row,:])))
    A[row,:] = A[row,:] * np.exp(1j * 0.5 * beta2 * L[row] * (2*np.pi*f)**2)
    a0[row,:] = np.fft.fftshift(np.fft.fft(np.fft.fftshift(A[row,:])))

ax = plt.figure().add_subplot(projection='3d')

for i in range(len(L)):                             #结合for循环绘制多张线图
    y=L[i]
    y=y*np.ones(256)                                #保证与x轴的点数一致，这一步非常重要
    ax.plot(t*1e12,y,a0[i,:]*np.conj(a0[i,:]))

#坐标及坐标轴相关设置
if Nc == 1:
    ax.set_title('Gaussian pulse (Chirp-free, Nc=1)',fontsize=18)
else:
    ax.set_title('Chirped Gaussian pulse (Nc=2)',fontsize=18)

ax.set_xlabel('time'),ax.set_ylabel('fiber length / m'),ax.set_zlabel('magnitude')
plt.show()