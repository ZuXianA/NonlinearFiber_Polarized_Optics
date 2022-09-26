import re
import requests
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt


def GetHtmlText(url):       # 获取指定url页面的网页源代码 因为数据就在其中 参见SiO2.html文件第790和792行
    try:
        r = requests.get(url=url,timeout=30)
        r.raise_for_status
        r.encoding = r.apparent_encoding
        r.close()
        return r.text
    except:
        print('产生异常')

def GetData(text):          # 利用正则表达式提取数据
    obj1 = re.compile(pattern=r'data_n_wl=(?P<wavelength>.*?);',flags=re.S)
    obj2 = re.compile(pattern=r'data_n=(.*?);',flags=re.S)

    wavelength = obj1.findall(string=text)          # 函数返回的是带括号的列表
    wavelength = ''.join(wavelength).strip('[]')    # 将列表转换为字符串并消除中括号
    wavelength = wavelength.split(',')              # 以逗号分割数据，方便后续转为 ndarry 类型
    
    value = obj2.findall(string=text)
    value = ''.join(value).strip('[]')
    value = value.split(',')
    return wavelength,value

def Sellmeier(x,B1,C1,B2,C2,B3,C3):         # 构造三阶表达式 共6个系数有待确定
    offset = 1                              # 此处需要传入波长的平方
    ssum = B1 * x / (x - C1**2) + B2 * x/ (x - C2**2) + B3 * x/ (x - C3**2)
    return ssum+offset

def main():
    url = r'https://refractiveindex.info/?shelf=main&book=SiO2&page=Malitson'
    html = GetHtmlText(url=url)

    wavelength,value = GetData(text=html)
    wavelengths = np.array(wavelength,dtype=np.float16)     #将列表转换为 ndarry 并定义数据类型（原始列表存的是字符串类型）
    refractions = np.array(value,dtype=np.float32)

    plt.figure(1)
    # 画出原始表达式对应的曲线
    refractions_index = Sellmeier(wavelengths**2,0.6961663,-0.0684043,0.4079426,-0.1162414,0.8974794,-9.896161)
    plt.scatter(wavelengths,np.sqrt(refractions_index),color='r',label="origin")
    plt.xlabel('wavelength / um'),plt.ylabel('index of refrection')
    #plt.title('original Sellmeier')

    
    # 开始非线性拟合 Ps.需要给系数确定上下限，否则得到的结果与原表达差别非常大 即 bounds 参数项
    popt,pcov = optimize.curve_fit(f=Sellmeier,xdata=wavelengths**2,ydata=refractions*2,bounds=([0.5,-0.07,0,-0.2,0.8,-10],[0.7,0,0.41,0,0.9,-5]))
    popt_test,pcov_test = optimize.curve_fit(f=Sellmeier,xdata=wavelengths**2,ydata=refractions*2)      # 此处不限定系数的上下限
    print('给参数确定上下限后多项式对应的系数为：',*popt)
    print('参数取值无限制多项式对应的系数为：',*popt_test)
    
    f = open(file='FitSellmeierBy101.txt',mode='w',encoding='utf-8')
    for d in popt:
        f.write(f'{d}\t')
    f.write('\n')
    for d in popt_test:
        f.write(f'{d}\t')    
    f.close()
    
    plt.scatter(x=wavelengths,y=np.sqrt(Sellmeier(wavelengths**2,*popt)),color='g',label='fit')
    plt.grid(),plt.legend()
    plt.show()

    plt.scatter(x=wavelengths,y=np.sqrt(Sellmeier(wavelengths**2,*popt_test)),color='b',label='fit*')
    


if __name__ == '__main__':
    main()