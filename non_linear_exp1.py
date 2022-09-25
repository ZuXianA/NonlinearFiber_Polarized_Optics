import re
import requests
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt


def GetHtmlText(url):
    try:
        r = requests.get(url=url,timeout=30)
        r.raise_for_status
        r.encoding = r.apparent_encoding
        r.close()
        return r.text
    except:
        print('产生异常')

def GetData(text):
    obj1 = re.compile(pattern=r'data_n_wl=(?P<wavelength>.*?);',flags=re.S)
    obj2 = re.compile(pattern=r'data_n=(.*?);',flags=re.S)

    wavelength = obj1.findall(string=text)
    wavelength = ''.join(wavelength).strip('[]')    # 字符串
    wavelength = wavelength.split(',')              # 列表
    
    value = obj2.findall(string=text)
    value = ''.join(value).strip('[]')
    value = value.split(',')
    return wavelength,value

def Sellmeier(x,B1,C1,B2,C2,B3,C3):
    offset = 1
    
    ssum = B1 * x / (x - C1**2) + B2 * x/ (x - C2**2) + B3 * x/ (x - C3**2)
    return ssum+offset

def main():
    url = r'https://refractiveindex.info/?shelf=main&book=SiO2&page=Malitson'
    html = GetHtmlText(url=url)

    wavelength,value = GetData(text=html)
    wavelengths = np.array(wavelength,dtype=np.float16)
    refractions = np.array(value,dtype=np.float32)

    # 画出原始表达式对应的曲线
    refractions_index = Sellmeier(wavelengths**2,0.6961663,-0.0684043,0.4079426,-0.1162414,0.8974794,-9.896161)
    plt.scatter(wavelengths,np.sqrt(refractions_index),color='r',label="origin")
    plt.xlabel('wavelength / um'),plt.ylabel('index of refrection')
    plt.title('original Sellmeier')

    # 开始非线性拟合
    popt,pcov = optimize.curve_fit(f=Sellmeier,xdata=wavelengths**2,ydata=refractions*2,bounds=([0.5,-0.07,0,-0.2,0.8,-10],[0.7,0,0.41,0,0.9,-5]))
    plt.scatter(x=wavelengths,y=np.sqrt(Sellmeier(wavelengths**2,*popt)),color='g',label='fit')
    plt.grid(),plt.legend()
    plt.show()
    
    print('多项式对应的系数为：',*popt)

'''    with open(file='sio2.html',mode='w',encoding='utf-8') as f:
        f.write(html)
'''
if __name__ == '__main__':
    main()