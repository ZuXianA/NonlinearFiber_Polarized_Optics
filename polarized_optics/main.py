import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['agg.path.chunksize'] = 10000  # 增加块大小


class PolarizationTrace(object):
    def __init__(self, max_theta):
        self.max_theta = max_theta

    def __ROT(self, theta):  # 定义一个旋转器，庞加莱球上是绕着 S3 右旋 2*theta 角
        return np.array([[np.cos(np.deg2rad(2 * theta)), -np.sin(np.deg2rad(2 * theta)), 0],
                         [np.sin(np.deg2rad(2 * theta)), np.cos(np.deg2rad(2 * theta)), 0],
                         [0, 0, 1]])

    def __ROTINV(self, mat):  # 定义一个旋转器的逆，用于计算米勒矩阵 R
        return np.linalg.inv(mat)

    def __RWP(self, delta):  # 定义零方位角相位延迟器，庞加莱球上是绕着 S1 右旋 delta 角
        return np.array([[1, 0, 0],
                         [0, np.cos(np.deg2rad(delta)), -np.sin(np.deg2rad(delta))],
                         [0, np.sin(np.deg2rad(delta)), np.cos(np.deg2rad(delta))]])

    def SingleWaveplate(self, state, waveplate='QWP'):
        if waveplate.lower() == 'qwp':
            wp = self.__RWP(90)  # 1/4 波片
        elif waveplate.lower() == 'hwp':
            wp = self.__RWP(180)  # 1/2 波片
        else:
            print('Error, please re-input waveplate.')
            return None

        Out = []
        theta_range = np.arange(0, self.max_theta, 1)

        for i in theta_range:
            rot = self.__ROT(i)
            res = rot @ wp @ self.__ROTINV(rot)
            res[np.isclose(a=res, b=0, atol=1.3e-16)] = 0

            Out.append(res @ state)

        return np.array(Out)

    def TwoWaveplate(self, state):
        Out = []
        qwp = self.__RWP(90)  # 1/4 波片
        hwp = self.__RWP(180)  # 1/2 波片
        theta_range = np.arange(0, self.max_theta, 1)

        for i in theta_range:
            rot1 = self.__ROT(i)
            for j in theta_range:
                rot2 = self.__ROT(j)
                res = rot2 @ hwp @ rot1 @ qwp @ self.__ROTINV(rot1) @ self.__ROTINV(rot2)
                res[np.isclose(a=res, b=0, atol=1.3e-16)] = 0

                Out.append(res @ state)

        return np.array(Out)

    def ThreeWaveplate(self, state):
        Out = []
        qwp1 = self.__RWP(90)  # 1/4 波片
        hwp = self.__RWP(180)  # 1/2 波片
        qwp2 = self.__RWP(90)  # 1/4 波片
        theta_range = np.arange(0, self.max_theta, 1)

        for i in theta_range:
            rot1 = self.__ROT(i)
            for j in theta_range:
                rot2 = self.__ROT(j)
                for k in theta_range:
                    rot3 = self.__ROT(k)
                    res = rot3 @ qwp2 @ rot2 @ hwp @ rot1 @ qwp1 @ self.__ROTINV(rot1) @ self.__ROTINV(
                        rot2) @ self.__ROTINV(rot3)
                    res[np.isclose(a=res, b=0, atol=1.3e-16)] = 0

                    Out.append(res @ state)
        return np.array(Out)


def PlotShow(dat, S, viz=True):
    name = '{}_{}_{}'.format(str(S[0]), str(S[1]), str(S[2]))
    # 绘制Poincaré球上的轨迹
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(dat[:, 0], dat[:, 1], dat[:, 2], color='blue', linewidth=2, label=name)
    ax.scatter(dat[0, 0], dat[0, 1], dat[0, 2], color='red', s=30, label='Start')
    ax.scatter(dat[-1, 0], dat[-1, 1], dat[-1, 2], color='green', s=30, label='End')
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_xlabel('S1')
    ax.set_ylabel('S2')
    ax.set_zlabel('S3')
    ax.legend()
    plt.tight_layout()
    plt.savefig('qwp+hwp+qwp_' + name, dpi=300)

    if viz:
        plt.show()


def demo1():
    # 三个输入偏振态
    S1 = np.array([1, 0, 0])
    S2 = np.array([0, 1, 0])
    S3 = np.array([0, 0, 1])

    pl = PolarizationTrace(max_theta=180)
    res1 = pl.SingleWaveplate(state=S1, waveplate='qwp')
    res2 = pl.SingleWaveplate(state=S2, waveplate='qwp')
    res3 = pl.SingleWaveplate(state=S3, waveplate='qwp')

    PlotShow(dat=res1, S=S1, viz=False)
    PlotShow(dat=res2, S=S2, viz=False)
    PlotShow(dat=res3, S=S3, viz=False)


def demo2():
    # 三个输入偏振态
    S1 = np.array([1, 0, 0])
    S2 = np.array([0, 1, 0])
    S3 = np.array([0, 0, 1])

    pl = PolarizationTrace(max_theta=180)
    res1 = pl.SingleWaveplate(state=S1, waveplate='hwp')
    res2 = pl.SingleWaveplate(state=S2, waveplate='hwp')
    res3 = pl.SingleWaveplate(state=S3, waveplate='hwp')

    PlotShow(dat=res1, S=S1, viz=False)
    PlotShow(dat=res2, S=S2, viz=False)
    PlotShow(dat=res3, S=S3, viz=False)


def demo3():
    # 三个输入偏振态
    S1 = np.array([1, 0, 0])
    S2 = np.array([0, 1, 0])
    S3 = np.array([0, 0, 1])

    pl = PolarizationTrace(max_theta=180)
    res1 = pl.TwoWaveplate(state=S1)
    res2 = pl.TwoWaveplate(state=S2)
    res3 = pl.TwoWaveplate(state=S3)

    PlotShow(dat=res1, S=S1, viz=False)
    PlotShow(dat=res2, S=S2, viz=False)
    PlotShow(dat=res3, S=S3, viz=False)


def demo4():
    # 三个输入偏振态
    S1 = np.array([1, 0, 0])
    S2 = np.array([0, 1, 0])
    S3 = np.array([0, 0, 1])

    pl = PolarizationTrace(max_theta=180)
    res1 = pl.ThreeWaveplate(state=S1)
    res2 = pl.ThreeWaveplate(state=S1)
    res3 = pl.ThreeWaveplate(state=S1)

    PlotShow(dat=res1, S=S1, viz=False)
    PlotShow(dat=res2, S=S2, viz=False)
    PlotShow(dat=res3, S=S3, viz=False)


def main():
    # Out1 = Trace(state=S3)
    # PlotShow(dat=Out1, S=S3, viz=False)
    # demo4()
    S1 = np.array([1, 0, 0])
    pl = PolarizationTrace(max_theta=180)
    res1 = pl.ThreeWaveplate(state=S1)
    PlotShow(dat=res1, S=S1, viz=False)


if __name__ == '__main__':
    main()
