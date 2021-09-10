import cv2
import mrcfile
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


def read_stk(path):
    mrc = mrcfile.open(path)
    size = mrc.data.shape
    data = []
    for number_channel in range(size[0]):
        data.append(mrc.data[number_channel, :, :])
    return data


"""Usage: show power spectrum distribution of image."""
if __name__ == '__main__':
    path = '/path/to/data'
    s_group = path.split('.')
    f_type = s_group[len(s_group) - 1]
    if f_type == 'png' or f_type == 'jpg':
        data = cv2.imread(path)
        data = cv2.cvtColor(data, cv2.COLOR_RGB2GRAY)
        data = np.array(data, dtype=np.float32)
    elif f_type == 'mrcs':
        data = read_stk(path)
        data = data[0]
    else:
        raise Exception("Error: file format is not supported.")

    f = np.abs(np.fft.fftshift(np.fft.fft2(data)))
    f = f ** 2
    f = np.log(f)

    fig = plt.figure()
    ax = Axes3D(fig)

    y = x = np.linspace(0, data.shape[0] - 1, data.shape[0])
    X, Y = np.meshgrid(x, y)
    plt.title('demo')
    plt.xlabel('X')
    plt.ylabel('Y')
    """Use the below line to specify the perspective angle you like."""
    # ax.view_init(0, 180)
    ax.plot_surface(X, Y, f, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    plt.legend()
    plt.show()
