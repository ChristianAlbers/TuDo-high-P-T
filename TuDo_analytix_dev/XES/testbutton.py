from matplotlib import pyplot as plt
from matplotlib.widgets import Button

def handler(*args, **kwargs):
    print('handled')

def testfn():
    f1 = plt.figure('f1')
    b1 = Button(f1.add_axes([0.4, 0.3, 0.1, 0.04]), 'Click!')
    b1.on_clicked(handler)
    return b1

f2 = plt.figure('f2')
b2 = Button(f2.add_axes([0.4, 0.3, 0.1, 0.04]), 'Click!')
b2.on_clicked(handler)

button = testfn()

plt.show()