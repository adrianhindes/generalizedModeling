import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button, RadioButtons


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')        
X= np.arange(-50,50,2)
Y=np.arange(-50,50,2)
X,Y = np.meshgrid(X,Y)


h0=5

Z = np.sqrt(((X-h0)**2+(Y-h0)**2)/(np.tan(np.pi/120)))            
l = ax.plot_surface(X,Y,Z, rstride=3, cstride=3)         
plt.axis('scaled')



axhauteur = plt.axes([0.2, 0.1, 0.65, 0.03])
shauteur = Slider(axhauteur, 'H', 0.5, 10.0, valinit=h0)

def update(val): 
    h = shauteur.val 
    ax.clear()
    l=ax.plot_surface(X,Y,np.sqrt(((X-h)**2+(Y-h)**2)/(np.tan(np.pi/120))),rstride=3, cstride=3)
    fig.canvas.draw_idle() 
    
shauteur.on_changed(update)

plt.show()