"""
Using Astropy 
Reader4
Start with ALFALFA data in an ASCII Table.
ASCII format: (f8.1,1x,f12.6,1x,f10.4,1x,f10.4)
This table has four columns: heliocentric velocity [km/s], radio 
frequency [MHz], flux density [mJy], and a subtracted polynomial baseline
fit.

The point of this program is to plot the spectrum then allow the user to 
perform the measurement of velocity width, HI flux, and redshift of the detection.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy import stats
from matplotlib.widgets import Cursor
from termcolor import colored
from numpy import trapz
import operator

clicks = []   # to store mouse click positions

# Methods & Classes:
def smooth(x,window_len=9,window='hanning'):
    """smooth the data using a window with requested size. 
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
        flat window will produce a moving average smoothing.
    output:
        the smoothed signal    
    options:
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    """
    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

class Click():
    def __init__(self, ax, func, button=1):
        self.ax = ax 
        self.func = func
        self.button = button
        self.press = False
        self.move = False
        self.c1=self.ax.figure.canvas.mpl_connect('button_press_event', self.onpress)
        self.c2=self.ax.figure.canvas.mpl_connect('button_release_event', self.onrelease)
        self.c3=self.ax.figure.canvas.mpl_connect('motion_notify_event', self.onmove)
    def onclick(self,event):
        if event.inaxes == self.ax:
            if event.button == self.button:
                self.func(event, self.ax)
    def onpress(self,event):
        self.press = True
    def onmove(self,event):
        if self.press:
            self.move = True
    def onrelease(self,event):
        if self.press and not self.move:
            self.onclick(event)
        self.press = False; self.move = False

def func(event,ax):
    print(event.xdata,event.ydata)
    ax.scatter(event.xdata,event.ydata)
    clicks.append((event.xdata,event.ydata))
    ax.figure.canvas.draw()


# DATA HANDLING
var = input("Please enter a file name: ")
print("You chose " + var)

f = open(var,'r')
#f = open('A131378.txt','r')
header = []
v_helio = []
rf = []
flux = []
baseline = []
for i in range(30):
	header_i = f.readline()
	header_i = header_i.strip()
	header.append(header_i)
for line in f:
	line = line.strip()
	columns = line.split()
	v_helio.append(float(columns[0]))
	rf.append(float(columns[1]))
	flux.append(float(columns[2]))
	baseline.append(float(columns[3]))
f.close()

# Smoothing, default window size is 9
var1 = input('Do you wish to Hanning smooth? (y/n) ')

if var1 == 'y':
    var2 = input('Window size? (odd integer, default is 9) ')
    if not var2:
        flux = smooth(flux)
    else:
        flux = smooth(flux,window_len=int(var2))
    v_helio = np.linspace(v_helio[0],v_helio[-1],len(flux))

# Create points class with (v_helio, flux)
points = []
for i in range(len(v_helio)):
	point = (v_helio[i],flux[i])
	points.append(point)

# PLOTTING
fig, ax = plt.subplots(figsize=(10, 5))
ax.grid()
ax.set_xlabel('Heliocentric Velocity [km/s]')
ax.set_xlim([v_helio[0],v_helio[-1]])
ax.invert_xaxis()
ax.set_ylabel('Flux Density [mJy]')
ax.set_title('Baseline Subtracted Spectrum')
ax.plot(v_helio,flux,c='r',lw=0.9)
ymin, ymax = ax.get_ylim()
click = Click(ax, func, button=1)
plt.draw()
print(colored('\nZoom to fit. \n','green'))
print(colored('Select left edge of emission. \n','green'))
plt.pause(1) 
input("Hit Enter To Close\n")
plt.close(fig)

xs = []
ys = []
if len(clicks) == 2:
    for i in range(len(points)):
        if points[i][0] > clicks[0][0] and points[i][0] < clicks[1][0]:
            xs.append(points[i][0])
            ys.append(points[i][1])
    slope, intercept, r_value, p_value, std_err = stats.linregress(xs,ys)
    line = slope*np.asarray(xs)+intercept
X = np.linspace(v_helio[0],v_helio[-1],num=len(v_helio))
line = slope*np.asarray(X)+intercept

# Second plot that shows left edge fit. 
fig, ax = plt.subplots(figsize=(10, 5))
ax.grid()
ax.set_xlabel('Heliocentric Velocity [km/s]')
ax.set_xlim([v_helio[0],v_helio[-1]])
ax.set_ylim(ymin,ymax)
ax.invert_xaxis()
ax.set_ylabel('Flux Density [mJy]')
ax.set_title('Baseline Subtracted Spectrum')
ax.plot(v_helio,flux,c='r',lw=0.9)
ax.plot(X,line,c='k',lw=1)
click = Click(ax, func, button=1)
print(colored('Zoom to fit. \n','green'))
print(colored('Select right edge of emission. \n','green'))
plt.pause(1)
input("Hit Enter To Close\n")
plt.close(fig)



x1s = []
y1s = []
if len(clicks) == 4:
    for i in range(len(points)):
        if points[i][0] > clicks[2][0] and points[i][0] < clicks[3][0]:
            x1s.append(points[i][0])
            y1s.append(points[i][1])
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x1s,y1s)
    line1 = slope1*np.asarray(x1s)+intercept1
X1 = np.linspace(v_helio[0],v_helio[-1],num=len(v_helio))
line1 = slope1*np.asarray(X1)+intercept1


fig, ax = plt.subplots(figsize=(10, 5))
ax.grid()
ax.set_xlabel('Heliocentric Velocity [km/s]')
ax.set_xlim([v_helio[0],v_helio[-1]])
ax.set_ylim(ymin,ymax)
ax.invert_xaxis()
ax.set_ylabel('Flux Density [mJy]')
ax.set_title('Baseline Subtracted Spectrum')
ax.plot(v_helio,flux,c='r',lw=0.9)
ax.plot(X,line,c='k',lw=1)
ax.plot(X1,line1,c='k',lw=1)


plt.pause(1)
input("Hit Enter To Close\n")
plt.close(fig)

left_bound = -intercept/slope
right_bound = -intercept1/slope1
zone_x = []
zone_y = []
for i in range(len(points)):
    if points[i][0] > left_bound and points[i][0] < right_bound:
        zone_x.append(points[i][0])
        zone_y.append(points[i][1])
"""
Velocity width considerations...

max_i, max_val = max(enumerate(zone_y), key = operator.itemgetter(1))
HM = max_val/2.
nearest = (np.abs(zone_y - HM)).argmin()
"""

HI_flux = trapz(zone_y,dx=5)
print(colored('HI Flux is ~','green')+colored(str(HI_flux),'green')+colored(' [mJy km/s]','green'))


"""
Things that need to be implemented:
Interactive updating plots not opening & reopening
MAINTAIN ZOOM LEVEL
Options to do other kinds of smoothing.
Writing the analyzed properties to header.
Errors with estimation.
Gaussian support.

FWHM support.

vr = c(1-ν/ν0)
D = vr/H0

"""
