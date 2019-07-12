import numpy as np
import numpy
from matplotlib.pyplot import figure, rc, axes, gca, plot, xlabel, ylabel, title, grid, savefig, show, subplots_adjust, axis, annotate, imshow, colorbar, close, xticks, yticks 
from numpy.fft.helper  import fftshift
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy import *
from scipy import stats
import matplotlib.ticker as plticker
rc('font',**{'family':'sans-serif','sans-serif':['Times New Roman']})

def baseline_als(y, lam, p, niter=10):
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

names = ['SDRspec_T_10.0_Nspectra_10_mjd_58450.5879774.npz', 'SDRspec_T_10.0_Nspectra_10_mjd_58450.5915733.npz', 'SDRspec_T_10.0_Nspectra_10_mjd_58450.6010006.npz', 'SDRspec_T_10.0_Nspectra_10_mjd_58450.603365.npz', 'SDRspec_T_10.0_Nspectra_10_mjd_58450.607042.npz', 'SDRspec_T_10.0_Nspectra_10_mjd_58450.5957141.npz', 'SDRspec_T_10.0_Nspectra_10_mjd_58450.6103586.npz']
lam = [9,7.8,7.8,7.8,7.8,8.8,7.9]
p = [0.03,0.02,0.02,0.02,0.02,0.03,0.01]

doppcor = [-0.255078,-0.191443,2.560924,5.293835,7.98487,-0.1212,2.108824]

ksamp = 1024
MHz = 1.e6
GHz = 1.e9
fsample = 2048000

for i in range(len(names)):
    xxx=load(names[i])
    tpower_total=xxx['tpower_total']
    spec_array = xxx['spec_array']
    elapsed_time_vec = xxx['elapsed_time_vec']
    MJD_vec = xxx['MJD_vec']
    RF_MHz = xxx['RF_MHz'].item()
    BW_MHz = xxx['BW_MHz'].item()
    T_per_spectrum = xxx['T_per_spectrum'].item()
    Nfft = xxx['Nfft'].item()
    Nspectra = xxx['Nspectra'].item()
    gain = xxx['gain'].item()
    channel_bandwidth_kHz = xxx['channel_bandwidth_kHz'].item()
    sample_rate = xxx['sample_rate_Hz'].item()
    pointing_string = xxx['pointing_string'].item()
    RF_Hz = 1.e6 * RF_MHz
    df = fsample / Nfft
    freqbins = arange(Nfft) - Nfft/2
    freqvec = (RF_MHz + freqbins * df) / MHz

    Nblocks_per_spectrum = int((T_per_spectrum * sample_rate) / Nfft)
    specave = average(spec_array,axis=0)
    specave[Nfft//2] = 0.5*(specave[Nfft//2-1]+specave[Nfft//2+1])
    for n in range(Nspectra):
        spec_array[n, Nfft//2] = 0.5*(spec_array[n,Nfft//2-1] + spec_array[n,Nfft//2+1])

    speclog = log10(specave)
    poly = baseline_als(speclog, lam=10.**lam[i], p=p[i], niter=15)
    sub_spec = speclog-poly
    corrspec = sub_spec/poly
    specsmooth = smooth(corrspec)
    new_freq = 1420.00026+freqvec
    print(new_freq)

    vel_l = []
    for y in range(len(new_freq)):
        vel = ((299792.*(1420.00026-new_freq[y]))/1420.00026) - doppcor[i]
        vel_l.append(vel)

    x = linspace(vel_l[0],vel_l[-1],len(specsmooth))
    print(names[i])
    print(pointing_string[11:])
    #print(vel_l)

    fig=figure()
    subplots_adjust(left=0.15)
    ax=fig.add_subplot(111)
    xlabel(r'$\rm Velocity\ (km\ /\ s)$')
    ylabel(r'$\rm \frac{T(\nu)}{T_{sys}}$ (K)')
    #ylabel(r'$\rm \log_{10} Spectrum\ (arbitrary\ units)$')
    #title(r'$\rm Spectrum\ vs\ Hanning\ Smoothed\ Spectrum$')
    title(r'$\rm Spec\ $'+str(pointing_string[11:]))
    annotate('%s'%names[i], xy=(0.65, 0.025), xycoords='figure fraction', ha='left', va='center', fontsize=7)
    #xticks((freqvec[0], 3*freqvec[0]/4, freqvec[0]/2, freqvec[0]/4, 0, freqvec[-1]/4, freqvec[-1]/2,
    #       3*freqvec[-1]/4, freqvec[-1]))
    plot(x,specsmooth,'r')
    show()









