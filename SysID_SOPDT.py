# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 13:59:32 2017

@author: Peter Nachtwey
Delta Computer Systems, Inc.
"""
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.integrate import odeint
from readCSV import readCSV


def difeq(y, t, k, t0, t1, c, dt ):
    """ generate estimated SOPDT solution
        y[0] = process value
        y[1] = rate of change of the process value"""
    t = max(t-dt,0)
    _u = control_interp(t)              # offset CO for dead time
    _dy2dt = (-(t0+t1)*y[1]-y[0]+k*_u+c)/(t0*t1)    # SOPDT dif Eq
    return np.array([y[1], _dy2dt])


def t1p2(p, aTime, aPV):
    """ find the values of a1 and a2 that minimize the ITAE3 """
    _k = p[0]                           # open loop extend gain
    _t0 = p[1]
    _t1 = p[2]
    _c = p[3]                           # output offset or bias
    _dt = p[4]
    _pv0 = [aPV[0], 0.0]                # initial process value and rate
    _aEV = odeint(difeq, _pv0, aTime, args=(_k, _t0, _t1, _c, _dt))
    _sse = np.sum((aPV-_aEV[:,0])**2)    
    # print("sse = {}".format(_sse))
    # plot_data(aTime, aPV, _aEV[:,0], aCO)
    return _sse 


def plot_data(aTimes, aPV, aEV, aCO):
    _fig, _ax0 = plt.subplots()
    _fig.set_size_inches(9.0, 6.0)
    _line0, = _ax0.plot(aTimes, aPV,
                      'c-', label='process variable')
    _line1, = _ax0.plot(aTimes, aEV,
                      'r--', label='estimated value')
    _ax0.set_title('Process and Estimated Values vs Time')
    _ax0.set_xlabel('time')             # units are data dependent 
    _ax0.set_ylabel('process and estimated values')

    _ax2 = _ax0.twinx()
    _line2, = _ax2.plot(aTimes, aCO ,'g-',label='control %')
    _ax2.set_ylabel('control %')
    if min(aCO) < 0:
        _ax2.set_ylim(-100.0, 100.) 
    else:        
        _ax2.set_ylim(0.0, 100.) 

    _lines = [_line0, _line1, _line2]
    _ax0.legend(_lines, [l.get_label() for l in _lines], loc='best')
    _fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    # path = "C://Users//peter//Desktop//Valve1.txt"
    path="c:\\Users\peter\Downloads\System RS Mode 1.csv"
    path="C:\\Users\Peter\OneDrive\Documents\mcd\Sysid\hotrod.txt"
    [aTime, aCO, aPV] = readCSV(path,'\t')
    N = len(aTime)
    aEV = np.zeros((N,2))               # estimated PV and PV'
    control_interp = interp1d(aTime, aCO, kind='linear',
                          bounds_error=False, fill_value='extrapolate')
    # initial guesses for the parameter array
    pv0 = [ 2,          # open loop gain.  PV change / %control output
           10,          # time constant 0
           10,          # time constant 1
           20,          # PV offset, ambient PV
           0]           # dead time 
    res = minimize(t1p2, pv0, args=(aTime, aPV), method='Nelder-Mead')
    pv0 = res.x         # do again to avvoid local minimum 
    res = minimize(t1p2, pv0, args=(aTime, aPV), method='Nelder-Mead')
    print(res)
    k = res.x[0]        # open loop gain.  PV change / %control output
    t0 = res.x[1]       # time constant 0
    t1 = res.x[2]       # time constant 1
    c = res.x[3]        # PV offset, ambient PV
    dt = res.x[4]       # dead time
    # initial process value and rate of change
    pv0 = [aPV[0], (aPV[1]-aPV[0])/(aTime[1]-aTime[0])]
    aEV = odeint(difeq, pv0, aTime, args=(k, t0, t1, c, dt ))
    plot_data(aTime, aPV, aEV[:,0], aCO)
    print("RMS error          = {:7.3f}".format(sqrt(res.fun/len(aTime))))
    print("The open loop gain = {:7.3f} PV/%CO".format(res.x[0]))
    print("Time constant 0    = {:7.3f}".format(res.x[1]))
    print("Time constant 1    = {:7.3f}".format(res.x[2]))
    print("Ambient PV         = {:7.3f} in PV units".format(res.x[3]))
    print("Dead time          = {:7.3f}".format(res.x[4]))
    print("Time units are the same as provided in input file")
    # calculate the controller ISA PID parameters
    tc = max(0.1*max(t0,t1),0.8*dt)     # closed loop time constant
    kc = (t0+t1)/(k*(tc+dt))            # controller gain %CO/error
    ti = t0+t1                          # integrator time constant 
    td = t0*t1/(t0+t1)                  # derivative time constant
    print("The closed loop time constant = {:7.3f}".format(tc))
    print("The controller gain           = {:7.3f} %CO/unit of error"
          .format(kc))
    print("The integrator time constant  = {:7.3f}".format(ti))
    print("The derivative time constant  = {:7.3f}".format(td))
    