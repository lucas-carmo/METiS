import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from IPython import get_ipython

get_ipython().run_line_magic('matplotlib', 'qt')


class mtsOut(object):
    def __init__(self, paths):
        self.path = paths                
        self.data = [pd.read_csv(path, delim_whitespace = True) for path in paths]
        
        for data in self.data:
            data.columns = map(str.lower, data.columns)
        
    def motions(self):
        
        for data in self.data:
            ax = plt.subplot(2,3,1)
            ax.plot(data.time, data.surge)
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Surge (m)')
            
            ax = plt.subplot(2,3,2)
            ax.plot(data.time, data.sway)
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Sway (m)')            
            
            ax = plt.subplot(2,3,3)
            ax.plot(data.time, data.heave)
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Heave (m)')            
            
            ax = plt.subplot(2,3,4)        
            ax.plot(data.time, data.roll * np.pi / 180)
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Roll (deg)')            
            
            ax = plt.subplot(2,3,5)
            ax.plot(data.time, data.pitch * np.pi / 180)
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Pitch (deg)')            
    
            ax = plt.subplot(2,3,6)
            ax.plot(data.time, data.yaw * np.pi / 180)
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Yaw (deg)')            
        
        
        
flPaths = ["..\\test\\Dir45_H2p0_T30p00_out.txt", "C:\\Users\\lucas.henrique\\Google Drive\\Doutorado\\1Testes_OC4\\OC4_regular\\METiS\\Dir45_H2p0_T30p00_out.txt"];
        
mts = mtsOut(flPaths)
mts.motions()



