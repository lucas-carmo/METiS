# spell-checker: disable
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from IPython import get_ipython


class mtsOut(object):
    def __init__(self, flPath):
        self.path = flPath
        data = pd.read_csv(self.path, delim_whitespace = True)
        data.columns = map(str.lower, data.columns)
        for label in data.columns:
            setattr(self, label, getattr(data, label))

# Example:
# mts = mtsOut('./test/verifyAero/verifyAero_case1_out.txt')
# for ii in range(1, 7): 
#     ax = plt.subplot(2,3,ii)
#     ax.plot(mts.time, getattr(mts, 'ad_hub_force_' + str(ii)))
#     ax.set_xlabel('Time (s)')
#     ax.set_ylabel('Aero load - ' + str(ii))
# plt.show()