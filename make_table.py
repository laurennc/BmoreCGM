import numpy as np


def make_cooling_table():
	hden_n_bins, hden_min, hden_max = 40, -6, 1.8
	T_n_bins, T_min, T_max = 51, 3, 8
	patt = '/Users/tardis/data/CoolingMaps/CoolingMapRun/g1q1_run%i'+'.dat'
	hden=np.linspace(hden_min,hden_max,hden_n_bins)
        T=np.linspace(T_min,T_max, T_n_bins)
        table = np.zeros((hden_n_bins,T_n_bins))
        for i in range(hden_n_bins):
                table[i,:]=[float(l.split()[2]) for l in open(patt%(i+1)) if l[0] != "#"]
        return hden,T,table



