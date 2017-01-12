import yt
import numpy as np
import matplotlib.pyplot as plt

def starting_cloud_region(ds):
    
    return ds.sphere

def find_GasMass_dans_TempRange(region,low=5.0,high=1e7):
    idx = np.where( (data['temperature'] <= high) & (data['temperature'] >= low) )[0]
    return np.sum(data['CellMass'])



