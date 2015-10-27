import numpy as np
from Chandra.Time import DateTime
import plot_aimpoint

# Get 99th percential absolute pointing radius

plot_aimpoint.opt = plot_aimpoint.get_opt()
asols = plot_aimpoint.get_asol()
# Last six months of data
asols = asols[asols['time'] > DateTime(-183).secs]
# center of box of range of data
mid_dy = (np.max(asols['dy']) + np.min(asols['dy'])) / 2.
mid_dz = (np.max(asols['dz']) + np.min(asols['dz'])) / 2.
# radius of each delta in mm (asol dy dz in mm)
dr = np.sqrt((asols['dy'] - mid_dy) ** 2 + (asols['dz'] - mid_dz) ** 2)
dr_99 = np.percentile(dr, 99)
dr_99_arcsec = dr_99 * 20
print "99th percentile radius of 6m data is {} arcsec".format(dr_99_arcsec)
