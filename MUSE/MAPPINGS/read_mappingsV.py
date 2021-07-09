# @Author: mparra
# @Date:   06-07-2021
# @Email:  maxime.parra@irap.omp.eu
# @Last modified by:   mparra
# @Last modified time: 09-07-2021


import os,sys
sys.path.append('/home/mparra/PARRA/Scripts/Python/MUSE/MAPPINGS/')
sys.path.append('/home/mparra/PARRA/Scripts/Python/MUSE/')
sys.path.append('/home/mparra/PARRA/Scripts/Python/MUSE/BPT/')
import argparse
import numpy as np

from mpdaf.obj import Image

import bpt_config as bptcfg
from bpt_utils import bpt_single
from mappings_utils import shock_model


'''
This script compares the current results of the MAPPINGSV shock table database to BPT diagrams

for now we keep a single model (default: Gutkin) and a single density (default: 1/cm^3) per iteration to ensure visibility. 

The Gutkin model is the only one available for now (see mappings_utils).

The carbon to dust ratio for the Gutkin model is fixed to high.

The range of shock speed displayed can be modified. The default is 100-250 km/s
Note that the Gutkin model ranges from 100 to 1000 km/s in steps of 25.
'''

# plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')

ap = argparse.ArgumentParser(description='Script to read and plot the output of the mappings V libraries (Allen et al. 2008)')
ap.add_argument("-m", "--models", nargs='*', help="Model to be compared to the data", default='Gutkin')
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='mappings', type=str)
ap.add_argument("-r", "--regions", nargs='*', help="Region files to be overlaid on the image", default=None)
ap.add_argument("-c", "--contours", nargs='?', help="Fits file to use to overlay contours", type=str)
ap.add_argument("-d", "--density", nargs='?', help="density to use in per cubic cm", type=float,default=1)
ap.add_argument("-vmin", "--minspeed", nargs='?', help="minimal shock speed velocity to be displayed", type=float,default=100)
ap.add_argument("-vmax", "--maxspeed", nargs='?', help="max shock speed velocity to be displayed", type=float,default=250)
args = ap.parse_args()

#renaming arguments
outdir=args.outdir
if not os.path.isdir(outdir):
    os.mkdir(outdir)

model= args.models
dens=args.density

bpt_type = bptcfg.type
lineratiomaps = [bptcfg.lineratio_paths["n2_ha"], bptcfg.lineratio_paths["s2_ha"], bptcfg.lineratio_paths["oI_ha"]]

#range of speeds used
vmin=args.minspeed
vmax=args.maxspeed
speeds = np.linspace(vmin,vmax,num=int((vmax-vmin)/25)+1)

#creating the models
shock_tables=shock_model(model)

#getting the abundances and their solar fraction
abunds=shock_tables.abund_float
abunds_frac=abunds/0.01524

# use_pos=np.where((abunds_frac<=1) & (abunds_frac>0.01))
#instead of the line above we restrict manually to a few values for less overlap
use_pos=np.array([4,7,10])
abunds_use=abunds[use_pos]

#creating the labels which we will use for the grids
labels_use0=np.array([r"$Z_{ISM}$ = "]*len(abunds_use))
labels_use1=np.char.add(abunds_use.astype(str),np.array([' i.e. ']*len(abunds_use)))
labels_use2=np.char.add(abunds_frac.round(decimals=3)[use_pos].astype(str),np.array([r" $Z/Z_\odot$"]*len(abunds_use)))
labels_use=np.char.add(labels_use1,labels_use2)
labels_use=np.char.add(labels_use0,labels_use)

#putting the models of interest in an array
mod_use=np.array([None]*len(abunds_use))
for i in range(len(abunds_use)):
    #getting the the 2d grids for the abundance, with the correct density
    mod_use[i]=shock_tables.abund_single(abunds_use[i]).ctd_high.d_lgrid(dens)
    
    #restricting the shock velocity to the ones inputted
    mod_use[i]=mod_use[i][mod_use[i]['shck_vel'].isin(speeds)]
    
title="BPT compared to a range of "+model+" shock models, with n="+str(dens)+" cm$^{-3}$ and shock speeds in ["+str(vmin)+","+str(vmax)+"] km/s"

#computations
ymap = Image(bptcfg.lineratio_paths["o3_hb"])
logy = np.log10(ymap.data.data)
colormap = ["Blues", "Greens", "Purples", "Oranges"]

for i, ratiomap in enumerate(lineratiomaps):
    
    bpt_type = i + 1
    bpt_fits, bpt_indexes, figure = bpt_single(ratiomap, logy, args.regions, args.contours, bpt_type, colormap,title=title,
                                                shock_mods=mod_use,mod_labels=labels_use,figsize='big')

    outfile = "%s/mapV_bpt%d_d%d.png" % (outdir, bpt_type,dens)
    figure.savefig(outfile, format="png", bbox_inches="tight", pad_inches=0.4)
    print("Results stored to %s" % outfile)
