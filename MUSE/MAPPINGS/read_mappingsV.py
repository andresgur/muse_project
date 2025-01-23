# @Author: Andrés Gúrpide <agurpide>
# @Date:   07-06-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 26-01-2022

import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
import sys
# sys.path.append('/home/etudiant/Documents/GitHub/muse_project/MUSE/utils')
import muse_utils as mu
import matplotlib.patches as mpatches
# sys.path.append('/home/etudiant/Documents/GitHub/muse_project/MUSE/BPT')
from bpt.bpt import BPT_1, BPT_2, BPT_3
from configparser import ConfigParser, ExtendedInterpolation
from astropy.io import fits
import glob
# assuming a gas to dust ratio of 0.1
abundances_logOH = {"0d0001": 6.61, "0d0002": 6.91, "0d0005": 7.30, "0d001": 7.61, "0d002": 7.91, "0d004": 8.21, "0d010": 8.61, "0d01524": 8.83}

parser = argparse.ArgumentParser(description='Script to read and plot the output of the mappings III libraries (Allen et al. 2008)')
parser.add_argument("-a", "--abundances", nargs='+', help="Abundances", default=["0d01524"], type=str)
parser.add_argument("-d", "--densities", nargs='+', help="Densities", default=[1], type=float)
parser.add_argument("-m", "--model", nargs='?', help="Model", default="Gutkin16", choices=["Gutkin16", "Allen2008"])
parser.add_argument("-t", "--type", nargs='?', help="Shock, precursor or shock + precursors (s, p, s+p)", default="s+p", choices=["s", "p", "s+p"])
parser.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='MAPPINGS/', type=str)
parser.add_argument("--config", nargs=1, help="Input config file", required=True)
args = parser.parse_args()

style_file = '/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle'

if os.path.isfile(style_file):
    plt.style.use(style_file)

bptcfg = ConfigParser(interpolation=ExtendedInterpolation())
# bpt_type = bptcfg.type
bptcfg.read("%s" %args.config[0])
lineratio_paths = bptcfg["Paths"]
lineratiomaps = [lineratio_paths["n2_ha"], lineratio_paths["s2_ha"], lineratio_paths["oI_ha"]]
grid_figure, grid_axes = plt.subplots(1, 3, sharey=True, gridspec_kw={"wspace": 0.05}, figsize=(25.6, 9.8))
figure_bpt1, ax_1 = plt.subplots(1)
figure_bpt2, ax_2 = plt.subplots(1)
figure_bpt3, ax_3 = plt.subplots(1)
patches = []
#speeds = [100.0, 150.0, 200.0, 250.0, 300.0]
speeds = [100.0, 150, 200.0, 250, 300.0]
#speeds = [200.0, 225, 250.0, 300.0, 350, 400, 450, 500, 550, 600, 700, 800, 850, 900, 1000]
markers = mu.get_markers_array(len(speeds))
if args.type == "s":
    folder = "tables/shock"
elif args.type == "p":
    folder = "tables/precursor"
elif args.type == "s+p":
    folder = "tables/shock_plus_precursor"
print("Using %s models" % folder)
models = []
abund_labels = []

c_to_dust = "1d00"
# c_to_dust = "0d26"

if args.model=="Gutkin16":
    for abund in args.abundances:
        if abund not in abundances_logOH:
            print("Warning: abundance %s not recognized, skipping" % abund)
            continue

        files = glob.glob("%s/%s_ISM%s_C%s.php" % (folder, args.model, abund, c_to_dust))
        if len(files) ==0:
            print("%s/%s_ISM%s_C%s.php does not exist" % (folder, args.model, abund, c_to_dust))
            continue
        models.append(glob.glob("%s/%s_ISM%s_C%s.php" % (folder, args.model, abund, c_to_dust))[0])
        abund_labels.append("12 + log(O/H) = %.2f" % abundances_logOH[abund])
        if len(models) == 0:
            raise ValueError('No models found with the current inputs')

elif args.model=="Allen2008":

    for abund in args.abundances:
        models.append(glob.glob("%s/%s_%s.php" % (folder, args.model, abund))[0])
        abund_labels.append(abund)

colors = mu.create_color_array(len(models), "viridis")
grid_axes[0].set_ylabel("log([OIII]/H$_\\beta$)")
grid_axes[0].set_xlabel("log([NII]/H$_\\alpha$)")
grid_axes[1].set_xlabel("log([SII]/H$_\\alpha$)")
grid_axes[2].set_xlabel("log([OI]/H$_\\alpha$)")
ax_1.set_ylabel("log([OIII]/H$_\\beta$)")
ax_1.set_xlabel("log([NII]/H$_\\alpha$)")
ax_2.set_ylabel("log([OIII]/H$_\\beta$)")
ax_2.set_xlabel("log([SII]/H$_\\alpha$)")
ax_3.set_ylabel("log([OIII]/H$_\\beta$)")
ax_3.set_xlabel("log([OI]/H$_\\alpha$)")
outdir = args.outdir
if not os.path.exists(outdir):
    os.makedirs(outdir)
outdir = outdir + "abund"
outdir += "_".join(args.abundances)
outdir += "_dens"
for dens in args.densities:
    outdir += "_%.1f" % dens
outdir += "_%s_%s_C%s" % (args.model, args.type, c_to_dust)

linestyles = mu.get_linestyle_array(len(models))

for model, color, abund_label in zip(models, colors, abund_labels):

    print("Processing models %s: %s" % (model, abund_label))
    for dens, ls in zip(args.densities, linestyles):
        print("Reading file %s" % model)

        data = np.genfromtxt(model, names=True)
        data_dens = data[np.where(data["dens"] == dens)]
        if len(data_dens) == 0:
            print("Warning: model %s does not have values with density %.2f" % (model, dens))
            continue
        else:
            patches.append(mpatches.Patch(color=color, ls="solid", label="%s" % (abund_label)))
        ratios_y = []
        ratios_1 = []
        ratios_2 = []
        ratios_3 = []
        marker_norm = 5
        for speed, marker in zip(speeds, markers):
            data_shck = data_dens[np.where(data_dens["shck_vel"] == speed)]
            # take every two values of magnetic field
            data_shck = data_shck[np.arange(0, len(data_shck) - 1, 2)]

            #data_shck = data_shck[np.arange(0, len(data_shck) - 1, 2)]
            print("Magnetic values considered")
            print(data_shck["mag_fld"])
            oIII_hb = data_shck["OIII_Hb"]
            ratio_y = np.log10(oIII_hb)
            nII_ha = data_shck["NII_Ha"]
            ratio1 = np.log10(nII_ha)
            sII_ha = data_shck["SII_Ha"]
            ratio2 = np.log10(sII_ha)
            oI_ha = data_shck["OI_Ha"]
            ratio3 = np.log10(oI_ha)
            ratios_1.append(ratio1)
            ratios_y.append(ratio_y)
            ratios_2.append(ratio2)
            ratios_3.append(ratio3)
            for ax, grid_ax, ratio in zip([ax_1, ax_2, ax_3], grid_axes, [ratio1, ratio2, ratio3]):
                ax.plot(ratio, ratio_y, ls="--", color=color, lw=1, markersize=3, marker=".")
                ax.text(ratio[0], ratio_y[0], s="%d" % speed, fontsize=18, clip_on=True)
                grid_ax.plot(ratio, ratio_y, ls="--", color=color, lw=1, markersize=3, marker=".")
                grid_ax.text(ratio[0], ratio_y[0], s="%d" % speed, fontsize=22, clip_on=True)

        for ax, ratios in zip(grid_axes, [ratios_1, ratios_2, ratios_3]):
            ax.plot(ratios, ratios_y, color=color, ls=ls, lw=1)
        ax_1.plot(ratios_1, ratios_y, color=color, ls=ls, lw=1)
        ax_2.plot(ratios_2, ratios_y, color=color, ls=ls, lw=1)
        ax_3.plot(ratios_3, ratios_y, color=color, ls=ls, lw=1)
# plot BPT lines and data points
y_axis_map = fits.open(lineratio_paths["o3_hb"])
print("WARNING: values below -1.0 for  [OIII]/Hb are beign ignored")
y_axis_map[0].data[np.where(np.log10(y_axis_map[0].data) < -1.0)] = np.nan
logy = np.log10(y_axis_map[0].data.data)
for bpt_type, ax, map_1, grid_ax in zip([1, 2, 3], [ax_1, ax_2, ax_3], lineratiomaps, grid_axes):
    x_axis_map = fits.open(map_1)
    logx = np.log10(x_axis_map[0].data)
    # bpt_diagram = BPT_diagram(bpt)
    if bpt_type == 1:
        bpt_diagram = BPT_1()
    elif bpt_type==2:
        bpt_diagram = BPT_2()
    elif bpt_type==3:
        bpt_diagram = BPT_3()

    invalid_pixels = np.where((np.isnan(logx)) | (np.isnan(logy)))
    inter_starform = np.sort(np.reshape(logx, np.size(logx)))
    if grid_ax is not None:
        range = inter_starform[np.where(inter_starform < bpt_diagram.limit)]
    ax.plot(range, bpt_diagram.pure_starform_crv(range), color="black", zorder=100)

    grid_ax.plot(range, bpt_diagram.pure_starform_crv(range), color="black", zorder=100)
    inter_def = logy[logy > bpt_diagram.int_inter[0]]
    inter_def = np.sort(inter_def[inter_def < bpt_diagram.int_inter[1]])
    ax.plot(bpt_diagram.intermediate_crv(inter_def), inter_def, color='black', zorder=100)
    grid_ax.plot(bpt_diagram.intermediate_crv(inter_def), inter_def, color='black', zorder=100)
    agnliner_def = logx[logx > bpt_diagram.agnliner_inter[0]]
    agnliner_def = np.sort(agnliner_def[agnliner_def < bpt_diagram.agnliner_inter[1]])
    # zorder 100 so it stands above the points
    ax.plot(agnliner_def, bpt_diagram.agnliner_crv(agnliner_def), color='red', zorder=100)
    grid_ax.plot(agnliner_def, bpt_diagram.agnliner_crv(agnliner_def), color='red', zorder=100)
    offset = 0.04
    valid_indexes = ~((np.isnan(logx)) | (np.isnan(logy)))
    min_logx = np.nanmin(logx[valid_indexes])
    max_logx = np.nanmax(logx[valid_indexes])
    min_logy = np.nanmin(logy[valid_indexes])
    max_logy = np.nanmax(logy[valid_indexes])
    out_x=abs(max_logx-min_logx) * offset
    out_y=abs(max_logy-min_logy) * offset
    ax.set_xlim(min_logx-out_x,max_logx+out_x)
    ax.set_ylim(-0.75, max_logy+out_y)
    ax.scatter(logx, logy, alpha=0.1, color="gray", zorder=-10, ls="None")
    grid_ax.set_xlim(np.nanmin(logx)-out_x,np.nanmax(logx)+out_x)
    grid_ax.set_ylim(-0.75, np.nanmax(logy)+out_y)
    grid_ax.scatter(logx, logy, alpha=0.1, color="gray", zorder=-10, ls="None")

grid_axes[0].legend(handles=patches)
#ax_1.set_xlim(-1.5, 1)
#ax_1.set_ylim(-1., 1.5)
ax_1.legend(handles=patches)

if not os.path.isdir(outdir):
    os.mkdir(outdir)
outfile = "%s.png" % args.model
grid_figure.savefig("%s/grid_fig%s" % (outdir, outfile), bbox_inches='tight')
# ax_1.text(-0.5, 0.7, "AGN", fontsize=20, weight='bold', clip_on=True)
# ax_1.text(-0.2, 0.1, "LINER", fontsize=20, weight='bold', clip_on=True)
# ax_1.text(-0.5, 0.35, "Inter", fontsize=20, weight='bold', clip_on=True)
# ax_1.text(-1.1, 0.2, "HII", fontsize=20, weight='bold', clip_on=True)
# ax_3.text(-1.1, 0.7, "AGN", fontsize=20, weight='bold', clip_on=True)
# ax_3.text(-0.6, 0.1, "LINER", fontsize=20, weight='bold', clip_on=True)
# ax_3.text(-0.6, -0.6, "Inter", fontsize=20, weight='bold', clip_on=True)
# ax_3.text(-2, -0.6, "HII", fontsize=20, weight='bold', clip_on=True)
#ax_2.legend(handles=patches)
#ax_3.legend(handles=patches)


figure_bpt1.savefig("%s/bpt1%s" % (outdir, outfile))
figure_bpt2.savefig("%s/bpt2%s" % (outdir, outfile))
figure_bpt3.savefig("%s/bpt3%s" % (outdir, outfile))
grid_figure.savefig("%s/grid_fig%s" % (outdir, outfile), bbox_inches='tight')
print("Results stored to %s" % outdir)
