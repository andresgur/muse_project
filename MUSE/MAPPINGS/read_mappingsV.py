# @Author: Andrés Gúrpide <agurpide>
# @Date:   07-06-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 01-12-2021

import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
import muse_utils as mu
import matplotlib.patches as mpatches
import bpt_config as bptcfg
from astropy.io import fits
import glob
# assuming a gas to dust ratio of 0.1
abundances_logOH = {"0d0001": 6.61, "0d0002": 6.91, "0d0005": 7.30, "0d001": 7.61, "0d002": 7.91, "0d004": 8.21, "0d010": 8.61, "0d01524": 8.83}


class BPT_diagram:
    """Wrapper for the BPT diagram.
    Parameters
    ----------
    index: int
        Index to keep track of each region (0, 1, 2 or 3)
    name: str
        Name of the region (for plotting purposes mainly i.e. HII, LINER, etc)
    color: str
        Color to assign to this region be used in plots
    """
    def __init__(self, index):
        self.index = index

        if self.index in [1,2,3]:
            self.y_axis="log([OIII]/H$_\\beta$)"
            if self.index == 1:
                self.x_axis="log([NII]/H$_\\alpha$)"
                self.agnliner_inter = (-0.24, 0.5)
                self.int_inter = (-0.61, 1)
            elif self.index == 2:
                self.x_axis="log([SII]/H$_\\alpha$)"
                self.agnliner_inter = (-0.22, 0.3)
                self.int_inter = (-1.1, 1)
            elif self.index == 3:
                self.x_axis="log([OI]/H$_\\alpha$)"
                self.agnliner_inter = (-0.9, 0.3)
                self.int_inter = (-0.25, 0.65)

        if self.index=='proj':
            self.x_axis="P1 (0.77*N2+0.54*S2+0.33*R3)"
            self.y_axis="P2 (-0.57*N2+0.82*S2)"

    def pure_starform_crv(self, logx):
        if self.index == 1:
            return 0.359 / (logx + 0.032) + 1.083
        elif self.index == 2:
            return 0.41 / (logx - 0.198) + 1.164
        elif self.index == 3:
            return 0.612 / (logx + 0.360) + 1.179

    def int_crv(self, logy):
        if self.index == 1:
            return -0.479 * logy ** 4 - 0.594 * logy ** 3 - 0.542 * logy**2-0.056*logy-0.143
        elif self.index == 2:
            return -0.943*logy**4-0.45*logy**3+0.408*logy**2-0.61*logy-0.025
        elif self.index == 3:
            return 18.664*logy**4-36.343*logy**3+22.238*logy**2-6.134*logy-0.283

    def agnliner_crv(self, logx):
        if self.index == 1:
            return 0.95 * logx + 0.56
        elif self.index == 2:
            return 1.89 * logx + 0.76
        elif self.index == 3:
            return 1.18 * logx + 1.3

    def proj_x(self,logN2,logS2,logR3):
        return 0.77 * logN2 + 0.54 * logS2 + 0.33* logR3

    def proj_y(self, logN2, logS2):
        return -0.57 * logN2 + 0.82 * logS2

    def cold_projcrv(self, logy):
        return 2.597 * logy**3 - 1.865 * logy**2 + 0.105 * logy - 0.435

    def warm_projcrv(self, logy):
        return 3.4 * logy**3 - 2.233 * logy**2 - 0.184 * logy - 0.172


parser = argparse.ArgumentParser(description='Script to read and plot the output of the mappings III libraries (Allen et al. 2008)')
parser.add_argument("-a", "--abundances", nargs='+', help="Abundances", default=["0d01524"], type=str)
parser.add_argument("-d", "--densities", nargs='+', help="Densities", default=[1], type=float)
parser.add_argument("-m", "--model", nargs='?', help="Model", default="Gutkin16", choices=["Gutkin16", "Allen2008"])
parser.add_argument("-t", "--type", nargs='?', help="Shock, precursor or shock + precursors (s, p, s+p)", default="s", choices=["s", "p", "s+p"])
args = parser.parse_args()

style_file = '/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle'

if.os.isfile(style_file):
    plt.style.use(style_file)

bpt_type = bptcfg.type
lineratiomaps = [bptcfg.lineratio_paths["n2_ha"], bptcfg.lineratio_paths["s2_ha"], bptcfg.lineratio_paths["oI_ha"]]
grid_figure, grid_axes = plt.subplots(1, 3, sharey=True, gridspec_kw={"wspace": 0.05}, figsize=(25.6, 9.8))
figure_bpt1, ax_1 = plt.subplots(1)
figure_bpt2, ax_2 = plt.subplots(1)
figure_bpt3, ax_3 = plt.subplots(1)
patches = []
#speeds = [100.0, 150.0, 200.0, 250.0, 300.0]
speeds = [100.0, 200.0, 300.0]
#speeds = [200.0, 225, 250.0, 300.0, 350, 400, 450, 500, 550, 600, 700, 800, 850, 900, 1000]
markers = mu.get_markers_array(len(speeds))
if args.type == "s":
    folder = "shock"
elif args.type == "p":
    folder = "precursor"
elif args.type == "s+p":
    folder = "shock_plus_precursor"
print("Using %s models" % folder)
models = []
abund_labels = []

c_to_dust = "1d00"

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
outdir = "abund"
outdir += "_".join(args.abundances)
outdir += "_dens"
for dens in args.densities:
    outdir += "_%.1f" % dens
outdir += "_%s_%s_C%s" % (args.model, args.type, c_to_dust)

linestyles = mu.create_linestyle_array(len(models))

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
                #grid_ax.text(ratio[0], ratio_y[0], s="%d" % speed, fontsize=18, clip_on=True)
                grid_ax.text(ratio[0], ratio_y[0], s="%d" % speed, fontsize=22, clip_on=True)

        for ax, ratios in zip(grid_axes, [ratios_1, ratios_2, ratios_3]):
            ax.plot(ratios, ratios_y, color=color, ls=ls, lw=1)
        ax_1.plot(ratios_1, ratios_y, color=color, ls=ls, lw=1)
        ax_2.plot(ratios_2, ratios_y, color=color, ls=ls, lw=1)
        ax_3.plot(ratios_3, ratios_y, color=color, ls=ls, lw=1)
# plot BPT lines and data points
y_axis_map = fits.open(bptcfg.lineratio_paths["o3_hb"])
print("WARNING: values below -1.0 for  [OIII]/Hb are beign ignored")
y_axis_map[0].data[np.where(np.log10(y_axis_map[0].data) < -1.0)] = np.nan
logy = np.log10(y_axis_map[0].data.data)
for bpt, ax, map_1, grid_ax in zip([1, 2, 3], [ax_1, ax_2, ax_3], lineratiomaps, grid_axes):
    x_axis_map = fits.open(map_1)
    logx = np.log10(x_axis_map[0].data)
    bpt_diagram = BPT_diagram(bpt)

    invalid_pixels = np.where((np.isnan(logx)) | (np.isnan(logy)))
    inter_starform = np.sort(np.reshape(logx, np.size(logx)))
    if grid_ax is not None:
        if bpt == 3:
            range = inter_starform[np.where(inter_starform < -0.360)]
        else:
            range = inter_starform
    ax.plot(range, bpt_diagram.pure_starform_crv(range), color="black", zorder=100)

    grid_ax.plot(range, bpt_diagram.pure_starform_crv(range), color="black", zorder=100)
    inter_def = logy[logy > bpt_diagram.int_inter[0]]
    inter_def = np.sort(inter_def[inter_def < bpt_diagram.int_inter[1]])
    ax.plot(bpt_diagram.int_crv(inter_def), inter_def, color='black', zorder=100)
    grid_ax.plot(bpt_diagram.int_crv(inter_def), inter_def, color='black', zorder=100)
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
ax_1.text(-0.5, 0.7, "AGN", fontsize=20, weight='bold', clip_on=True)
ax_1.text(-0.2, 0.1, "LINER", fontsize=20, weight='bold', clip_on=True)
ax_1.text(-0.5, 0.35, "Inter", fontsize=20, weight='bold', clip_on=True)
ax_1.text(-1.1, 0.2, "HII", fontsize=20, weight='bold', clip_on=True)
ax_3.text(-1.1, 0.7, "AGN", fontsize=20, weight='bold', clip_on=True)
ax_3.text(-0.6, 0.1, "LINER", fontsize=20, weight='bold', clip_on=True)
ax_3.text(-0.6, -0.6, "Inter", fontsize=20, weight='bold', clip_on=True)
ax_3.text(-2, -0.6, "HII", fontsize=20, weight='bold', clip_on=True)
#ax_2.legend(handles=patches)
#ax_3.legend(handles=patches)


figure_bpt1.savefig("%s/bpt1%s" % (outdir, outfile))
figure_bpt2.savefig("%s/bpt2%s" % (outdir, outfile))
figure_bpt3.savefig("%s/bpt3%s" % (outdir, outfile))
grid_figure.savefig("%s/grid_fig%s" % (outdir, outfile.replace(".png", ".pdf")), bbox_inches='tight')
print("Results stored to %s" % outdir)
