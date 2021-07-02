# @Author: Andrés Gúrpide <agurpide>
# @Date:   07-06-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 30-06-2021

import argparse
import numpy as np
import matplotlib.pyplot as plt
import re
import muse_utils as mu
import matplotlib.patches as mpatches
import bpt_config as bptcfg
from astropy.io import fits


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

        if self.index == 1:
            self.pure_starform_max = -0.032
            self.agnliner_inter = (-0.24, 0.5)
            self.int_inter = (-0.61, 1)
        elif self.index == 2:
            self.pure_starform_max=0.198
            self.agnliner_inter = (-0.22, 0.3)
            self.int_inter = (-1.1, 1)
        elif self.index == 3:
            self.pure_starform_max=-0.36
            self.agnliner_inter = (-0.9, 0.3)
            self.int_inter = (-0.25, 0.65)


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


parser = argparse.ArgumentParser(description='Script to read and plot the output of the mappings III libraries (Allen et al. 2008)')
parser.add_argument("mappings", nargs='+', help="Input mappings file")
args = parser.parse_args()
plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')

bpt_type = bptcfg.type
lineratiomaps = [bptcfg.lineratio_paths["n2_ha"], bptcfg.lineratio_paths["s2_ha"], bptcfg.lineratio_paths["oI_ha"]]
grid_figure, grid_axes = plt.subplots(1, 3, sharey=True, gridspec_kw={"wspace": 0.1}, figsize=(12.8, 7.7))
figure_bpt1, ax_1 = plt.subplots(1)
figure_bpt2, ax_2 = plt.subplots(1)
figure_bpt3, ax_3 = plt.subplots(1)
figure_balmer_decrement, balmer_decrement_ax = plt.subplots(1)
patches = []
speeds = ["100", "150", "200", "250", "300"]
markers = mu.get_markers_array(len(speeds))
models = ["M", "R", "J", "P", "Q", "T", "U", "V", "L", "S"]
colors = mu.create_color_array(len(models), "Dark2")
labels = ["Solar abundance, n = 1 cm$^{-3}$", "Twice solar abundance, n = 1 cm$^{-3}$", "Dopita et al. 2005, n = 1 cm$^{-3}$", "SMC, n = 1 cm$^{-3}$", "LMC, n = 1 cm$^{-3}$",
          "Solar abundance, n = 0.01 cm$^{-3}$", "Solar abundance, n = 0.1 cm$^{-3}$", "Solar abundance, n = 10 cm$^{-3}$", "Solar abundance, n = 100 cm$^{-3}$", "Solar abundance, n = 1000 cm$^{-3}$"]
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
out_file = "_".join(models)
# Get all models starting with a given letter
print(args.mappings)
for model, color, label in zip(models, colors, labels):
    input_models = [i for i in args.mappings if i.startswith(model)]
    if len(input_models) > 0:
        patches.append(mpatches.Patch(color=color, ls="solid", label=label))
    else:
        continue

    print("Processing %d models %s: %s" % (len(input_models), model, label))
    linestyles = mu.create_linestyle_array(len(input_models))
    for mp, ls in zip(input_models, linestyles):
        print("Reading file %s" % mp)
        b_value = re.search("b.", mp)[0]
        n_density = re.search("n\d", mp)[0]
        nomenclature = mp[0]

        data = np.genfromtxt(mp, skip_header=10, dtype=[('ID', 'U8'), ('Atom', 'U8'), ('Species', 'U8'), ('Wavelength', 'f8'),
                    ('100', '<f8'), ('125', '<f8'), ('150', '<f8'), ('175', '<f8'), ('200', '<f8'), ('225', '<f8'), ('250', '<f8'), ('275', '<f8'),
                    ('300', '<f8'), ('325', '<f8'), ('350', '<f8'), ('375', '<f8'), ('400', '<f8'), ('425', '<f8'), ('450', '<f8'), ('475', '<f8'), ('500', '<f8'),
                    ('525', '<f8'), ('550', '<f8'), ('575', '<f8'), ('600', '<f8'), ('625', '<f8'), ('650', '<f8'), ('675', '<f8'), ('700', '<f8'), ('725', '<f8'),
                    ('750', '<f8'), ('775', '<f8'), ('800', '<f8'), ('825', '<f8'), ('850', '<f8'), ('875', '<f8'), ('900', '<f8'), ('925', '<f8'), ('950', '<f8'),
                    ('975', '<f8'), ('1000', '<f8')], usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32))

        halpha_flux = data[(data["Atom"] == "H") & (data["Species"] == "I") & (data["Wavelength"] == 6562.8000)]
        #hbeta_flux = np.where((data["Atom"] == "H") & (data["Species"] == "I") & (data["Wavelength"] == 4861.3200)) # hbeta is set to 1
        nII_flux = data[(data["Atom"] == "N") & (data["Species"] == "II") & (data["Wavelength"] == 6583.3400)]
        sII_I_flux = data[(data["Atom"] == "S") & (data["Species"] == "II") & (data["Wavelength"] == 6716.3100)]
        sII_II_flux = data[(data["Atom"] == "S") & (data["Species"] == "II") & (data["Wavelength"] == 6730.6800)]
        oI_flux = data[(data["Atom"] == "O") & (data["Species"] == "I") & (data["Wavelength"] == 6300.2000)]
        oIII_flux = data[(data["Atom"] == "O") & (data["Species"] == "III") & (data["Wavelength"] == 5006.7700)]
        ratios_y = []
        ratios_1 = []
        ratios_2 = []
        ratios_3 = []
        halphas = []
        for speed, marker in zip(speeds, markers):
            ratio1 = np.log10(nII_flux[speed] / halpha_flux[speed])
            ratio_y = np.log10(oIII_flux[speed])
            ax_1.scatter(ratio1, ratio_y, color=color, s=float(speed), marker=marker)
            grid_axes[0].scatter(ratio1, ratio_y, color=color, s=float(speed), marker=marker)
            ratios_1.append(ratio1)
            ratios_y.append(ratio_y)
            ratio2 = np.log10((sII_I_flux[speed] + sII_II_flux[speed]) / halpha_flux[speed])
            ax_2.scatter(ratio2, ratio_y, color=color, s=float(speed), marker=marker)
            grid_axes[1].scatter(ratio2, ratio_y, color=color, s=float(speed), marker=marker)
            ratios_2.append(ratio2)
            ratio3 = np.log10(oI_flux[speed] / halpha_flux[speed])
            ax_3.scatter(ratio3, ratio_y, color=color, s=float(speed), marker=marker)
            grid_axes[2].scatter(ratio3, ratio_y, color=color, s=float(speed), marker=marker)
            ratios_3.append(ratio3)
            halphas.append(halpha_flux[speed])
            balmer_decrement_ax.scatter(speed, halpha_flux[speed], color=color, s=float(speed), marker=marker)
        for ax, ratios in zip(grid_axes, [ratios_1, ratios_2, ratios_3]):
            ax.plot(ratios, ratios_y, color=color, ls=ls)
        ax_1.plot(ratios_1, ratios_y, color=color, ls=ls)
        ax_2.plot(ratios_2, ratios_y, color=color, ls=ls)
        ax_3.plot(ratios_3, ratios_y, color=color, ls=ls)
        balmer_decrement_ax.plot(speeds, halphas, color=color, ls=ls)
# plot BPT lines and data points
y_axis_map = fits.open(bptcfg.lineratio_paths["o3_hb"])
logy = np.log10(y_axis_map[0].data)
for bpt, ax, map_1, grid_ax in zip([1, 2, 3], [ax_1, ax_2, ax_3], lineratiomaps, grid_axes):
    x_axis_map = fits.open(map_1)
    logx = np.sort(np.log10(x_axis_map[0].data))
    bpt_diagram = BPT_diagram(bpt)
    inter_starform = np.sort(logx[logx < bpt_diagram.pure_starform_max])
    ax.plot(inter_starform, bpt_diagram.pure_starform_crv(inter_starform), color="black", zorder=100)
    grid_ax.plot(inter_starform, bpt_diagram.pure_starform_crv(inter_starform), color="black", zorder=100)
    inter_def = logy[logy > bpt_diagram.int_inter[0]]
    inter_def = np.sort(inter_def[inter_def < bpt_diagram.int_inter[1]])
    ax.plot(bpt_diagram.int_crv(inter_def), inter_def, color='black', zorder=100)
    grid_ax.plot(bpt_diagram.int_crv(inter_def), inter_def, color='black', zorder=100)
    agnliner_def = logx[logx > bpt_diagram.agnliner_inter[0]]
    agnliner_def = np.sort(agnliner_def[agnliner_def < bpt_diagram.agnliner_inter[1]])
    # zorder 100 so it stands above the points
    ax.plot(agnliner_def, bpt_diagram.agnliner_crv(agnliner_def), color='red', zorder=100)
    grid_ax.plot(agnliner_def, bpt_diagram.agnliner_crv(agnliner_def), color='red', zorder=100)
    out_x=abs(np.nanmax(logx)-np.nanmin(logx)) * 0.1
    out_y=abs(np.nanmax(logy)-np.nanmin(logy)) * 0.08
    ax.set_xlim(np.nanmin(logx)-out_x,np.nanmax(logx)+out_x)
    ax.set_ylim(np.nanmin(logy)-out_y,np.nanmax(logy)+out_y)
    ax.scatter(logx, logy, alpha=0.1, color="gray", zorder=-10, ls="None")
    grid_ax.set_xlim(np.nanmin(logx)-out_x,np.nanmax(logx)+out_x)
    grid_ax.set_ylim(np.nanmin(logy)-out_y,np.nanmax(logy)+out_y)
    grid_ax.scatter(logx, logy, alpha=0.1, color="gray", zorder=-10, ls="None")
grid_axes[0].legend(handles=patches, fontsize=18)
ax_1.legend(handles=patches)
#ax_2.legend(handles=patches)
#ax_3.legend(handles=patches)
balmer_decrement_ax.legend(handles=patches)
outfile = out_file.replace(".txt#", "") + ".png"
figure_bpt1.savefig("bpt1_%s" % outfile)
figure_bpt2.savefig("bpt2_%s" % outfile)
figure_bpt3.savefig("bpt3_%s" % outfile)
figure_balmer_decrement.savefig("balemer_decrement.png")
grid_figure.savefig("grid_fig%s.png" % outfile)
print("Results stored to %s" % outfile)
