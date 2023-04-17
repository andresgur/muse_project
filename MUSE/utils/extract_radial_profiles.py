import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpdaf.obj import Image, Cube
from regions import Regions
import os
import astropy.units as u
from regions import RectanglePixelRegion
from regions import PixCoord


plt.style.use('~/.config/matplotlib/stylelib/paper.mplstyle')

def create_plot(image):
    """Creates image with nrectangles"""
    #image.unmask()
    center = 0.5 * np.array([image.wcs.naxis2 - 1, image.wcs.naxis1 - 1])
    fig, ax = plt.subplots(subplot_kw={"projection": image.wcs.wcs})

    for index in range(nrectangles):
        angle = np.degrees(central_angles[index])
        center_rec = center + (max_r / 2 * np.sin(central_angles[index]), max_r / 2 * np.cos(central_angles[index]))
        RectanglePixelRegion(PixCoord(x=center_rec[1], y=center_rec[0]), width, max_r, angle=(angle - 90) * u.deg).plot(ax=ax,
                                 facecolor="none", edgecolor="white", lw=2)
        plt.text(center_rec[1], center_rec[0], s=index, color="white", fontsize=22)
    image.plot(ax=ax, colorbar="v", scale='linear', show_xlabel=False, show_ylabel=False, zscale=True, cmap="cividis")
    im = ax.images
    im[-1].colorbar.ax.set_ylabel("[O III]$\lambda$5007/H$\\beta$", fontsize=18)
    im[-1].colorbar.ax.yaxis.set_label_position('right')
    ax.set_xlabel("Ra")
    ax.set_ylabel("Dec", labelpad=-1)
    basename = os.path.basename(image.filename)
    outname = basename.replace(".fits", "_n_%d_max_%d_w%d.pdf" % (nrectangles, max_r, width))
    plt.savefig(outname)
    plt.close()
    print("Stored sector plot to %s" % outname)


ap = argparse.ArgumentParser(description='Extracts radial profiles from the input Images. All images are assumed to be the same size')
ap.add_argument("input_images", nargs="+", help="List of images to extract the profiles from")
ap.add_argument("--region", nargs=1, help='Region file to cut the image', default=None)
ap.add_argument("-r", "--radius", nargs="?", help='Maximum radius for the rectangles. Default uses the whole image', default=None, type=int)
ap.add_argument("-n", "--nrect", nargs="?", help='Number of rectangles. Defaults to 4', default=4, type=int)
ap.add_argument("-s", "--step", nargs="?", help="Radial step in pixels. Defaults to 5", default=5, type=int)
ap.add_argument("-w", "--width", nargs="?", help="Rectangle width over which to average. Default 5", default=5, type=int)
ap.add_argument("--offset", nargs="?", help='Offset angle in deg. Default 0', default=0, type=float)
args = ap.parse_args()

nrectangles = args.nrect

offset = args.offset / 180 * np.pi

width = args.width
central_angles = 2 * np.pi * np.arange(0, nrectangles) / nrectangles + offset
print("Splitting image into %d rectangles with width = %.d pixels" % (nrectangles, width))

img = Image(args.input_images[0])
start = [img.wcs.naxis2 - 1, img.wcs.naxis1 - 1]
center = 0.5 * np.array([img.wcs.naxis2 - 1, img.wcs.naxis1 - 1])
max_img_r = np.abs(start[0] - center[0])
max_r = args.radius if args.radius is not None else max_img_r
step = args.step
half_step = step / 2
half_width = args.width / 2
radii = np.arange(0, max_r, step)

if args.region is not None:
    reg = Regions.read(args.region[0], format="ds9")

for image_file in args.input_images:
    image = Image(image_file)
    basename = os.path.basename(image_file)
    if args.region is not None:
        image = image.subimage((reg[0].center.dec.to(u.deg).value,
                               reg[0].center.ra.to(u.deg).value),
                               size=reg[0].radius.to(u.arcsec).value, unit_size=reg[0].radius.to(u.arcsec).unit)
    create_plot(image)

    # mask the entire image
    for index in range(nrectangles):
        means = []
        center_rec = center
        angle = np.degrees(central_angles[index])

        for i, radius in enumerate(radii[:-1]):
            center_rec = center + (half_step * np.sin(central_angles[index]), half_step* np.cos(central_angles[index])) + ((i) * step * np.sin(central_angles[index]), i* step * np.cos(central_angles[index]))
            image.mask_region(center_rec, (half_width, half_step), unit_center=None, unit_radius=None,
                              posangle=angle - 90, inside=False)

            avg_r = np.mean(image.data)
            means.append(avg_r)
            image.unmask() # unmask for next radius

        outputs = np.array([radii[:-1], means])
        np.savetxt(basename.replace(".fits", "_%d.dat" % index), outputs.T, header="r\tmean", fmt="%.4f")

import glob
files = glob.glob(basename.replace(".fits", "") + "_[0-%d].dat" % (nrectangles - 1 ))

for file in files:
    data = np.genfromtxt(file, names=True)
    plt.plot(data["r"] * img.get_step(unit=u.arcsec)[0], data["mean"],"-o", label=file)

plt.legend()
plt.xlabel("Radial distance (arcsec)")
plt.ylabel("Flux")
plt.savefig(basename.replace(".fits", "_profile_%d.png" % nrectangles))
