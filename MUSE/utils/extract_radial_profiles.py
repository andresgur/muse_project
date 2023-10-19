import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpdaf.obj import Image, Cube
from regions import Regions
import os
import astropy.units as u
from regions import RectanglePixelRegion
from regions import PixCoord
import glob

def get_center(region):
    reg = Regions.read(region, format="ds9")
    return reg[0].center.dec.to(u.deg).value, reg[0].center.ra.to(u.deg).value


plt.style.use('~/.config/matplotlib/stylelib/paper.mplstyle')

def create_plot(image):
    """Creates image with nrectangles"""
    #image.unmask()
    fig, ax = plt.subplots(subplot_kw={"projection": image.wcs.wcs})

    for index in range(nrectangles):
        angle = np.degrees(central_angles[index])
        center_rec = center + (max_r / 2 * np.sin(central_angles[index]), max_r / 2 * np.cos(central_angles[index]))
        RectanglePixelRegion(PixCoord(x=center_rec[1], y=center_rec[0]), width, max_r, angle=(angle - 90) * u.deg).plot(ax=ax,
                                 facecolor="none", edgecolor="white", lw=2)
        plt.text(center_rec[1], center_rec[0], s=index, color="white", fontsize=22)
    image.plot(ax=ax, colorbar="v", scale='linear', show_xlabel=False, show_ylabel=False, zscale=True, cmap="cividis")
    im = ax.images
    im[-1].colorbar.ax.set_ylabel("Flux (%s)" % image.primary_header["BUNIT"], fontsize=18)
    im[-1].colorbar.ax.yaxis.set_label_position('right')
    ax.set_xlabel("Ra")
    ax.set_ylabel("Dec", labelpad=-1)
    basename = os.path.basename(image.filename)
    outname = basename.replace(".fits", "_n_%d_max_%d_w%d_o%.1f.pdf" % (nrectangles, max_r, width, offsetdeg))
    plt.savefig(outname)
    plt.close()
    print("Stored sector plot to %s" % outname)


ap = argparse.ArgumentParser(description='Extracts radial profiles from the input Images. All images are assumed to be the same size')
ap.add_argument("input_images", nargs="+", help="List of images to extract the profiles from")
ap.add_argument("--region", nargs=1, help='Region file for the center if not use the center of the image. Only circles in fk5 (ds9) allowed.', default=None)
ap.add_argument("-r", "--radius", nargs="?", help='Maximum radius for the rectangles in pixels. Default uses the whole image', default=None, type=int)
ap.add_argument("-n", "--nrect", nargs="?", help='Number of rectangles. Defaults to 4', default=4, type=int)
ap.add_argument("-s", "--step", nargs="?", help="Radial step in pixels. Defaults to 5", default=5, type=int)
ap.add_argument("-w", "--width", nargs="?", help="Rectangle width over which to average in pixels. Default 5", default=5, type=int)
ap.add_argument("--offset", nargs="?", help='Offset angle in deg. Default 0', default=0, type=float)
args = ap.parse_args()

nrectangles = args.nrect
offsetdeg = args.offset
offset = offsetdeg / 180 * np.pi

width = args.width
central_angles = 2 * np.pi * np.arange(0, nrectangles) / nrectangles + offset
print("Splitting image into %d rectangles with width = %.d pixels" % (nrectangles, width))

step = args.step
half_step = step / 2
half_width = args.width / 2


if args.region is not None:
    dec_ra = get_center(args.region[0])
else:
    dec_ra = None
for image_file in args.input_images:
    image = Image(image_file, ext=1)
    basename = os.path.basename(image_file)
    if dec_ra is None:
        center = 0.5 * np.array([image.wcs.naxis2 - 1, image.wcs.naxis1 - 1])
    else:
        center = image.wcs.sky2pix(dec_ra, unit=u.deg)[0] # this returns an array of centers or cooridnates, we just use one boviously
    start = [image.wcs.naxis2 - 1, image.wcs.naxis1 - 1]
    max_img_r = np.abs(start[0] - center[0])
    max_r = args.radius if args.radius is not None else max_img_r
    radii = np.arange(0, max_r, step)

    create_plot(image)

    for index in range(nrectangles):
        means = []
        stds = []
        angle = np.degrees(central_angles[index])

        for i, radius in enumerate(radii[:-1]):
            center_rec = center + (half_step * np.sin(central_angles[index]), half_step* np.cos(central_angles[index])) + ((i) * step * np.sin(central_angles[index]), i* step * np.cos(central_angles[index]))
            image.mask_region(center_rec, (half_width, half_step), unit_center=None, unit_radius=None,
                              posangle=angle - 90, inside=False)

            avg_r = np.mean(image.data)
            std_r = np.std(image.data)
            means.append(avg_r)
            stds.append(std_r)
            image.unmask() # unmask for next radius
        # store central pixel
        outputs = np.array([radii[:-1] + half_step, means, stds])
        np.savetxt(basename.replace(".fits", "_%d.dat" % index), outputs.T, header="r_px\tmean\tstd", fmt="%.4f")

    files = glob.glob(basename.replace(".fits", "") + "_[0-%d].dat" % (nrectangles - 1 ))

    plt.figure()
    for file in files:
        data = np.genfromtxt(file, names=True)
        plt.errorbar(data["r_px"] * image.get_step(unit=u.arcsec)[0], data["mean"],yerr=data["std"], fmt="-o", label=file)

    plt.legend()
    plt.xlabel("Radial distance (arcsec)")
    plt.ylabel("Flux (%s)" % image.primary_header["BUNIT"])
    plt.savefig(basename.replace(".fits", "_profile_%d.png" % nrectangles))
    plt.close()
