import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpdaf.obj import Image, Cube
from regions import Regions
import os
import astropy.units as u
from regions import CirclePixelRegion
from regions import PixCoord


def wrap_to_pi(angle):
    # Wrap the angle from 0 to 2*pi
    angle = np.mod(angle, 2*np.pi)
    # Convert the angle to the range from -pi to pi
    angle[angle > np.pi] = angle[angle > np.pi] - 2 * np.pi
    angle[angle < -np.pi] = angle[angle < -np.pi] + 2 * np.pi
    return angle


def create_plot(image):
    """Creates image with sectors"""
    image.unmask()
    image.mask = True * np.ones(image.data.shape)

    r, theta = image.wcs.coord(polar=True, spaxel=True)

    radial_expr = (r < max_r)

    for index in range(sectors):
        if start_angles[index] > np.pi/2 and end_angles[index] < -np.pi/2: # this is the special case where we cross the -pi +pi boundary
            sector_condition = radial_expr & ((theta < end_angles[index]) | (theta > start_angles[index]))
        else:
            sector_condition = radial_expr & (theta < end_angles[index]) & (theta > start_angles[index])

        image.mask[sector_condition] = False

    fig, ax = plt.subplots(subplot_kw={"projection": image.wcs.wcs})
    x, y  = 0.5 * np.array([image.wcs.naxis1 - 1, image.wcs.naxis2 - 1])
    # plot circles
    [CirclePixelRegion(PixCoord(x=x, y=y), radius=radius).plot(ax=ax, facecolor='none', edgecolor='red', lw=0.5) for radius in radii[1:]]
    plt.imshow(image.data)
    plt.xlabel("Ra")
    plt.ylabel("Dec")
    basename = os.path.basename(image.filename)
    outname = basename.replace(".fits", "_sectors_%d_max_%d.png" % (sectors, max_r))
    plt.savefig(outname)
    plt.close()
    print("Stored sector plot to %s" % outname)


ap = argparse.ArgumentParser(description='Extracts radial profiles from the input Images. All images are assumed to be the same size')
ap.add_argument("input_images", nargs="+", help="List of images to extract the profiles from")
ap.add_argument("--region", nargs=1, help='Region file to cut the image')
ap.add_argument("-r", "--radius", nargs="?", help='Maximum radius for the radial profiles. Default None', default=None, type=int)
ap.add_argument("-s", "--sectors", nargs="?", help='Number of radial directions or "sectors". Defaults to 4', default=4, type=int)
ap.add_argument("--start_angle", nargs="?", help='Start angle offset in deg. Default 0', default=0, type=float)
args = ap.parse_args()

sectors = args.sectors

starting_angle = args.start_angle / 180 * np.pi

width = 2 * np.pi / sectors / 1.2
sections = 2 * np.pi * np.arange(0, sectors) / sectors + starting_angle + np.pi #  + img.wcs.get_rot()
#sections = np.mod(sections + np.pi, 2 * np.pi) - np.pi
start_angles, end_angles = sections - width / 2, sections + width / 2
start_angles[start_angles < 0] = 2 * np.pi + start_angles[start_angles < 0]
end_angles[end_angles > 2*np.pi] = end_angles[end_angles > 2*np.pi] - 2 * np.pi
sections = wrap_to_pi(sections)
start_angles = wrap_to_pi(start_angles)
end_angles = wrap_to_pi(end_angles)

end_angles = np.mod(end_angles + np.pi, 2 * np.pi) - np.pi
width_deg = width / (2 * np.pi) * 360
print("Splitting image into %d sectors with width = %.2f deg" % (sectors, width_deg))

img = Image(args.input_images[0])
start = [img.wcs.naxis2 - 1, img.wcs.naxis1 - 1]
center = 0.5 * np.array([img.wcs.naxis2 - 1, img.wcs.naxis1 - 1])
max_img_r = np.hypot(np.abs(start[0] - center[0]), np.abs(start[1] - center[1]))
max_r = args.radius if args.radius is not None else max_img_r
radii = np.arange(0, max_r, 5)

if args.region is not None:
    reg = Regions.read(args.region, format="ds9")


for image_file in args.input_images:
    image = Image(image_file)
    basename = os.path.basename(image_file)
    if args.region is not None:
        image = image.subimage((reg[0].center.dec.to(u.deg).value, reg[0].center.ra.to(u.deg).value), size=reg[0].radius,
                       unit_size=reg[0].radius.unit)
    create_plot(image)
    r, theta = image.wcs.coord(polar=True, spaxel=True)
    # mask the entire image
    for index in range(sectors):
        means = []

        for i, radius in enumerate(radii[:-1]):
            radial_expr = (r > radii[i]) & (r < radii[i + 1])

            image.mask = True * np.ones(image.data.shape)

            if start_angles[index] > np.pi/2 and end_angles[index] < -np.pi/2: # this is the special case where we cross the -pi +pi boundary
                sector_condition = radial_expr & ((theta < end_angles[index]) | (theta > start_angles[index]))
            else:
                sector_condition = radial_expr & (theta < end_angles[index]) & (theta > start_angles[index])
            image.mask[sector_condition] = False #mask_selection(sector)
            image.mask[np.isnan(image.data)] = True
            avg_r = np.mean(image.data)
            means.append(avg_r)
        outputs = np.array([radii[:-1], means])
        np.savetxt(basename.replace(".fits", "%d.dat" % index), outputs.T, header="r\tmean", fmt="%.4f")

import glob
files = glob.glob(basename.replace(".fits", "[0-%d].dat" % (sectors - 1 )))

for file in files:
    data = np.genfromtxt(file, names=True)
    plt.plot(data["r"] * img.get_step(unit=u.arcsec)[0], data["mean"], label=file)

plt.legend()
plt.xlabel("Radial distance (arcsec)")
plt.ylabel("Flux")
plt.savefig(basename.replace(".fits", "_profile_%d.png" % sectors))
