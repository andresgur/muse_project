# @Author: Andrés Gúrpide <agurpide>
# @Date:   11-12-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 11-12-2021

from photutils.aperture import CircularAperture, CircularAnnulus, EllipticalAperture, SkyCircularAperture, SkyCircularAnnulus, SkyEllipticalAperture

def get_image_filter(header):
    if "FILTER" in header:
        return header["FILTER"]
    elif "FILTNAM1" in header:
        return header["FILTNAM1"]
    else:
        return header["FILTER1"]
    return filter



def region_to_aperture(region):
    """Convert region object to aperture object."""

    region_type = type(region).__name__
    if "Pixel" in region_type:
        source_center = (region.center.x, region.center.y)
        if region_type == 'CirclePixelRegion':
            return CircularAperture(source_center, r=region.radius)
        elif region_type == "CircleAnnulusPixelRegion":
            return CircularAnnulus(source_center, r_in=region.inner_radius, r_out=region.outer_radius)
        elif region_type == "EllipsePixelRegion":
            # to be tested
            return EllipticalAperture(source_center, a=region.width / 2, b=region.height / 2, theta=region.angle)
    elif "Sky" in region_type:
        center = region.center.fk5
        if region_type == "CircleSkyRegion":
            return SkyCircularAperture(center, r=region.radius)
        elif region_type == "EllipseSkyRegion":
            return SkyEllipticalAperture(center, a=region.width / 2, b=region.height / 2, theta=region.angle)
        elif region_type == "CircleAnnulusSkyRegion":
            return SkyCircularAnnulus(center, r_in=region.inner_radius, r_out=region.outer_radius)
    else:
        TypeError("Error region not implemented")
        return None
