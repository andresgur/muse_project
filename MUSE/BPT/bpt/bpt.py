# @Author: Andrés Gúrpide <agurpide> and Maxime Parra
# @Date:   07-07-2023
# @Email:  a.gurpide-lasheras@soton.ac.uk
# !/usr/bin/env python3

import numpy as np
from abc import ABC, abstractmethod

class BPT:
    """Wrapper for the regions of the BPT diagram.
    Parameters
    ----------
    index: int
        Index to keep track of each region (0, 1, 2 or 3)
    name: str
        Name of the region (for plotting purposes mainly i.e. HII, LINER, etc)
    color: str
        Color to assign to this region be used in plots
    """
    def __init__(self):
        self.y_axis="log([OIII]/H$_\\beta$)"
        # initialize common attributes
        self.logx = None
        self.logy = None
        # limits for the lines to be drawn
        self.agnliner_inter = None # x-axis values
        self.int_inter = None # y-axis values
        # limit to avoid divergence in star forming line
        self.limit = None

    @abstractmethod
    def pure_starform_crv(self, logx):
        """
        Abstract method for the pure star-forming curve.
        Parameters:
            logx (ndarray): The logarithm of the x-values.
        Returns:
            ndarray: The y-values on the pure star-forming curve.
        """
        pass

    @abstractmethod
    def intermediate_crv(self, logy):
        """
        Abstract method for the intermediate curve
        Parameters:
            logy (ndarray): The logarithm of the y-values.
        Returns:
            ndarray: The x-values on the intermediate region curve.
        """
        pass
    @abstractmethod
    def agn_crv(self, logx):
        """
        Abstract method for the AGN curve
        Parameters:
            logx (ndarray): The logarithm of the x-values.
        Returns:
            ndarray: The y-values on the intermediate region curve.
        """
        pass

class BPT_1(BPT):
    def __init__(self):
        super().__init__()
        self.x_axis="log([NII]/H$_\\alpha$)"
        self.limit = -0.023
        self.int_inter = (-0.65, 1)
        # values from Section 4 Law+2021
        self.agnliner_inter = (-0.24, 0.5)
        self.x_axis="log([NII]/H$_\\alpha$)"

    def pure_starform_crv(self, logx):
        return 0.438 / (logx + 0.023) + 1.222
    def intermediate_crv(self, logy):
        return -0.390 * logy ** 4 - 0.582 * logy ** 3 - 0.637 * logy**2- 0.048*logy - 0.119
    def agnliner_crv(self, logx):
        return 0.95 * logx + 0.56

class BPT_2(BPT):
    def __init__(self):
        super().__init__()
        self.limit = 0.324
        self.int_inter = (-1.1, 1.10)
        # values from Section 4 Law+2021
        self.agnliner_inter = (-0.22, 0.3)
        self.x_axis="log([SII]/H$_\\alpha$)"
    def pure_starform_crv(self, logx):
        return 0.648 / (logx - 0.324) + 1.349
    def intermediate_crv(self, logy):
        return -1.107 * logy**4 - 0.489* logy**3+0.580* logy**2-0.579* logy - 0.043
    def agnliner_crv(self, logx):
        return 1.89 * logx + 0.76

class BPT_3(BPT):
    def __init__(self):
        super().__init__()
        self.limit = -0.124
        self.int_inter = (-0.25, 0.6)
        # values from Section 4 Law+2021
        self.agnliner_inter = (-0.9, 0.3)
        self.x_axis="log([OI]/H$_\\alpha$)"

    def pure_starform_crv(self, logx):
        return 0.884 / (logx + 0.124) + 1.291
    def intermediate_crv(sefl, logy):
        return 19.021 * logy**4 - 36.452* logy**3 + 21.741* logy**2 - 5.821* logy - 0.328
    def agnliner_crv(self, logx):
        return 1.18 * logx + 1.3
