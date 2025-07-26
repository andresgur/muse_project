# -*- coding: utf-8 -*
# !/usr/bin/env python3




class Line:
    """This class is the basic for defining a line. A line is defined by its name, its wavelength and the reference line to which it is attached if the ratio has to be constrained.

    """
    def __init__(self, name, wave, ref=None, fit=False, low=0, up=None, th=None, index=None):
        self.index = index
        self.name = name
        self.wave = wave
        self.fit = fit
        self.ref = ref
        if ref is not None:
            self.low = low
            self.up = up
            self.th = th
    def __str__(self):
        return self.name


class Lines:
    """This class enables to deal with lines.  A dictionary stored in lines will contain informations on each lines.

    """

    def append(self, line):
        self.lines[line.name] = line
        self.lines[line.name].index = self.index
        self.index += 1

    def __init__(self):
        self.index = 0
        self.lines = {}
        # ref always to the higher wavelength line
        self.append(Line("HeII4686", 4685.682))
        self.append(Line('HBETA', 4861.32))
        self.append(Line('OIII4959', 4958.911))
        self.append(Line('OIII5007', 5006.843, ref='OIII4959', low=2.5, up=3.3, th=2.94))
        self.append(Line("HeI5875", 5875))
        self.append(Line("NII5755", 5754.8))
        self.append(Line('OI6300', 6300.3))
        self.append(Line("OI6363", 6363.78, ref="OI6300", low=0.303, up=0.4, th=0.3289)) # theoreitical is 3.04 (OI6300/OI6363 = 3.04) so we use the inverse ratio
        self.append(Line('HALPHA', 6563)) # ref='HBETA', low=2.65, th=2.85, up=6))
        self.append(Line('NII6548', 6548.05))
        self.append(Line('NII6583', 6583.45, ref='NII6548', low=2.5, up=3.3, th=3.))
        self.append(Line("HeI6678", 6678.15))
        # there ratios are the inverse of 0.4 and 1.5 respectively
        
        self.append(Line('SII6716', 6716.44))
        self.append(Line('SII6731', 6730.82, ref='SII6716', low=0.65, up=3, th=1.)) # this is from 1.53 to 0.33 (see Osterbroks & Ferland 2006, p. 123)
        self.append(Line("HeI7065", 7065.18))
        self.append(Line("ArIII7135", 7135.79))
        self.append(Line("SIII9069", 9069.27))