# -*- coding: utf-8 -*
# !/usr/bin/env python3




class Line:
    """This class is the basic for defining a line. A line is defined by its name, its wavelength and the reference line to which it is attached if the ratio has to be constrained.

    """
    def __init__(self, name, wave, ref=None, fit=False, low=0, up=None, th=None, index=None, label=None):
        self.index = index
        self.name = name
        self.wave = wave
        self.fit = fit
        self.ref = ref
        if ref is not None:
            self.low = low
            self.up = up
            self.th = th
        self.label = label
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
        self.append(Line("HeII4686", 4685.682, label=rf"HeII$\lambda$4686"))
        self.append(Line("FeIII4701", 4701.55, label=rf"[Fe III]$\lambda$4701"))
        self.append(Line("ArIV4711", 4711.26, label=rf"[Ar IV]$\lambda$4711"))
        self.append(Line("ArIV4740", 4740.12, label=rf"[Ar IV]$\lambda$4740", ref="ArIV4711", low=0.25, up=30, th=1))
        self.append(Line('HBETA', 4861.32, label=rf'H$\beta$'))
        self.append(Line('OIII4959', 4958.911, label=fr'[O III]$\lambda$4959')) # theoreitical is 3.04 (OIII5007/OIII4959 = 3.04) so we use the inverse ratio
        self.append(Line('OIII5007', 5006.843, ref='OIII4959', low=2.5, up=3.5, th=2.94, label=r'[O III]$\lambda$5007'))
        self.append(Line("NII5755", 5754.8, label=rf"[N II]$\lambda$5755"))
        self.append(Line("HeI5875", 5875, label=rf"[He I]$\lambda$5875"))
        self.append(Line('OI6300_sky', 6300.3, label=rf"[O I]$\lambda$6300"))
        self.append(Line('OI6300', 6300.3, label=rf"[O I]$\lambda$6300"))
        self.append(Line("SIII6312", 6312.06, label=rf"[S III]$\lambda$6312"))
        self.append(Line("OI6363_sky", 6363.78, ref="OI6300_sky", low=0.103, up=0.7, th=0.3289)) # theoreitical is 3.04 (OI6300/OI6363 = 3.04) so we use the inverse ratio
        self.append(Line("OI6363", 6363.78, ref="OI6300", low=0.25, up=0.5, th=0.3289, label=rf"[O I]$\lambda$6363")) # theoreitical is 3.04 (OI6300/OI6363 = 3.04) so we use the inverse ratio
        self.append(Line('NII6548', 6548.05, label=rf"[N II]$\lambda$6548"))
        self.append(Line('HALPHA', 6563, label=rf"H$\alpha$")) # ref='HBETA', low=2.65, th=2.85, up=6))
        self.append(Line('NII6583', 6583.45, ref='NII6548', low=2.4, up=3.6, th=3., label=rf"[N II]$\lambda$6583"))
        self.append(Line("HeI6678", 6678.15, label=rf"He I$\lambda$6678"))
        # there ratios are the inverse of 0.4 and 1.5 respectively
        
        self.append(Line('SII6716', 6716.44, label=rf"[S II]$\lambda$6716"))
        self.append(Line('SII6731', 6730.82, ref='SII6716', low=0.65, up=3, th=1.)) # this is from 1.53 to 0.33 (see Osterbroks & Ferland 2006, p. 123)
        self.append(Line("HeI7065", 7065.18, label=rf"He I$\lambda$7065"))
        self.append(Line("ArIII7135", 7135.79, label=rf"[Ar III]$\lambda$7135"))
        self.append(Line("SIII9069", 9069.27, label=rf"[S III]$\lambda$9069"))