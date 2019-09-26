#!/usr/bin/python
#coding: utf-8

import os
import sys
from StarRW import *

def readPar(filename):
        l = [i.split() for i in open(filename, "r").readlines() if not i.startswith("C")]
        fpar = []
        for i in l:
                d = {"num":int(i[0]), "psi":float(i[1]), "theta":float(i[2]), "phi":float(i[3]), "shx":float(i[4]), "shy":float(i[5]), "mag":int(i[6]), "film":int(i[7]), "df1":float(i[8]), "df2":float(i[9]), "angast":float(i[10]), "occ":float(i[11]), "logp":float(i[12]), "sigma":float(i[13]), "score":float(i[14]), "change":float(i[15])}
                fpar.append(d)
        return fpar

if len(sys.argv) != 4:
        print "python par2star.py input.par ref.star output.star"
        sys.exit()

par = readPar(sys.argv[1])
star = Star(sys.argv[2])

psiindex = star.find_parameter_index("_rlnAnglePsiPrior")
thetaindex = star.find_parameter_index("_rlnAngleTiltPrior")
xindex = star.find_parameter_index("_rlnOriginX")
yindex = star.find_parameter_index("_rlnOriginY")
apixindex = star.find_parameter_index("_rlnDetectorPixelSize")
magindex = star.find_parameter_index("_rlnMagnification")
star.add_parameter("_rlnAngleRotPrior")
phiindex = star.find_parameter_index("_rlnAngleRotPrior")

cont = []
APIX = star.content[0][apixindex]/star.content[0][magindex]*10000
print "APIX=%f"%APIX

for i in range(len(par)):
        star.content[i][psiindex] = par[i]["psi"]
        star.content[i][thetaindex] = par[i]["theta"]
        star.content[i][phiindex] = par[i]["phi"]
        star.content[i][xindex] = -par[i]["shx"]/APIX
        star.content[i][yindex] = -par[i]["shy"]/APIX

star.write(sys.argv[3])
