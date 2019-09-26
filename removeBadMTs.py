#!/usr/bin/python
#coding: utf-8

import os
import sys
import optparse
import numpy as np
from EMAN2 import *
from StarRW import *

def setupParserOptions():
        parser = optparse.OptionParser()
        parser.add_option("-f", dest="fpar", type="string", metavar="FILE", help="FREALIGN par file of particles")
        parser.add_option("-s", dest="stack", type="string", metavar="FILE", help="Particle stack file in mrcs format")
        parser.add_option("--star", dest="star", type="string", metavar="FILE", help="Relion star file of particles")
        options, args = parser.parse_args()
        if len(args) > 1: parser.error("Unknown command-line options: %s"%str(args))
        if len(sys.argv) < 4:
                parser.print_help()
                sys.exit()
        params = {}
        for i in parser.option_list:
                if isinstance(i.dest, str): params[i.dest] = getattr(options, i.dest)
        return params

def readPar(filename):
        l = [i.split() for i in open(filename, "r").readlines() if not i.startswith("C")]
        fpar = []
        for i in l:
                d = {"num":int(i[0]), "psi":float(i[1]), "theta":float(i[2]), "phi":float(i[3]), "shx":float(i[4]), "shy":float(i[5]), "mag":int(i[6]), "film":int(i[7]), "df1":float(i[8]), "df2":float(i[9]), "angast":float(i[10]), "occ":float(i[11]), "logp":float(i[12]), "sigma":float(i[13]), "score":float(i[14]), "change":float(i[15])}
                fpar.append(d)
        return fpar

def main(params):
        #file sthg
        par = readPar(params["fpar"])
        particles = EMData.read_images(params["stack"])
        fout2 = file("%s-BadMTremoved"%params["fpar"], "w")
        star = Star(params["star"])
        newcont = []
        #do sthg
        remove_MT = []
        philist = [i.split() for i in open("philist.txt", "r").readlines()]
        for l in philist:
                if float(l[2]) < 0: remove_MT.append(int(l[0]))
	im = EMData(params["stack"],0)
	nx = im.get_xsize()
        count = 0
        for i in range(len(par)):
                fpar = par[i]
                if fpar["film"] not in remove_MT:
                        newcont.append(star.content[i])
                        FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
                        fout2.write(FORMAT%(count+1, fpar["psi"], fpar["theta"], fpar["phi"], fpar["shx"], fpar["shy"], fpar["mag"], fpar["film"], fpar["df1"], fpar["df2"], fpar["angast"], fpar["occ"], fpar["logp"], fpar["sigma"], fpar["score"], fpar["change"]))
                        count += 1
        write_star(params["star"]+"_BadMTremoved", star.header, newcont)
        fout2.close()
        print "Allocating space for BadMTremoved-particles.mrc ...\n"
        MTavgstack = EMData(nx,nx,count)
        MTavgstack.write_image("BadMTremoved-particles.mrc")
        print "Done allocating"
        count = 0
        for i in range(len(par)):
                fpar = par[i]
                if fpar["film"] not in remove_MT:
                        region = Region(0,0,count,nx,nx,1)
                        particles[i].write_image("BadMTremoved-particles.mrc",0,EMUtil.get_image_ext_type("mrc"), False, region, EMUtil.EMDataType.EM_FLOAT, True)
                        count += 1
		if i % 100 == 0: print "working on ptcl %d\t\r"%i

if __name__ == "__main__":
        params = setupParserOptions()
        main(params)
        sys.exit()
