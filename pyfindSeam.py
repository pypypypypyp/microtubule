#!/usr/bin/python

import os
import sys
from EMAN2 import *
from scipy.ndimage import zoom
import optparse
import numpy as np
import datetime

def setupParserOptions():
        parser = optparse.OptionParser()
        parser.add_option("-r", dest="ref", type="string", metavar="FILE", help="Reference volume in mrc format")
        parser.add_option("-f", dest="fpar", type="string", metavar="FILE", help="FREALIGN par file of super particles")
        parser.add_option("-s", dest="sp", type="string", metavar="FILE", help="Super particle stack file in mrcs format")
        parser.add_option("--pf", dest="pf", type="int", metavar="13/14", help="N of protofilament")
        parser.add_option("--apix", dest="apix", type="float", metavar="FLOAT", help="Pixel size of micrograph")
        parser.add_option("--proc", dest="proc", type="int", metavar="INT", help="Process number for multi-processing")
        parser.add_option("--bin", dest="bin", type="int", metavar="INT", help="Binning factor")
        parser.add_option("--num", dest="num", type="int", metavar="INT", help="Number of ptcls")
        parser.add_option("--orad", dest="orad", type="float", metavar="FLOAT", help="Outer radius in Angstrom")
        parser.add_option("--debug", dest="debug", action="store_true", default=False, help="Save debug files")
        options, args = parser.parse_args()
        if len(args) > 1: parser.error("Unknown command-line options: %s"%str(args))
        if len(sys.argv) < 2:
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

def cc(img1, img2):
        img_prod = np.fft.fft2(img1) * np.fft.fft2(img2).conj()
        cc_image = np.fft.fftshift(np.fft.ifft2(img_prod))
        return cc_image.real

def maskimg(img, edgewidth, vertical):
        h, w = img.shape
        img[:vertical, :] = 0
        img[h-vertical:, :] = 0
        img[:, :edgewidth] = 0
        img[:, w-edgewidth:] = 0
        return img

def numCPUs():
        cmd = 'cat /proc/cpuinfo |grep processor |wc'
        d = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        lines = d.stdout.readlines()
        lines = re.split('\s+', lines[0])
        number_of_procs = int(lines[1])
        return number_of_procs

def mymedian(mylist):
	import copy
	n = len(mylist)
        flag1 = False
        flag2 = False
        for i in mylist:
                if 0 <= i and i < 10: flag1 = True
                elif 350 < i and i <= 360: flag2 = True
        if flag1 and flag2:
                for i in range(n):
                        if 180 < mylist[i] and mylist[i] <= 360:
                                mylist[i] -= 360
	if n%2 == 1:
		return np.median(mylist)
	else:
		mylist2 = copy.copy(mylist)
		mylist2.sort()
		return np.median(mylist2[1:])

def main(params):
        #prepare directory
        dirname = "MTSuperPtcl_Part%d_SeamSearch"%params["proc"]
        os.system("mkdir %s"%dirname)
        os.chdir(dirname)
        os.system("mv ../%s ."%params["fpar"])
        os.system("ln -s ../%s ."%params["ref"])
        os.system("mv ../%s ."%params["sp"])
        fout = file("perptcl_stats.txt", "w")
        fout2 = file("philist.txt", "w")
        fout.write("#PTCL\t#FILM\tROT\tCC2\tCC3\tDIFFCC\n\n")
        #read reference 3d volume
        e = EMData(params["ref"])
        ZOOM = float(params["bin"])
        befclip = EMNumPy.em2numpy(e)
        #read par file
        fpar = readPar(params["fpar"])
        pf = params["pf"]
        apix = params["apix"]
        #constants
        SIZE = befclip.shape[0]
        ORAD = params["orad"]/apix/ZOOM
        GRID = 6
        MONOMER = 41.0
        TWIST = -360./pf
        MONOMER = MONOMER/apix/ZOOM
        RISE = MONOMER*3/pf
        DIMER = MONOMER*2
        #mainloop for each super particles
        clipw = int(SIZE-0.9*(SIZE-max(SIZE-2*(ORAD+MONOMER), ORAD*2)))/2
        clipw = clipw + 4-clipw%4
        clip = befclip[befclip.shape[0]/2-clipw:befclip.shape[0]/2+clipw, befclip.shape[1]/2-clipw:befclip.shape[1]/2+clipw, befclip.shape[2]/2-clipw:befclip.shape[2]/2+clipw]
        REF = EMNumPy.numpy2em(zoom(clip, 1./ZOOM, order=1))
        ccs_perfilm = []
        rots_perfilm = []
        prevfilm = fpar[0]["film"] #previous film number
        for i in range(params["num"]):
                ccs = []
                rots = []
                e = EMData(params["sp"], i)
                img = EMNumPy.em2numpy(e)
                img = zoom(img, 1/ZOOM, order=1)
                for j in range(pf):
                        #read i-th super particle
                        size = img.shape[0]
                        #read parameters
                        rot = fpar[i]["phi"]
                        theta = fpar[i]["theta"]
                        correspondingProjRot = (rot+(j-GRID)*TWIST)%360
                        rots.append(correspondingProjRot)
                        #transform reference volume
                        t = Transform()
                        t.set_rotation({"type":"spider", "phi":correspondingProjRot, "theta":theta})
                        REFrot = REF.process("xform", {"transform":t})
                        projbef = np.sum(EMNumPy.em2numpy(REFrot), axis=0)
                        befsize = projbef.shape[0]
                        proj = np.zeros((size, size))
                        proj[size/2-befsize/2:size/2+befsize/2, size/2-befsize/2:size/2+befsize/2] = projbef
                        shifty = int(round((rot-correspondingProjRot)/TWIST))*RISE
                        #set shifty in range(-MONOMER, MONOMER)
                        shifty = shifty%DIMER
                        if shifty > MONOMER: shifty -= DIMER
                        shifty = shifty*np.sin(theta/180.*np.pi)
                        #round shifty to int
                        sy =  int(round(shifty))
                        img2 = np.zeros(img.shape)
                        img3 = np.zeros(img.shape)
                        sy2 = int(round(shifty+MONOMER))
                        #shift images
                        if sy >= 0:
                                img2[:, sy:] = img[:, :img2.shape[1]-sy]
                        if sy < 0:
                                img2[: ,:sy] = img[:, -sy:]
                        if sy2 >= 0:
                                img3[:, sy2:] = img[:, :img2.shape[1]-sy2]
                        if sy2 < 0:
                                img3[: ,:sy2] = img[:, -sy2:]
                        #mask
                        maskwidth = int(ORAD+MONOMER)
                        img2 = maskimg(img2, maskwidth, size/2-int(ORAD))
                        img3 = maskimg(img3, maskwidth, size/2-int(ORAD))
                        clip = max(img2.shape)
                        img2 = img2[img2.shape[0]/2-clip/2:img2.shape[0]/2+clip/2, img2.shape[1]/2-clip/2:img2.shape[1]/2+clip/2]
                        img3 = img3[img3.shape[0]/2-clip/2:img3.shape[0]/2+clip/2, img3.shape[1]/2-clip/2:img3.shape[1]/2+clip/2]
                        if params["debug"]:
                                EMNumPy.numpy2em(img2).write_image("cc2.mrcs", -1)
                                EMNumPy.numpy2em(img3).write_image("cc3.mrcs", -1)
                                EMNumPy.numpy2em(proj).write_image("proj.mrcs", -1)
                        #calculate cc
                        xlim1, xlim2 = proj.shape[0]/2-int(4/ZOOM), proj.shape[1]/2+int(4/ZOOM)
                        ylim1, ylim2 = proj.shape[0]/2-int(4/ZOOM), proj.shape[1]/2+int(4/ZOOM)
                        if params["debug"]:
                                EMNumPy.numpy2em(cc(img2, proj)).write_image("cc_1.mrcs", -1)
                                EMNumPy.numpy2em(cc(img3, proj)).write_image("cc_2.mrcs", -1)
                        cc2 = cc(img2, proj)[xlim1:xlim2, ylim1:ylim2].max()
                        cc3 = cc(img3, proj)[xlim1:xlim2, ylim1:ylim2].max()
                        ccpart = cc2 - cc3
                        ccs.append(ccpart)
                        fout.write("%d\t%d\t%.1f\t%.1f\t%.1f\t%.1f\n"%(i+1, fpar[i]["film"], correspondingProjRot, cc2, cc3, ccpart))
                fout.write("\n")
                if prevfilm != fpar[i]["film"] or i == params["num"]-1:
                        #diffcc_abs = np.sum(np.abs(np.array(ccs_perfilm)), axis=0)
                        #diffcc = np.sum(np.array(ccs_perfilm), axis=0)
                        #diffcc_abs = list(diffcc_abs)
                        #index = diffcc_abs.index(max(diffcc_abs))
                        diffcc = np.sum(np.array(ccs_perfilm), axis=0)
                        absdiffcc = list(np.abs(diffcc))
                        index = absdiffcc.index(max(absdiffcc))
                        phi = mymedian(np.array(rots_perfilm)[:, index])
                        phi %= 360
                        fout2.write("%d\t%.2f\t%.1f\n"%(prevfilm, phi, diffcc[index]))
                        ccs_perfilm = []
                        rots_perfilm = []
                        ccs_perfilm.append(ccs)
                        rots_perfilm.append(rots)
                elif prevfilm == fpar[i]["film"]:
                        ccs_perfilm.append(ccs)
                        rots_perfilm.append(rots)
                prevfilm = fpar[i]["film"] #previous film number
                print "proc %d: %d/%d done."%(params["proc"], i+1, params["num"])
        fout.close()
        fout2.close()
        print "%s => Done."%dirname
        done = file("done", "w")
        done.write("DONE")
        done.close()

if __name__ == "__main__":
        params = setupParserOptions()
        main(params)
        sys.exit()

