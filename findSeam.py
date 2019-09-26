#!/usr/bin/python
#coding: utf-8

import os
import sys
from EMAN2 import *
from scipy.ndimage import zoom
import optparse
import numpy as np
import datetime
import re
import time
from getpass import getpass

def setupParserOptions():
        parser = optparse.OptionParser()
        parser.add_option("-r", dest="ref", type="string", metavar="FILE", help="Reference volume in mrc format")
        parser.add_option("-f", dest="fpar", type="string", metavar="FILE", help="FREALIGN par file of super particles")
        parser.add_option("-s", dest="sp", type="string", metavar="FILE", help="Super particle stack file in mrcs format")
        parser.add_option("--pf", dest="pf", type="int", metavar="13/14", help="Number of protofilament")
        parser.add_option("--apix", dest="apix", type="float", metavar="FLOAT", help="Pixel size of micrograph")
        parser.add_option("--bin", dest="bin", type="int", metavar="INT", help="Binning factor")
        parser.add_option("--orad", dest="orad", type="float", metavar="FLOAT", help="Outer radius in Angstrom")
        parser.add_option("--servers", dest="servers", type="string", metavar="SERVER_NAMES", help="Names of servers you want to use (except for the mother server)")
        parser.add_option("--username", dest="username", type="string", metavar="USERNAME", help="Your name for log in to servers")
        parser.add_option("--continue", dest="contdir", type="string", metavar="DIRECTORY", help="Directory where you want to continue processing")
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

def numCPUs():
        cmd = 'cat /proc/cpuinfo |grep processor |wc'
        d = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        lines = d.stdout.readlines()
        lines = re.split('\s+', lines[0])
        number_of_procs = int(lines[1])
        return number_of_procs

def numCPUs_multi(server_name, cmd):
        cmd = cmd%(server_name, 'cat /proc/cpuinfo |grep processor |wc >> tmp.out')
        os.system(cmd)
        time.sleep(3)
        cont = open("tmp.out", "r").read().split()
        number_of_procs = int(cont[0])
        os.remove("tmp.out")
        if server_name.startswith("gpu"): return number_of_procs*3
        elif server_name.startswith("dhc"): return number_of_procs*3
        elif server_name == "dl380": return 48
        else: return number_of_procs

def main(params):
        if params["contdir"] is None:
                #prepare directory
                dt_now = str(datetime.datetime.now())
                rand = "_".join(dt_now.split())
                dirname = "tmp.%s"%rand
                os.system("mkdir %s"%dirname)
                os.chdir(dirname)
                os.system("ln -s ../%s ."%params["fpar"])
                os.system("ln -s ../%s ."%params["ref"])
                os.system("ln -s ../%s ."%params["sp"])
                os.system("cp ../pyfindSeam.py .")
                #read par file
                fpar = readPar(params["fpar"])
                nMTs = fpar[-1]["film"]
                imgs = EMData().read_images(params["sp"])
                #servers
                if params["servers"] is None: nServers = 1
                else:
                        servers = params["servers"].split(",")
                        nServers = len(servers)+1
                        NAME = params["username"]
                cmd_sshlogin = 'ssh %s@%s "cd %s;%s" &'%(NAME, "%s", os.getcwd(), "%s")
                #split into multi-processes
                CPUs = [numCPUs()]
                if nServers != 1: CPUs = CPUs + [numCPUs_multi(i, cmd_sshlogin) for i in servers]
                nCPUs = sum(CPUs)
                charge = nMTs/nCPUs+1 #the number of MTs each process will be in charge of
                MTrange = [0, 0]
                prev = 0
                print "nCPU:%d, charge:%d"%(nCPUs, charge)
                servers_willuse = []
                current_server = os.uname()[1]
                if nServers > 1: server_names = [current_server] + servers
                else: server_names = [current_server]
                for i in range(nServers):
                        servers_willuse += [server_names[i]]*CPUs[i]
                for i in range(nCPUs):
                        start = MTrange[1]+1
                        end = start+charge-1
                        if start > nMTs:
                                i -= 1
                                break
                        if end > nMTs: end = nMTs
                        MTrange = [start, end]
                        fpar_part = []
                        ind = prev
                        for pars in fpar[prev:]:
                                if MTrange[0] <= pars["film"] and pars["film"] <= MTrange[1]:
                                        fpar_part.append(pars)
                                else:
                                        prev = ind
                                        break
                                ind += 1
                        st, ed = fpar_part[0]["num"], fpar_part[-1]["num"]
                        writePar(fpar_part, i+1)
                        for j in range(st-1, ed):
                                imgs[j].write_image("MTSuperPtcl_Part%d.mrcs"%(i+1), -1)
                        com4debug = "python pyfindSeam.py -r %s -f %s -s %s --pf=%d --apix=%f --proc=%d --bin=%d --orad=%f --num=%d --debug &"%(params["ref"], "MTSuperPtcl_Part%d.par"%(i+1), "MTSuperPtcl_Part%d.mrcs"%(i+1), params["pf"], params["apix"], i+1, params["bin"], params["orad"], ed-st+1)
                        com4nodebug = "python pyfindSeam.py -r %s -f %s -s %s --pf=%d --apix=%f --proc=%d --bin=%d --orad=%f --num=%d &"%(params["ref"], "MTSuperPtcl_Part%d.par"%(i+1), "MTSuperPtcl_Part%d.mrcs"%(i+1), params["pf"], params["apix"], i+1, params["bin"], params["orad"], ed-st+1)
                        if params["debug"]:
                                if servers_willuse[i] == current_server:
                                        os.system(com4debug)
                                else:
                                        os.system(cmd_sshlogin%(servers_willuse[i], "conda activate pfp >> tmp.out; "+com4debug))
                        else:
                                if servers_willuse[i] == current_server:
                                        os.system(com4nodebug)
                                else:
                                        os.system(cmd_sshlogin%(servers_willuse[i], "conda activate pfp >> tmp.out; "+com4nodebug))
                print "%d subprocesses are running in total."%(i+1)
                del imgs
                print "Please wait for pyfindSeam.py to find the seam location of each MT ..."
        if params["contdir"] is not None: os.chdir(params["contdir"])
        while True:
                dirs = [i for i in os.listdir(os.getcwd()) if i.startswith("MTSuperPtcl_Part") and i.endswith("SeamSearch")]
                flag = True
                for dir in dirs:
                        if not os.path.isfile("%s/done"%dir): flag = False
                if flag: break
                time.sleep(60)
        list2 = [(re.search("[0-9]+", x).group(), x) for x in dirs]
        list2.sort(cmp = lambda x, y: cmp(int(x[0]), int(y[0])))
        dirs = [x[1] for x in list2]
        cont = []
        for dir in dirs:
                cont += open(dir+"/philist.txt", "r").readlines()
        open("philist.txt", "w").write("".join(cont))

def writePar(fpar_part, procid):
        fout = file("MTSuperPtcl_Part%d.par"%procid, "w")
        for i in range(len(fpar_part)):
                fpar = fpar_part[i]
                fout.write("%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n"%(fpar["num"], fpar["psi"], fpar["theta"], fpar["phi"], fpar["shx"], fpar["shy"], fpar["mag"], fpar["film"], fpar["df1"], fpar["df2"], fpar["angast"], fpar["occ"], fpar["logp"], fpar["sigma"], fpar["score"], fpar["change"]))
        fout.close()

if __name__ == "__main__":
        params = setupParserOptions()
        main(params)
        sys.exit()

