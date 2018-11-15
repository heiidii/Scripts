import os
import sys


execfile("Inputs.py")
cmd.load(native+".pdb")

def openbytraj():
 for idx in range(start,stop,step): 
   curname = basename %(trajId,idx)
   cmd.load(curname)
   argname = curname.split(".pdb")[0]
   print "/%s//A" %argname,"1ACB_b//E"
   cmd.align("/%s//A" %argname,"/1ACB_b//E")

def openfromfile():
   f = open(filename,'r')
   lines=f.readlines()
   f.close()
   count=0
   for line in lines[start:stop:step]:
      if line.find('#') != -1: continue
      if count>countmax:
         print count,line
         break
      curfile = line.split()[0]+".pdb"
      cmd.load(curfile)
      argname = curfile.split(".pdb")[0].split('/')[1]
      #print "/%s//A" %argname,"1_b//E"
      cmd.align("/%s" %argname,"/%s" %native)
      count+=1

if mode=="file": openfromfile()
#elif mode=="traj": openbytraj()
else: "Need a mode: Use file or traj"
