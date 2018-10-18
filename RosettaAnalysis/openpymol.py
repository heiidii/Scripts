import os
import sys


execfile("Inputs.py")
cmd.load("1ACB_b.pdb")

def openbytraj():
 for idx in range(start,stop,step): 
   curname = basename %(trajId,idx)
   cmd.load(curname)
   argname = curname.split(".pdb")[0]
   print "/%s//A" %argname,"1ACB_b//E"
   cmd.align("/%s//A" %argname,"/1ACB_b//E")

def openfromfile():
   tempfilename = "temp_%d" %trajId
   command = "awk -v ref_rep=\"%d\" '$2==ref_rep' %s >%s" %(trajId,filename,tempfilename)
   print command
   os.system(command)
   os.system("cat %s" %tempfilename)
   f = open(tempfilename,'r')
   lines=f.readlines()
   f.close()
   count=0
   for line in lines[start:stop:step]:
      if count>countmax:
         print count,line
         break
      curfile = line.split()[0]+".pdb"
      cmd.load(curfile)
      argname = curfile.split(".pdb")[0]
      print "/%s//A" %argname,"1ACB_b//E"
      cmd.align("/%s//A" %argname,"/1ACB_b//E")
      count+=1

if mode=="file": openfromfile()
#elif mode=="traj": openbytraj()
else: "Need a mode: Use file or traj"
