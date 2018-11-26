import os
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import glob
import mpl_toolkits


def plotgraph(xvals,yvals,pltinfo,bHold=False,bshow=False):
	print pltinfo
	if not bHold:
        	plt.title(pltinfo["title"]+"\n"+pltinfo["key"])
        	plt.xlabel(pltinfo["xlabel"])
        	plt.ylabel(pltinfo["ylabel"])
	if "ylow" in pltinfo:
        	plt.ylim(ymin=pltinfo["ylow"])
	if "yhigh" in pltinfo:
		plt.ylim(ymax=pltinfo["yhigh"])
	if "xlow" in pltinfo:
		plt.xlim(xmin=pltinfo["xlow"])
	if "xhigh" in pltinfo:
		plt.xlim(xmax=pltinfo["xhigh"])
	if "key" in pltinfo:
        	plt.scatter(xvals,yvals,alpha=pltinfo["alpha"],edgecolor=pltinfo["edgecolor"],s = pltinfo["scale"],color=pltinfo["color"],marker=pltinfo["marker"],label=pltinfo["key"])#,ecolor='xkcd:lavendar',capthick=4,capsize=2)
		plt.legend()
	else:
		plt.scatter(xvals,yvals,color=pltinfo["color"],marker=pltinfo["marker"])
        if (not bshow) and (not bHold): 
		#do not show and do not hold
		plt.savefig("scatter"+pltinfo["outfile"],dpi=600)
        if not bHold and (not bshow):
		plt.close()
                plt.clf()
	elif not bHold and ( bshow):
		plt.show()
	else:
		plt.hold(True)
	
	

def plotgraphs():
	#execfile(infofile)
	print pltinfo
        bLastFile=False
	ifile = 0
	for curfile in pltinfo["files"]:
		print curfile
		mapheader = dict()
		mapheader = getmapping(curfile)
		xindex=-1
		yindex=-1
		if pltinfo["xfield"] in mapheader:
			xindex = mapheader[pltinfo["xfield"]]
			print "found xfield"
		else:
			print "xfield not found"
			print mapheader
			exit()

		if pltinfo["yfield"] in mapheader:
			yindex = mapheader[pltinfo["yfield"]]
			print "found xfield",pltinfo["xfield"]
                else:
                        print "yfield not found",pltinfo["yfield"]
                        print mapheader
			exit()
		f = open(curfile,'r')
		lines = f.readlines()
		f.close()
		xvals=[]
		yvals=[]
		for line in lines[i0+1:]:
                	data = line.split()
                	if line.find("total_score") != -1: continue
                	addstring = []
			xvals.append(float(data[xindex]))
               		yvals.append(float(data[yindex]))
		print ifile
		temp = pltinfo["outfile"]
		temptitle = pltinfo["outfile"]
		if not pltinfo["hold"]:
			#Generating different plots for each file	
			pltinfo["outfile"] = temp + "_"+str(ifile)
			#pltinfo["title"] = temptitle + "\n"+ curfile
		else:
			pltinfo["color"] = pltinfo["colors"][ifile]
		pltinfo["key"] = pltinfo["keys"][ifile]
		if ifile==(len(pltinfo["files"])-1) and pltinfo["hold"]: pltinfo["hold"]=False
        	plotgraph(xvals,yvals,pltinfo,pltinfo["hold"],pltinfo["show"])
		pltinfo["outfile"]=temp
		pltinfo["title"]=temptitle
		ifile+=1

def printmapping(sfile):
        outf = open("MapHeader.txt",'w')
        mapheader = dict()
	mapheader = getmapping(sfile)
        i=0
        for entry in mapheader:
                outf.write("#%s\t%d\n" %(entry, mapheader[entry]))
	outf.close()
	
def getmapping(sfile):
	print "mapping %s ..." %sfile
        f = open(sfile,'r')
        lines = f.readlines()
        f.close()
        print "Lines in file ",len(lines)
        line0=lines[i0]
        headersplit = line0.split()
        print headersplit
        mapheader = dict()
        i=0
        for entry in headersplit:
                mapheader[entry] = i
                i+=1
                print entry,i
	return mapheader


def filterpdb(sfile,basename):
	print "filtering %s ..." %sfile
	f = open(sfile,'r')
	lines = f.readlines()
	f.close()
	outf = open("FilteredFile_%s.txt" %basename,'w')
	outfdata = open("FilteredFileData_%s.txt" %basename,'w')
	print "Lines in file ",len(lines)
	line0=lines[i0]
	outfdata.write(line0)
	headersplit = line0.split()
	mapheader = dict()
        i=0
	for entry in headersplit:
		mapheader[entry] = i
		i+=1
		print entry,i
	for field in fields:
                minfield = fields[field][0]
                maxfield = fields[field][1]
                print field,minfield,maxfield
		outf.write("#%s\t%f\t%f\n" %(field, minfield, maxfield))
                if field in mapheader:
                	print "field found in header"
			fields[field].append(mapheader[field])
			print mapheader[field]
	for line in lines[i0+1:]:
		data = line.split()
		if line.find("I_sc") != -1: continue
		setTrue = False
		addstring = []
		for field in fields:
			minfield = fields[field][0]
			maxfield = fields[field][1]
			index = fields[field][2]
			val = float(data[index])
			if val<maxfield and val>=minfield:
				setTrue = True
				#print field,val
				addstring.append(str(val))
			else:
				setTrue = False
				#print field,"False"
				break
		if setTrue:
			pdbfilename=data[mapheader["description"]]
			print pdbfilename+".pdb"
			
			outf.write(pdbfilename+'\t'+'\t'.join(addstring)+"\n")
			outfdata.write(line)
			if createlink:
				fullpath = outputdir + "/"+pdbfilename
				print pdbfilename,filtereddir
				parentfilename = pdbfilename[:-6]
				print parentfilename
				if os.path.exists(outputdir + "/"+pdbfilename+'.pdb.gz'):
					cmd="gunzip %s.pdb.gz" %fullpath
					os.system(cmd)
				if os.path.exists(outputdir + "/"+pdbfilename+".pdb"):
					cmd='cp %s/%s.pdb %s/.' %(outputdir,pdbfilename,filtereddir)
					print cmd
					os.system(cmd)
				if os.path.exists(pdbfilename+".pdb"):
                                        cmd = 'cp %s.pdb %s/.' %(pdbfilename,filtereddir)
					print cmd
					os.system(cmd)
				if os.path.exists(parentfilename+".pdb"):
                                        cmd = 'cp %s.pdb %s/.' %(parentfilename,filtereddir)
                                        print cmd
                                        os.system(cmd)
				
	outf.close()
	outfdata.close()

def mergefiles(scorefile1,scorefile2,type="parentchild"):
	mapping1 = dict()
	mapping2 = dict()
	mapping1 = getmapping(scorefile1)
	mapping2 = getmapping(scorefile2)
	if not "description" in mapping1:
		print "description field not found"
	if not "description" in mapping2:
                print "description field not found"
	sf1 = open(scorefile1,'r')
	if type=="parentchild":
		outf=open("Mergedfile.sc",'w')
		newfields=""
		for field in mergedata:
			newfields += field + '\t'
		sf11=open(scorefile1,'r')
        	headerlines1 = sf11.readlines()[:i0+1]
        	sf11.close()
        	header = headerlines1[i0].split('\n')[0] + '\t' +newfields + '\n'
        	outf.write(headerlines1[0]+header)
		i=0
		for line1 in sf1:
			if line1.find("description")!=-1 or line1.find("SEQUENCE")!=-1:
				i+=1
				continue
			data1 = line1.split()
			filename1 = data1[mapping1["description"]]
			print filename1
			i+=1
			sf2 = open(scorefile2,'r')
			for line2 in sf2:
				if line2.find("description")!=-1 or line2.find("SEQUENCE")!=-1:
                                	continue
				data2 = line2.split()
				filename2 = data2[mapping2["description"]]
                        	print filename2
				if filename1.find(filename2):
					outline = line1.split('\n')[0]
					for field in mergedata:
						outline += '\t' + data2[mapping2[field]]
					outf.write(outline+'\n')
					break
			sf2.close()
		outf.close()
			
	sf1.close()

#fields = {"I_sc":[-160,-130], "temp_level": [1,2], "rms": [15,25]}

i0=1
parser = argparse.ArgumentParser()
parser.add_argument('--copy',default=False,action='store_true', help='copy files')
parser.add_argument('--file', help='scorefile name or graph info file name')
parser.add_argument('--plotscatter',default=False,action='store_true', help='plot graphs')
parser.add_argument('--file1', help='scorefile name')
parser.add_argument('--file2', help='scorefile name')
parser.add_argument('--mergefiles',default=False,action='store_true', help='merge files')
parser.add_argument('--printheader',default=False,action='store_true', help='print header columns')
args = parser.parse_args()
createlink = args.copy
#fields = {"interface_delta_X":[-10,-4]}
data = ["interaction_energy"]
#basename="i_sc_les_than_neg4"
fields = {"interaction_energy":[-40,-28.5],"substrate_rmsd": [0,100]}
basename="HighResolution"
outputdir = "output"
pltinfo=dict()
for field in fields:
	basename+=field+"_"
filtereddir = basename + "filtereddir"
os.system('mkdir -p %s' %filtereddir)
if args.file:
	printmapping(args.file)
if args.file1:
	print "file1 mapping"
	printmapping(args.file1)
if args.file2:
        print "file1 mapping"
        printmapping(args.file2)
if not (args.printheader or args.mergefiles or args.plotscatter):
	filterpdb(args.file,basename)
elif args.mergefiles:
	mergedata=["interaction_energy","substrate_rmsd","glycan_rmsd"]
	mergefiles(args.file1,args.file2)
if args.plotscatter:
	execfile(args.file)
	print pltinfo
	plotgraphs()
