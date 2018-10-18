import os
import sys




def getpdbsbyrmsd():
	print "filtering by rmsds"
	

def getpdbsbytemperature():
	print "filtering by temeprature"

def getpdbsbyscore():
	print "filtering by score"

def printmapping(sfile):
	print "mapping %s ..." %sfile
        f = open(sfile,'r')
        lines = f.readlines()
        f.close()
        outf = open("MapHeader.txt",'w')
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
                outf.write("#%s\t%d\n" %(entry, i))

	outf.close()

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
			print field,minfield,maxfield,index,val
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
			#print pdbfilename
			
			outf.write(pdbfilename+'\t'+'\t'.join(addstring)+"\n")
			outfdata.write(line)
	outf.close()
	outfdata.close()


#fields = {"I_sc":[-160,-130], "temp_level": [1,2], "rms": [15,25]}

i0=1
fields = {"interaction_energy":[-40,-20]}
data = ["interaction_energy"]
basename=""
for field in fields:
	basename+=field+"_"
#printmapping(sys.argv[1])
filterpdb(sys.argv[1],basename)
