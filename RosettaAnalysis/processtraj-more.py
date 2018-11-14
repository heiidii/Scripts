import sys
import os
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import glob
import mpl_toolkits
#.basemap import Basemap as Basemap
import cv2

basename= "Traj0split"
command= "ProteinConnectivityForSampling.linuxgccdebug"
n1=1
n2=102
reslist = [15]#,188] #W138, Cys188

def cleanfiles():
 for i in range(1,101):
   filename = "%s%04d" %(basename,i)
   print filename
   cmd  = "~/RosettaDevelopment/Rosetta/tools/protein_tools/scripts/clean_pdb.py %s A " %(filename)
   print cmd
   os.system(cmd)

def getgraphs():
 for i in range(1,102):
   filename = "%s%04d_A.pdb" %(basename,i)
   print filename
   cmd = "%s -in:file:native %s " %(command, filename)
   print cmd
   os.system(cmd)

def getxy(res):
	y=res%10
        xtemp=res-y
        if xtemp>10:
        	x=xtemp/10
        else:
        	x=0
	return x,y

def makemovieforres(res):
 for i in range(n1,n2):
   pngname = "Fig_%s_%d" %(basename,res)
   print pngname
   #cmd="ffmpeg -r %d -i %s_\%04d.png -vcodec mpeg4 -y movie_%d_fr%d.mp4 " %(framerate,pngname,res,framerate)

def makemovieforimages(basenamestr,i1,i2,step,outname):
	vid = None #cv2.VideoWriter(outname,-1,1,(width,height))
	size=None
	format = "XVID"
	for j in range(i1,i2+1,step):
		imgname = basenamestr %j
		if not os.path.exists(imgname): continue
		print imgname
                img = cv2.imread(imgname)
		if vid is None:
                	if size is None:
                    		size = img.shape[1], img.shape[0]
                vid = cv2.VideoWriter(outname, cv2.VideoWriter_fourcc(*format), float(5), size)
            	if size[0] != img.shape[1] and size[1] != img.shape[0]:
                	img = resize(img, size)
		vid.write(img)

	cv2.destroyAllWindows()
	vid.release()


def getneighbors(res):
 for i in range(n1,n2):
   filename = "PrunedSym_%s%04d_A.out" %(basename,i)
   print filename
   outname = "Fig_%s_%d_%04d.png" %(basename,res,i)
   getneighborsforfile(res,filename,outname)

def getneighborsforfile(res,filename,outname):
   Gnew = nx.Graph()
   print filename
   if os.path.exists(filename):
	G=nx.read_edgelist(filename,nodetype=int)
        L=G.neighbors(res)
        lneighbors = list(L)
        #print lneighbors.sort()
        G.add_node(res)
        G.nodes[res]['neighborsize']=len(lneighbors)
	posdict = dict()
	sizedict = dict()
        x,y = getxy(res)
	posdict[res]=[x,y]
	sizedict[res]=len(lneighbors)
        #print posdict[res]
	for resnn in lneighbors:
		Lnext = G.neighbors(resnn)
                lneighborsnext = list(Lnext)
                Gnew.add_edge(res,resnn)
		x,y = getxy(resnn)
		posdict[resnn]=[x,y]
		sizedict[resnn]=len(lneighborsnext)
		#print posdict[resnn]
	Gnew.nodes[resnn]['neighborsize']=len(lneighborsnext)
	#print "Drawing network"
	nx.draw_networkx(Gnew,pos=posdict,with_labels=True,font_size=11)
	#nodelist=sizedict.keys(),node_size=[ v/(max(sizedict.values()))*100 for v in sizedict.values()])
	#plt.ylim(-1,11)
	#plt.xlim(-2,27)
	plt.xticks([])
        plt.yticks([])
	plt.savefig(outname,dpi=300)
	#plt.savefig("Fig_%s_%d_%04d.png" %(basename,res,i), dpi=300)
	#plt.xticks([])
	#plt.yticks([])
	#plt.show()
	Gnew.clear()
	plt.close()
        plt.clf()
   else:
	print "file does not exist ",filename

def readweightedgraph(filename):
   G = nx.Graph()
   #print filename
   if os.path.exists(filename):
        G=nx.read_weighted_edgelist(filename,nodetype=int)
   else: 
	print "File not found"
   return G

def getprunedneighbors(G,res,distThreshold):
	nneighbors = 0
	if G.has_node(res):
        	L=G.neighbors(res)
        	lneighbors = list(L)
        	for resnn in lneighbors:
			if G[res][resnn]['weight']>distThreshold: continue
			nneighbors +=1
	return nneighbors

def getprunedgraph(Gw,distThreshold):
        nneighbors = 0
	Gpruned = nx.Graph()
        for n,nbrs in Gw.adjacency():
		for nbr,eattr in nbrs.items():
			edgeweight = eattr['weight']
			if edgeweight <= distThreshold:
                        	Gpruned.add_edge(n,nbr,weight=edgeweight)
        return Gpruned


#getgraphs()
framerate=2
singlefile="PrunedSym_Bound.out"

def gethistogramforfile(file,bins='auto',res=-1):
   f=open(file,'r')
   lines=f.readlines()
   f.close()
   nlines = len(lines)
   valarray = np.zeros((nlines,1))
   i=0
   for line in lines:
      values = [float(t) for t in line.split()]
      if res!=-1:
	if values[0]!=float(res) or values[1]!=float(res): continue
      valarray[i,0]=values[2]
      i+=1
   n, bins, patches = plt.hist(x=valarray,bins=bins, color='#0504aa',alpha=0.7, rwidth=0.85)
   print n
   return n,bins

def getdistancesforresidue(file,res,bins='auto'):
   f=open(file,'r')
   lines=f.readlines()
   f.close()
   nlines = len(lines)
   valarray = np.zeros((nlines,1))
   vallist=[]
   i=0
   for line in lines:
      values = [float(t) for t in line.split()]
      #print values
      if values[0]==float(res) or values[1]==float(res):
      	vallist.append(values[2])
      	i+=1
   
   n, bins, patches = plt.hist(x=vallist,bins=bins, color='#0504aa',alpha=0.7, rwidth=0.85)
   #print n,bins
   return vallist,i,n,bins

def getdistancehistogramforresidue(files,res):
   #hist, center = np.histogram()
   valarray, numvals,sumfreq,bins = getdistancesforresidue(files[0],res)
   sumsqfreq = np.square(sumfreq)
   sumdistance = np.sum(valarray)
   sumsqdistance = np.sum(np.square(valarray))
   sumedgecount = numvals
   #print sumfreq,bins
   for file in files[1:]:
        valarray, numvals,freq,bins = getdistancesforresidue(file,res,bins)
        sumfreq += freq
        sumsqfreq += np.square(freq)
        sumdistance += np.sum(valarray)
        sumsqdistance += np.sum(np.square(valarray))
        sumedgecount += numvals
        #print freq,bins
   fN =  float(len(files))
   avgfreq = sumfreq/fN
   sdfreq = np.sqrt(abs(np.square(sumfreq) - fN*sumsqfreq))/fN
   print avgfreq
   print sdfreq
   avgdistance = sumdistance/float(sumedgecount)
   sddistance = np.sqrt(abs(np.square(sumdistance) - float(sumedgecount)*sumsqdistance))/float(sumedgecount)
   plt.clf()
   plt.xlabel('Distance')
   plt.ylabel('Frequency')
   centers = (bins[1:] + bins[:-1])*0.5
   plt.errorbar(centers,avgfreq,yerr=sdfreq,color='xkcd:sky blue',marker='o',ecolor='xkcd:pink',capthick=4,capsize=2)
   #plt.show()
   plt.savefig("DistanceDistribution_Res" + str(res)+".png" ,dpi=600)
   print res,avgdistance,sddistance   

def getdistancehistogram(files):
   #hist, center = np.histogram()
   sumfreq,bins = gethistogramforfile(files[0])
   sumsqfreq = np.square(sumfreq)
   for file in files[1:]:
   	freq,bins = gethistogramforfile(file,bins)
   	sumfreq += freq
	sumsqfreq += np.square(freq)
   fN = float(len(files))
   avgfreq = sumfreq/fN
   sdfreq = np.sqrt(abs(np.square(sumfreq) - fN*sumsqfreq))/fN
   print avgfreq
   print sdfreq
   centers = (bins[1:] + bins[:-1])*0.5
   plt.clf()
   plt.xlabel('Distance')
   plt.ylabel('Frequency')
   #plt.plot(centers,avgfreq,color)
   #plt.bar(centers,avgfreq, alpha=0.65)
   plt.errorbar(centers,avgfreq,yerr=sdfreq,color='xkcd:sky blue',marker='o',ecolor='xkcd:pink',capthick=4,capsize=2)
   plt.savefig("DistanceDistribution.png",dpi=600)


def GetNegihborsForResidue():
   for res in reslist:
	outfile = "DistancePruned_%d.png" %res
	#getneighborsforfile(res,singlefile,outfile)
        getneighbors(res)

def getmapping(maplist):
	newlist=[]
	for entry in maplist:
		if entry in mappednames:
			newlist.append(mappednames[entry])
		else:
			print "Entry not found",entry
			newlist.append(str(entry))
	return newlist

def plotresiduebargraph(xvals,yvals,pltinfo,bshow=False):
	plt.title(pltinfo["title"])
        plt.xlabel(pltinfo["xlabel"])
        plt.ylabel(pltinfo["ylabel"])
        
        plt.bar(range(len(yvals)),yvals,color='xkcd:sky blue')#,marker='o',ecolor='xkcd:lavendar',capthick=4,capsize=2)
        mappedlist = getmapping(xvals)
        plt.xticks(range(len(xvals)),mappedlist,fontsize=7,rotation=70)
	if bshow: plt.show()
        else: plt.savefig(pltinfo["outfile"],dpi=600)
	plt.close()
	plt.clf()

def plotresidueplotgraph(xvals,yvals,pltinfo,bshow=False):
        plt.title(pltinfo["title"])
        plt.xlabel(pltinfo["xlabel"])
        plt.ylabel(pltinfo["ylabel"])

        plt.plot(range(len(yvals)),yvals,color='xkcd:sky blue',marker='o')#,ecolor='xkcd:lavendar',capthick=4,capsize=2)
        mappedlist = getmapping(xvals)
        plt.xticks(range(len(xvals)),mappedlist,fontsize=7,rotation=70)
        if bshow: plt.show()
        else: plt.savefig("plt"+pltinfo["outfile"],dpi=600)
        plt.close()
        plt.clf()

def plotresidueerrorgraph(xvals,yvals,yerr,pltinfo,bshow=False):
        plt.title(pltinfo["title"])
        plt.xlabel(pltinfo["xlabel"])
        plt.ylabel(pltinfo["ylabel"])

        plt.errorbar(range(len(yvals)),yvals,yerr=yerr,color='xkcd:sky blue',marker='o',ecolor='xkcd:pink',capthick=4,capsize=2)
        mappedlist = getmapping(xvals)
        plt.xticks(range(len(xvals)),mappedlist,fontsize=7,rotation=70)
        if bshow: plt.show()
        else: plt.savefig("yerr"+pltinfo["outfile"],dpi=600)
        plt.close()
        plt.clf()

def getmappingpdb(inlist):
	newlist=[]
	for entry in inlist:
		newlist.append(getmappingpdbentry(entry))
	return newlist

def getmappingpdbentry(inres):
	return inres+offset_ros2pdb

def writepropstofile(xvals,yvals,info):
	print "Writing file",info["outfile"]
	f=open(info["outfile"],'w')
        mappedlist = getmappingpdb(xvals)
	#f.write("#"+info["title"]+"\n")
	for i in range(0,len(yvals)):
		f.write(str(xvals[i])+'\t'+str(yvals[i])+'\t'+str(mappedlist[i])+'\n')
	f.close()

def getresiduesinrange(reslist,vals,low,high):
	range_reslist = []
	range_vals = []
	for j in range(0,len(vals)):
                if vals[j]>=low and vals[j]<high:
                        range_vals.append(vals[j])
                        range_reslist.append(reslist[j])
	return range_reslist,range_vals

def getresiduesinrange_2(reslist,vals,sdlist,low,high):
        range_reslist = []
        range_vals = []
	range_sdlist = []
        for j in range(0,len(vals)):
                if vals[j]>=low and vals[j]<high:
                        range_vals.append(vals[j])
                        range_reslist.append(reslist[j])
			range_sdlist.append(sdlist[j])
        return range_reslist,range_vals,range_sdlist

def GetDistancePrunedNeighbors(files,distThreshold,sdmorethan,reslist=range(1,240)):
	sum_nn_list= [0.0 for t in range(0,len(reslist))]
        count_nn_list = [0.0 for t in range(0,len(reslist))]
	sumsq_nn_list= [0.0 for t in range(0,len(reslist))]
	avg_nn_list= [0.0 for t in range(0,len(reslist))]
	sd_nn_list= [0.0 for t in range(0,len(reslist))]
        for curfile in files:
		G = readweightedgraph(curfile)
		for res in reslist:
			indres = reslist.index(res)
			nneighbors = getprunedneighbors(G,res,distThreshold)
			sum_nn_list[indres]  += nneighbors
			sumsq_nn_list[indres] += nneighbors*nneighbors
			count_nn_list[indres] += 1
	for i in range(0,len(reslist)):
		if count_nn_list[i] !=0:
			avg_nn_list[i] = sum_nn_list[i]/count_nn_list[i]
			sd_nn_list[i] = (abs(sum_nn_list[i]*sum_nn_list[i] - count_nn_list[i] * sumsq_nn_list[i]))**0.5/count_nn_list[i]
		else:
			avg_nn_list[i]=0
			sd_nn_list[i]=0
	sortedavg_indices = np.argsort(avg_nn_list)
	sortedsd_indices = np.argsort(sd_nn_list)
	sortedavg_reslist = [reslist[t] for t in sortedavg_indices]
	sortedavg_avg = [avg_nn_list[t] for t in sortedavg_indices]
	sortedavg_sd = [sd_nn_list[t] for t in sortedavg_indices]
	sortedsd_reslist = [reslist[t] for t in sortedsd_indices]
        sortedsd_avg = [avg_nn_list[t] for t in sortedsd_indices]
        sortedsd_sd = [sd_nn_list[t] for t in sortedsd_indices]
        pltinfo = dict()
        pltinfo["xlabel"]="Residue"
        pltinfo["ylabel"]=r'Number of neighbors within %3.1f $\AA$' %distThreshold
	pltinfo["title"]=""
        pltinfo["outfile"]="SortedAllresidues_distance%3.1f.png" %distThreshold
        plotresiduebargraph(sortedavg_reslist,sortedavg_avg,pltinfo)
	plotresidueerrorgraph(sortedavg_reslist,sortedavg_avg,sortedavg_sd,pltinfo)
	pltinfo.clear()

	info = dict()
        info["outfile"]="%sAvg_distance%3.1f.txt" %(prefix,distThreshold)
        writepropstofile(sortedavg_reslist,sortedavg_avg,info)

	alllimits=[[min(sortedavg_avg),6],
		   [6,8],
		   [8,10],
		   [10,12],
		   [12,14],
		   [14,max(sortedavg_avg)]]
	for limits in alllimits:
		range_reslist , range_avg , range_sd= getresiduesinrange_2(sortedavg_reslist,sortedavg_avg,sortedavg_sd,limits[0],limits[1])
		pltinfo["xlabel"]="Residue"
		pltinfo["ylabel"]=r'Average Number of Neighbors between %d and %d within %3.1f $\AA$' %(limits[0],limits[1],distThreshold)
		pltinfo["title"]=""
		pltinfo["outfile"]="RangeAvgNN_%d-%d_distance%3.1f.png" %(limits[0],limits[1],distThreshold)
		plotresiduebargraph(range_reslist,range_avg,pltinfo)
        	plotresidueerrorgraph(range_reslist,range_avg,range_sd,pltinfo)
        	pltinfo.clear()

	info = dict()
        info["outfile"]="%sSD_distance%3.1f.txt" %(prefix,distThreshold)
        writepropstofile(sortedsd_reslist,sortedsd_sd,info)
	#SD sorted
	low=sdmorethan
	high=max(sortedsd_sd)
        sortedsd_reslist_morethan, sortedsd_sd_morethan = getresiduesinrange(sortedsd_reslist,sortedsd_sd,low,high)
	pltinfo["xlabel"]="Residue"
	pltinfo["ylabel"]=r'SD of neighbors within %3.1f $\AA$' %distThreshold
	pltinfo["title"]="SD > %3.1f " %sdmorethan
	pltinfo["outfile"]="%sSD_morethan%3.1f_distance%3.1f.png" %(prefix,sdmorethan,distThreshold)
	plotresiduebargraph(sortedsd_reslist_morethan,sortedsd_sd_morethan,pltinfo)
	plotresidueplotgraph(sortedsd_reslist_morethan,sortedsd_sd_morethan,pltinfo)

def WriteGraphForPymol(Gv,outname):
        nodeall = [node for node in Gv.nodes(data=True)]
	outf = open(outname, 'w')
	for node in nodeall:
		#print node[0]
		curlabel = getmappingpdbentry(node[0])
		outf.write(str(node[0])+"\t"+str(curlabel)+"\n")
		print str(node[0])+"\t"+str(curlabel)+"\n"
	outf.close()


def VisualizeGraph(Gv,outname,bLabels=False,bshow=False,rescur=-1):
	elarge = [(u, v) for (u, v, d) in Gv.edges(data=True) if d['weight'] > 4.0]
	esmall = [(u, v) for (u, v, d) in Gv.edges(data=True) if d['weight'] <= 4.0]
	#for res in reslist:
        #	indres = reslist.index(res)
        #	nneighbors = getprunedneighbors(G,res,distThreshold)
	pos = nx.spring_layout(Gv)  # positions for all nodes
	plt.clf()

	# nodes
	ncolor='r'
	#ncolor=['red' for t in range(0,len())]
	if rescur != -1:
		ncolor = []
		for node in Gv:
			if rescur==node:
				ncolor.append('blue')
			else:
				ncolor.append('red')
	nx.draw_networkx_nodes(Gv, pos,node_color=ncolor, node_size=200)

	# edges
	nx.draw_networkx_edges(Gv, pos, edgelist=elarge,width=4,alpha=0.5)
	nx.draw_networkx_edges(Gv, pos, edgelist=esmall,width=4, alpha=0.5, edge_color='b', style='dashed')

	# labels
	if bLabels:
		curlabels = dict()
		for key in pos:
			if key in mappednames:
				curlabels[key] = mappednames[key]
			else:
				curlabels[key] = str(key)
		nx.draw_networkx_labels(Gv, pos, font_size=4,labels=curlabels, font_family='sans-serif')

	plt.axis('off')
	if bshow:
		plt.show()
	else:
		plt.savefig(outname,dpi=600)
	#plt.close()
	#plt.clf()
	

def CreateGlobalMapping():
	f=open(mappingfile,'r')
	lines=f.readlines()
	f.close()
	for line in lines:
		mappednames[int(line.split()[0])]=line.split()[2]
	return mappednames

#distfile = sys.argv[1]
#files=glob.glob('Distance*.out')
#getdistancehistogram(files)
residue=15
files = ['DistancesTraj0split0001_A.out']#,'DistancesTraj0split0002_A.out','DistancesTraj0split0004_A.out','DistancesTraj0split0100_A.out']
#getdistancehistogramforresidue(files,residue)

distThlist=[6.0,7.0,8.0,9.0,10.0]
sdThlist=[1.5,2.0,2.5]

offset_ros2pdb=+1
prefix ="md"
mappingfile='resmappingmd.out'
mappednames = dict()
mappednames = CreateGlobalMapping()
G=nx.Graph()

basename = "DistancesTraj0split"
i1=2
i2=102
step=2
files = ["%s%04d_A.out" %(basename,i) for i in range(i1,i2,step)]

ressubgraph=[137,214,138] #GLY_140
bSF = False
singlefile = 'Distances1ACB_b_start.out'
outnamesinglefile = "1ACB_b"
for dist in [5.0]:
   for res in ressubgraph:
        if bSF:
		print res
		curfile = singlefile
         	Gv = readweightedgraph(curfile)
         	Gpruned = getprunedgraph(Gv,dist)
         	outname = "NetworkPrunedLabeled_%2.1f_%s.png" %(dist,outnamesinglefile)
         	#VisualizeGraph(Gpruned,outname,True,True)
         	Gcc = sorted(nx.connected_component_subgraphs(Gpruned), key=len, reverse=True)

         	for k in range(0,3):
                	outname = "Subnetwork%d_%2.1f_%s" %(k,dist,outnamesinglefile)
                	#VisualizeGraph(Gcc[k],outname + ".png",True)
			WriteGraphForPymol(Gcc[k],outname + ".txt")
         	outname = "SubnetworkRes%d_%2.1f_%s" %(res,dist,outnamesinglefile)
         	m=0
         	for sgraph in Gcc:
                	if sgraph.has_node(res):
                        	#print i,m,len(sgraph)
                        	print str(i)+'\t'+str(len(sgraph))+'\t'+str(m)+'\n'
                        	#VisualizeGraph(sgraph,outname + ".png",True,False,res)
				WriteGraphForPymol(Gcc[k],outname + ".txt")
                	m+=1
		continue
	print res
	outsgraph = open("SubnetworkRes%d_%3.1f.txt" %(res,dist),'w')
	outsgraph.write('#image\tsubgraphlen\tsubgraphid\n')
        outsgraph0 = open("Subnetwork0_%3.1f.txt" %(dist),'w')
        outsgraph0.write('#image\tsubgraphlen\tsubgraphid\n')
	for i in range(i1,i2,step):
	 print i
	 curfile = "%s%04d_A.out" %(basename,i)
	 Gv = readweightedgraph(curfile)
	 Gpruned = getprunedgraph(Gv,dist)
	 #outname = "NetworkPruned_%2.1f_%04d.png" %(dist,i)
	 #VisualizeGraph(Gpruned,outname,False)
	 Gcc = sorted(nx.connected_component_subgraphs(Gpruned), key=len, reverse=True)
    	  
	 for k in range(0,3):
	 	outname = "Subnetwork%d_%2.1f_%04d" %(k,dist,i)
		if k==0:
			outsgraph0.write(str(i)+'\t'+str(len(Gcc[k]))+'\t'+str(k)+'\n')
         	#VisualizeGraph(Gcc[k],outname + ".png",True)
		WriteGraphForPymol(Gcc[k],outname + ".txt")
	 outname = "SubnetworkRes%d_%2.1f_%04d" %(res,dist,i)
         m=0
	 for sgraph in Gcc:
		if sgraph.has_node(res):
			#print i,m,len(sgraph)
			outsgraph.write(str(i)+'\t'+str(len(sgraph))+'\t'+str(m)+'\n')
         		#VisualizeGraph(sgraph,outname+ ".png",True,False,res) 
			WriteGraphForPymol(Gcc[k],outname + ".txt")
	 	m+=1
	#makemovieforimages("NetworkPruned_%2.1f" %dist + "_%04d.png",i1,i2,step,"MovingNetwork_%2.1f_%04d-%04d.avi" %(dist,i1,i2))
	outsgraph.close()
	outsgraph0.close()

for dt in distThlist:
	for sdt in sdThlist:
		continue
 		GetDistancePrunedNeighbors(files,dt,sdt)
