from pyrosetta import *
from rosetta import *
pyrosetta.init('-include_sugars')


#Setting up mover for viewing structures
pm = PyMOLMover() # 1. Different from tutorial
import argparse
from rosetta.core.pose import *
# pose_from_saccharide_sequence
from rosetta.core.scoring import *


def CalculateRMSDCA(pose_mod,pose_ref):
	#print "rmsd ca"
	rmsd = CA_rmsd(pose_mod,pose_ref)
	#print rmsd
	return rmsd

def CalculateRMSDProteinOnly():
	print "rmsd protein only"

def CalculateRMSDHeavyAtoms(pose_mod,pose_ref):
	#print "rmsd heavy atoms"
	rmsd = rms_at_corresponding_heavy_atoms(pose_mod,pose_ref)
	#print rmsd
	return rmsd

def CalculateRMSDSelection():
	print "rmsd selection"

def CalculateRMSDAll():
	rmsdca = CalculateRMSDCA()
	rmsdha = CalculateRMSDProteinOnly()

def GetRMSDForFile(filename,pose_ref,rmsdtype):
	#print "Getting rmsd for file", filename
	if os.path.exists(filename):
		pose_compare = pose_from_pdb(filename)
		rmsd='-1.0'
		if rmsdtype=='rmsd_ca':
                	rmsd = CalculateRMSDCA(pose_compare,pose_ref)
        	if rmsdtype=='rmsd_ha':
                	rmsd = CalculateRMSDHeavyAtoms(pose_compare,pose_ref)
		return rmsd

def GetRMSDForList(filelist,pose_ref,rmsdtype,fout):
	f=open(filelist,'r')
	for line in f:
		filename=line.split()[0]
		if os.path.exists(filename):
			rmsd = GetRMSDForFile(filename,pose_ref,rmsdtype)
			fout.write(filename+'\t%06.4f\n' %rmsd)

if __name__ ==  "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--type_calc', help='type of rmsd calculation, options are rmsd_ca')
	parser.add_argument('--native', help='reference pdb',default='native.pdb')
	parser.add_argument('--traj', action='store_false', help='process list of files')
	parser.add_argument('--fileonly', action='store_true', help='process single file')
	parser.add_argument('--singlefile', help='pdb file name for rmsd calculation', default='input.pdb')
	parser.add_argument('--filelist', help='list of pdb files for rmsd calculation', default='list.txt')
	
	args = parser.parse_args()
	print args
	
	pose_ref = Pose()
	if os.path.exists(args.native):
		print "Native file: ",args.native, " used for comparison"
		pose_ref = pose_from_pdb(args.native)
	else:
		print "Cannot calcualte rmsd without native structure"
		exit()
	if args.fileonly:
		filename=args.singlefile
		if os.path.exists(filename):
			rmsd = GetRMSDForFile(filename,pose_ref,args.type_calc)
			print "\nNative: ",args.native, "\nComparison: ",args.singlefile," \nRMSD: ",rmsd
		else:
			print "provide filename for calculation with --singlefile option. Use --help for more details."
	elif args.traj:
		filelist=args.filelist
		if os.path.exists(filelist):
			outfilename = filelist + "_rmsd.out"
			fout = open(outfilename,'w')
			fout.write('#%s\n' %args)
			GetRMSDForList(filelist,pose_ref,args.type_calc,fout)
			fout.close()
		else:
			print "provide filename for list of pdbs with --filelist option. Use --help for more details."
	else:
		print "option does not exist"
