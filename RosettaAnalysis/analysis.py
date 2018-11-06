from pyrosetta import *
from rosetta import *
pyrosetta.init('-include_sugars')


#Setting up mover for viewing structures
pm = PyMOLMover() # 1. Different from tutorial
import argparse
from rosetta.core.pose import pose_from_saccharide_sequence
from rosetta.core.scoring import *


def CalculateRMSDCA(pose_mod,pose_ref):
	print "rmsd ca"
	rmsd = CA_rmsd(pose_mod,pose_ref)
	print rmsd

def CalculateRMSDProteinOnly():
	print "rmsd protein only"

def CalculateRMSDHeavyAtoms(pose_mod,pose_ref):
	print "rmsd heavy atoms"
	rmsd = rms_at_corresponding_heavy_atoms(pose_mod,pose_ref)
	print rmsd

def CalculateRMSDSelection():
	print "rmsd selection"

def CalculateRMSDAll():
	CalculateRMSDCA()
	CalculateRMSDProteinOnly()

def GetRMSDForFile(filename,pose_ref,rmsdtype):
	print "Getting rmsd for file", filename
	if os.path.exists(filename):
		pose_compare = pose_from_pdb(filename)
	if rmsdtype=='rmsd_ca':
                CalculateRMSDCA(pose_compare,pose_ref)
        if rmsdtype=='rmsd_ha':
                CalculateRMSDHeavyAtoms(pose_compare,pose_ref)	

def GetRMSDForList(filelist,pose_ref,rmsdtype):
	f=open(filelist,'r')
	for line in f:
		filename=line.split()[0]
		if os.path.exists(filename):
			GetRMSDForFile(filename,pose_ref,rmsdtype)
			 

if __name__ ==  "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--type_calc', help='type of rmsd calculation, options are rmsd_ca')
	parser.add_argument('--native', help='reference pdb',default='native.pdb')
	parser.add_argument('--traj', action='store_false', help='process list of files')
	parser.add_ardument('--file', action='store_true', help='process single file')
	parser.add_argument('--singlefile', help='pdb file name for rmsd calculation', default='input.pdb')
	parser.add_argument('--filelist', help='list of pdb files for rmsd calculation', default='list.txt')
	
	args = parser.parse_args()
	print args
	pose_ref = pose()
	if os.path.exists(args.native):
		print "Native file: ",args.native, " used for comparison"
		pose_ref = pose_from_pdb(args.native)
	else:
		print "Cannot calcualte rmsd without native structure"
		exit()
	if args.file:
		filename=args.singlefile
		if os.path.exists(singlefile):
			GetRMSDForFile(filename,pose_ref,args.type_calc)
		else:
			print "provide filename for calculation with --singlefile option. Use --help for more details."
	if args.traj:
		filelist=args.filelist
		if os.path.exists(filelist):
			GetRMSDForList(filelist,pose_ref,args.type_calc)
		else:
			print "provide filename for list of pdbs with --filelist option. Use --help for more details."
