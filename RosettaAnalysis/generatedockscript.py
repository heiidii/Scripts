import os
import sys
import argparse

templatepdbflags='''
-in
        -file
                -s T12_GalNAc.pdb
                -extra_res_fa NGA_ideal.params
-out:path:pdb /home/aferna45/T12/docking/T12_docking/global_docking/randoms/out_%03d/

-packing
        -ex1
        -ex2
        -no_optH false
        -flip_HNQ true
        -ignore_ligand_chi true





-multiple_processes_writing_to_one_directory true
-mistakes
        -restore_pre_talaris_2013_behavior true
-nstruct 1
'''
templateflags='''
-in
        -file
                -s startfrom.pdb
                -extra_res_fa NGA_ideal.params
-out:path:pdb /home/aferna45/T12/docking/T12_docking/global_docking/randoms/out_%03d/output_pdbs
-out:pdb_gz
-out:file:scorefile random_global_docking.sc
        
-packing
        -ex1
        -ex2
        -no_optH false
        -flip_HNQ true
        -ignore_ligand_chi true





-multiple_processes_writing_to_one_directory true
-in:file:native startfrom.pdb
-mistakes
        -restore_pre_talaris_2013_behavior true
-nstruct 2000
'''
templateDump='''
<ROSETTASCRIPTS>

                <SCOREFXNS>
                        <ScoreFunction name="ligand_soft_rep" weights="ligand_soft_rep">
                        </ScoreFunction>
                        <ScoreFunction name="hard_rep" weights="ligand">
                        </ScoreFunction>
                </SCOREFXNS>

                <RESIDUE_SELECTORS>
                        <Chain name="Protein_and_Glycan" chains="A,X" />
                        <Chain name="Glycan" chains="X"/>
                        <Chain name="ProteinOnly" chains="A" />
                </RESIDUE_SELECTORS>

                <LIGAND_AREAS>
                        <LigandArea name="inhibitor_dock_sc" chain="X" cutoff="6.0" add_nbr_radius="true" all_atom_mode="false"/>
                        <LigandArea name="inhibitor_final_sc" chain="X" cutoff="6.0" add_nbr_radius="true" all_atom_mode="false"/>
                        <LigandArea name="inhibitor_final_bb" chain="X" cutoff="7.0" add_nbr_radius="false" all_atom_mode="true" Calpha_restraints="0.3"/>
                </LIGAND_AREAS>

                <INTERFACE_BUILDERS>
                        <InterfaceBuilder name="side_chain_for_docking" ligand_areas="inhibitor_dock_sc"/>
                        <InterfaceBuilder name="side_chain_for_final" ligand_areas="inhibitor_final_sc"/>
                        InterfaceBuilder name="backbone" ligand_areas="inhibitor_final_bb" extension_window="3"/>
                </INTERFACE_BUILDERS>

                <MOVEMAP_BUILDERS>
                        <MoveMapBuilder name="docking" sc_interface="side_chain_for_docking" minimize_water="false"/>
                        <MoveMapBuilder name="final" sc_interface="side_chain_for_final"  minimize_water="false"/>
                </MOVEMAP_BUILDERS>

                <SCORINGGRIDS ligand_chain="X" width="15">
                        <ClassicGrid grid_name="classic" weight="1.0"/>
                </SCORINGGRIDS>

                <MOVERS>
                        <StartFrom name="start" chain="X">
                                <Coordinates x="%s" y="%s" z="%s"/>
                        </StartFrom>
                        <DumpPdb name="dumppdb" fname="startfrom.pdb" scorefxn="ligand_soft_rep"/>
                        <DumpPdb name="dumppdb2" fname="transform1.pdb" scorefxn="ligand_soft_rep"/>
                        <Transform name="transform" chain="X" box_size="20.0" move_distance="0.2" angle="20" cycles="500" repeats="1" temperature="5" initial_perturb="1"/>
                        <HighResDocker name="high_res_docker" cycles="12" repack_every_Nth="3" scorefxn="ligand_soft_rep" movemap_builder="docking"/>
                        <FinalMinimizer name="final" scorefxn="hard_rep" movemap_builder="final"/>

		 </MOVERS>

                <PROTOCOLS>
                        <Add mover_name="start"/>
                        <Add mover_name="dumppdb"/>
                </PROTOCOLS>


</ROSETTASCRIPTS>
'''
templateGlobal='''
<ROSETTASCRIPTS>

                <SCOREFXNS>
                        <ScoreFunction name="ligand_soft_rep" weights="ligand_soft_rep">
                        </ScoreFunction>
                        <ScoreFunction name="hard_rep" weights="ligand">
                        </ScoreFunction>
                </SCOREFXNS>

                <RESIDUE_SELECTORS>
                        <Chain name="Protein_and_Glycan" chains="A,X" />
                        <Chain name="Glycan" chains="X"/>
                        <Chain name="ProteinOnly" chains="A" />
                </RESIDUE_SELECTORS>

                <SIMPLE_METRICS>
                        <RMSDMetric name="rmsd" custom_type="default"  residue_selector="Protein_and_Glycan" use_native="true"/>
                        <RMSDMetric name="rmsdligand" custom_type="ligand_ha" rmsd_type="rmsd_all_heavy" residue_selector="Glycan" use_native="true"/>
                        <RMSDMetric name="rmsdprotein_sc" custom_type="protein_sc_ha" rmsd_type="rmsd_sc_heavy" residue_selector="ProteinOnly" use_native="true"/>
                </SIMPLE_METRICS>

                <LIGAND_AREAS>
                        <LigandArea name="inhibitor_dock_sc" chain="X" cutoff="6.0" add_nbr_radius="true" all_atom_mode="false"/>
                        <LigandArea name="inhibitor_final_sc" chain="X" cutoff="6.0" add_nbr_radius="true" all_atom_mode="false"/>
                        <LigandArea name="inhibitor_final_bb" chain="X" cutoff="7.0" add_nbr_radius="false" all_atom_mode="true" Calpha_restraints="0.3"/>
                </LIGAND_AREAS>

                <INTERFACE_BUILDERS>
                        <InterfaceBuilder name="side_chain_for_docking" ligand_areas="inhibitor_dock_sc"/>
                        <InterfaceBuilder name="side_chain_for_final" ligand_areas="inhibitor_final_sc"/>
                        InterfaceBuilder name="backbone" ligand_areas="inhibitor_final_bb" extension_window="3"/>
                </INTERFACE_BUILDERS>

                <MOVEMAP_BUILDERS>
                        <MoveMapBuilder name="docking" sc_interface="side_chain_for_docking" minimize_water="false"/>
                        <MoveMapBuilder name="final" sc_interface="side_chain_for_final"  minimize_water="false"/>
                </MOVEMAP_BUILDERS>

                <SCORINGGRIDS ligand_chain="X" width="15">
                        <ClassicGrid grid_name="classic" weight="1.0"/>
                </SCORINGGRIDS>

                <MOVERS>
                        <StartFrom name="start" chain="X">
                                <Coordinates x="%s" y="%s" z="%s"/>
                        </StartFrom>
                        <DumpPdb name="dumppdb" fname="startfrom.pdb" scorefxn="ligand_soft_rep"/>
                        <DumpPdb name="dumppdb2" fname="transform1.pdb" scorefxn="ligand_soft_rep"/>
                        <Transform name="transform" chain="X" box_size="20.0" move_distance="0.2" angle="20" cycles="500" repeats="1" temperature="5" initial_perturb="1"/>
                        <HighResDocker name="high_res_docker" cycles="12" repack_every_Nth="3" scorefxn="ligand_soft_rep" movemap_builder="docking"/>
                        <FinalMinimizer name="final" scorefxn="hard_rep" movemap_builder="final"/>
                        <InterfaceScoreCalculator name="add_scores" chains="X" scorefxn="hard_rep" native="startfrom.pdb"/>
                        <RunSimpleMetrics name="rmsdmetrics" metrics="rmsd,rmsdligand,rmsdprotein_sc" />
                </MOVERS>

                <PROTOCOLS>
                        Add mover_name="start"/>
                        Add mover_name="dumppdb"/>
                        <Add mover_name="transform"/>
                        <Add mover_name="dumppdb2"/>
                        <Add mover_name="high_res_docker"/>
                        <Add mover_name="final"/>
                        <Add mover_name="add_scores"/>
                        <Add mover_name="rmsdmetrics"/>
                        <Add mover_name="add_scores"/>
                </PROTOCOLS>


</ROSETTASCRIPTS>

'''

def writefilewithcoordinates(x,y,z,templ,filename):
	f=open(filename,'w')
	f.write(templ %(x,y,z))
	f.close()


def writefileswithcoordinatesfromfile(inputfile,basename):
	f=open(inputfile,'r')
	lines=f.readlines()
	f.close()
	outfilename = basename
	i=0
	for line in lines:
		
		dat = line.split()
		if len(dat)==3:
			i+=1		
			dirpath = '%s_%03d' %(basename,i)
			if mkdir:
				os.system('mkdir -p %s' %dirpath)
			if copyflags:
				f1 = open(dirpath + "/flags",'w')
				f1.write(templateflags %(i))
				f1.close()
				f2 = open(dirpath + "/pdb_dump_flags",'w')
                                f2.write(templatepdbflags %(i))
                                f2.close()			
			if dump:
				writefilewithcoordinates(dat[0],dat[1],dat[2], templateDump, dirpath + "/dump.xml")
			writefilewithcoordinates(dat[0],dat[1],dat[2], templateGlobal, dirpath + "/" +outfilename + ".xml")
			if rundump:
				cmdrundump = 'cd %s; %s -parser:protocol dump.xml @pdb_dump_flags ' %(dirpath,RosettaBinary)
				print cmdrundump
				os.system("cp NGA_ideal.params %s/." %dirpath)
				os.system("cp T12_GalNAc.pdb %s/." %dirpath)	
				os.system(cmdrundump)
				if os.path.exists('startfrom.pdb'):
					print "file created %03d" %i
                                	#os.system("mv startfrom.pdb %s/." %dirpath)	

parser = argparse.ArgumentParser()
parser.add_argument('--file', help='type of rmsd calculation, options are rmsd_ca')
parser.add_argument('--output', help='list of pdb files for rmsd calculation', default='list.txt')
mkdir=True
copyflags=True
dump=True
rundump=True
RosettaBinary="/home/saipooja/RosettaDevelopmentGlycosylation/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease"
args = parser.parse_args()

print args	
inputfile=args.file
basename=args.output
writefileswithcoordinatesfromfile(inputfile,basename)
