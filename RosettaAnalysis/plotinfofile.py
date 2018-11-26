pltinfo = dict()
pltinfo["hold"] = True
pltinfo["show"] = False
pltinfo["color"] = "xkcd:sky blue"
pltinfo["marker"] = "."
pltinfo["scale"]=50
pltinfo["edgecolor"]='none'
pltinfo["outfile"]="hOGT_Native_Peptide_Protocols_HRs"
pltinfo["xfield"] = "substrate_rmsd"
pltinfo["yfield"] = "interaction_energy"
pltinfo["ylow"] = -35
#pltinfo["xlow"]=0.6
pltinfo["yhigh"]=-10
#pltinfo["xhigh"]=1.8
pltinfo["alpha"]=0.7
pltinfo["xlabel"]="Substrate Rmsd (Angstrom)"
pltinfo["ylabel"]="Interaction Energy (REU)"
pltinfo["title"]="hOGT native peptide substrate benchmark"
pltinfo["files"]=["HighResolution_NoRandomize/score.sc","HighResolution_WithBackbone/score.sc"]#,"LowResolution/score.sc"]#,"Full/score.sc","LowResolution_WithBackbone/LowResolution_WithBackbonescore.sc"]
pltinfo["keys"]=["HR(SC Refinement)","HR(BB+SC Refinement)"]#, "LR + HR(SC Refinement)"]#,"Randomize + LR + HR(BB+SC)","LR + HR(BB+SC Refinement)"]
#pltinfo["files"]=["LowResolution_WithBackbone/LowResolution_WithBackbonescore.sc"]
#pltinfo["keys"]=["LR + HR(BB+SC Refinement)"]
#pltinfo["files"]=["HighResolution_NoRandomize/score.sc","Full/score.sc"]
#pltinfo["keys"]=["HR(SC Refinement)","Randomize + LR + HR(BB+SC)"]
pltinfo["colors"]=["r","b","g","c","y"]