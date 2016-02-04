import subprocess, re, glob, os
import numpy as np
import Tkinter
from rosetta import *
rosetta.init()
from toolbox import cleanATOM
from toolbox import mutate_residue
from rosetta.protocols.analysis import *
from rosetta.protocols.vip import *


import app.pyrosetta_toolkit.modules.tools.loops as loop_tools
from app.pyrosetta_toolkit.window_main import global_variables



#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' } 

#not allowing mutation to cys            
rev_resdict = { 'A':'ALA', 'D':'ASP', 'E':'GLU', 'F':'PHE', \
                'G':'GLY', 'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', \
                'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN', 'R':'ARG', \
                'S':'SER', 'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'} 

def makeFoldxBuildModel(mutant, name, output_pdb):
  output = open('individual_list.txt', 'w')
  output.write(mutant + ';')
  output.close()
  
  to_file = '<TITLE>FOLDX_runscript;\n'+\
            '<JOBSTART>#;\n'+\
            '<PDBS>#;\n'+\
            '<BATCH>list.txt;\n'+\
            '<COMMANDS>FOLDX_commandfile;\n'+\
            '<BuildModel>test.out,individual_list.txt;\n'+\
            '<END>#;\n'+\
            '<OPTIONS>FOLDX_optionfile;\n'+\
            '<Temperature>298;\n'+\
            '<R>#;\n'+\
            '<pH>7;\n'+\
            '<IonStrength>0.050;\n'+\
            '<water>-CRYSTAL;\n'+\
            '<metal>-CRYSTAL;\n'+\
            '<VdWDesign>2;\n'+\
            '<OutPDB>' + output_pdb + ';\n'+\
            '<pdb_hydrogens>;\n'+\
            '<END>#;\n'+\
            '<JOBEND>#;\n'+\
            '<ENDFILE>#;\n'

  output = open(name, 'w')
  output.write(to_file)
  output.close()

def makeFoldxPositionScan(mutant, name, output_pdb):
  output = open('individual_list.txt', 'w')
  output.write(mutant + ';')
  output.close()

  to_file = '<TITLE>FOLDX_runscript;\n'+\
            '<JOBSTART>#;\n'+\
            '<PDBS>#;\n'+\
            '<BATCH>list.txt;\n'+\
            '<COMMANDS>FOLDX_commandfile;\n'+\
            '<PositionScan>#,' + mutant + ';\n'+\
            '<END>#;\n'+\
            '<OPTIONS>FOLDX_optionfile;\n'+\
            '<Temperature>298;\n'+\
            '<R>#;\n'+\
            '<pH>7;\n'+\
            '<IonStrength>0.050;\n'+\
            '<water>-CRYSTAL;\n'+\
            '<metal>-CRYSTAL;\n'+\
            '<VdWDesign>2;\n'+\
            '<OutPDB>' + output_pdb + ';\n'+\
            '<pdb_hydrogens>false;\n'+\
            '<END>#;\n'+\
            '<JOBEND>#;\n'+\
            '<ENDFILE>#;\n'

  output = open(name, 'w')
  output.write(to_file)
  output.close()
  
def runFoldxSimpleMutator(mutant, pdbs):
  output = open('list.txt', 'w')
  fns = ''
  for pdb in pdbs:
    fns += pdb + '\n'
  output.write(fns)
  output.close()

  m=mutant[2:-1]
  here=mutant[-1]

  name=rev_resdict[mutant[-1]]+mutant[2:-1]+'_'+pdbs
  print("name of file",str(name))
  
  initial_pose=pose_from_pdb(pdbs)
  res = initial_pose.residue(int(mutant[2:-1]))
  #Set up ScoreFunction.
  sf = get_fa_scorefxn()

  pose2=mutate_residue( initial_pose, int(mutant[2:-1]),mutant[-1],10,sf )
  print("mutant", mutant)
  print("pdbs", pdbs)

  pose2.dump_pdb(str(name))

def makeFoldxStability(name):
  to_file = '<TITLE>FOLDX_runscript;\n'+\
            '<JOBSTART>#;\n'+\
            '<PDBS>#;\n'+\
            '<BATCH>list.txt;\n'+\
            '<COMMANDS>FOLDX_commandfile;\n'+\
            '<Stability>Stability.txt;\n'+\
            '<END>#;\n'+\
            '<OPTIONS>FOLDX_optionfile;\n'+\
            '<Temperature>298;\n'+\
            '<R>#;\n'+\
            '<pH>7;\n'+\
            '<IonStrength>0.050;\n'+\
            '<water>-CRYSTAL;\n'+\
            '<metal>-CRYSTAL;\n'+\
            '<VdWDesign>2;\n'+\
            '<OutPDB>false;\n'+\
            '<pdb_hydrogens>;\n'+\
            '<END>#;\n'+\
            '<JOBEND>#;\n'+\
            '<ENDFILE>#;\n'

  output = open(name, 'w')
  output.write(to_file)
  output.close()

def runFoldxStability(name, pdbs):
  output = open('list.txt', 'w')
  fns = ''
  for pdb in pdbs:
    fns += pdb + '\n'
  output.write(fns)
  output.close()
  
  name = 'run_' + name + '.foldx'
  makeFoldxStability(name)
  command = '/home/agm854/local/bin/foldx3b6 -runfile ' + name
  status = -1
  while status < 0:
    status = subprocess.call(command, shell=True)
  print('\n\nStability exit status = ' + str(status))
  
def makeFoldxAnalyzeComplex(name):
  to_file = '<TITLE>FOLDX_runscript;\n'+\
            '<JOBSTART>#;\n'+\
            '<PDBS>#;\n'+\
            '<BATCH>list.txt;\n'+\
            '<COMMANDS>FOLDX_commandfile;\n'+\
            '<AnalyseComplex>#;\n'+\
            '<END>#;\n'+\
            '<OPTIONS>FOLDX_optionfile;\n'+\
            '<Temperature>298;\n'+\
            '<R>#;\n'+\
            '<pH>7;\n'+\
            '<IonStrength>0.050;\n'+\
            '<water>-CRYSTAL;\n'+\
            '<metal>-CRYSTAL;\n'+\
            '<VdWDesign>2;\n'+\
            '<OutPDB>false;\n'+\
            '<pdb_hydrogens>;\n'+\
            '<END>#;\n'+\
            '<JOBEND>#;\n'+\
            '<ENDFILE>#;\n'

  output = open(name, 'w')
  output.write(to_file)
  output.close()

def getscore(pose):
  sf = get_fa_scorefxn()
  post_pre_packing_score = sf(pose)

  return(post_pre_packing_score)

def minimize(pose):
  sf = get_fa_scorefxn()

   # Set up MoveMap.
  mm = MoveMap()
  mm.set_bb(True)

  minmover = MinMover(mm, sf, 'dfpmin_armijo_nonmonotone', 10, True) ## I don't know the meaning of dfpmin,10,True. I saw it somewhere and used it

  print("finished pack")
  minmover.apply(pose)
 
  print("finish moved")
  return(pose)	


def runFoldxAnalyzeComplex(name, pdbs):
  print(pdbs)
  print("name:", name)
  output = open('list.txt', 'w')
  fns = ''
  for pdb in pdbs:
    fns += pdb + '\n'
  output.write(fns)
  output.close()

  #if we dont care about the number of chains
 
  wt_pose=pose_from_pdb(pdbs[0])
  mut_pose=pose_from_pdb(pdbs[1])

  chain_wt=wt_pose.split_by_chain()

  wt_a=getscore(chain_wt[1]) 
  wt_b=getscore(chain_wt[2])  

  
  print("seq",mut_pose.sequence())
  print("info", mut_pose.pdb_info())

  chain_mut=mut_pose.split_by_chain()

  mut_a=getscore(chain_mut[1]) 
  print("finsiehd mut a")
  mut_b=getscore(chain_mut[2])
  print("finished mut b") 
  print("done getting indivd scores")

  sf = get_fa_scorefxn()
  #analysis.analyze_interface(wt_pose,sf)
  interface_mover = InterfaceAnalyzerMover(1, False, sf, False, False, False, False )

  #wt_pose=minimize(wt_pose)

  #print("done min complex")
  interface_mover.apply(wt_pose)
  #print("done infer mover apply") 

  wt_interface=interface_mover.get_interface_dG()


  print("wt done")
 
  #mut_pose=minimize(mut_pose)
  interface_mover_mut = InterfaceAnalyzerMover(1, False, sf, False, False, False, False )
  interface_mover_mut.apply(mut_pose)

  mut_interface=interface_mover_mut.get_interface_dG()

  #chains_mut=mut_pose.split_by_chain()
  

  print("mut done")
  

  indiv_energy1 = 'Indiv_energies_AnalyseComplex_' + name[:-8] + '.fxout'
  indiv_energy2 = 'Indiv_energies_AnalyseComplex_' + name[:-8] + '.wt.fxout'

  #wild type first
  org_file=open(indiv_energy2,'w')
  line0="Pdb\tchian\tenergy\n"
  line1=name + '\t'+'A'+'\t'+ str(wt_a)+'\n'
  line2=name + '\t'+'B'+'\t'+ str(wt_b)+'\n'
  org_file.write(line0)
  org_file.write(line1)
  org_file.write(line2)
  org_file.close()

  mut_file=open(indiv_energy1,'w')
  line1=name+'\t'+'A'+'\t'+str(mut_a)+'\n'
  line2=name+'\t'+'B'+'\t'+str(mut_b)+'\n'
  mut_file.write(line0)
  mut_file.write(line1)
  mut_file.write(line2)
  mut_file.close()
 
  inter_energy1 = 'Interaction_AnalyseComplex_' + name[:-8] + '.fxout'
  inter_energy2 = 'Interaction_AnalyseComplex_' + name[:-8] + '.wt.fxout'

  mut_inter=open(inter_energy1,"w")
  line1=name+'\t'+'A'+'\t'+'B'+'\t'+str(mut_interface)+'\n'
  mut_inter.write(line1)
  mut_inter.close()

  org_inter=open(inter_energy2,"w")
  line1=name+'\t'+'A'+'\t'+'B'+'\t'+str(wt_interface)+'\n'
  org_inter.write(line1)
  org_inter.close()
  

def runFoldxAnalyzeComplex_single_chain(name, pdbs):
  print(pdbs)
  print("name:", name)
  output = open('list.txt', 'w')
  fns = ''
  for pdb in pdbs:
    fns += pdb + '\n'
  output.write(fns)
  output.close()

  #if we dont care about the number of chains
 
  wt_pose=pose_from_pdb(pdbs[0])
  mut_pose=pose_from_pdb(pdbs[1])

  wt_score=getscore(wt_pose) 
  mut_score=getscore(mut_pose)

  indiv_energy1 = 'Indiv_energies_AnalyseComplex_' + name[:-8] + '.fxout'
  indiv_energy2 = 'Indiv_energies_AnalyseComplex_' + name[:-8] + '.wt.fxout'

  #wild type first
  org_file=open(indiv_energy2,'w')
  line0="Pdb\tchian\tenergy\n"
  line1=name + '\t'+'A'+'\t'+ str(wt_score)+'\n'
  #line2=name + '\t'+'B'+'\t'+ str(wt_b)+'\n'
  org_file.write(line0)
  org_file.write(line1)
  #org_file.write(line2)
  org_file.close()

  mut_file=open(indiv_energy1,'w')
  line1=name+'\t'+'A'+'\t'+str(mut_score)+'\n'
  #line2=name+'\t'+'B'+'\t'+str(mut_b)+'\n'
  mut_file.write(line0)
  mut_file.write(line1)
  #mut_file.write(line2)
  mut_file.close()
 
  inter_energy1 = 'Interaction_AnalyseComplex_' + name[:-8] + '.fxout'
  inter_energy2 = 'Interaction_AnalyseComplex_' + name[:-8] + '.wt.fxout'

  mut_inter=open(inter_energy1,"w")
  line1=name+'\t'+'A'+'\t'+'B'+'\t'+str(0.0)+'\n'
  mut_inter.write(line1)
  mut_inter.close()

  org_inter=open(inter_energy2,"w")
  line1=name+'\t'+'A'+'\t'+'B'+'\t'+str(0.0)+'\n'
  org_inter.write(line1)
  org_inter.close()

def runFoldxAnalyzeComplex_single(name, pdbs):
  print(pdbs)
  print("name:", name)
  output = open('list.txt', 'w')
  fns = ''
  for pdb in pdbs:
    fns += pdb + '\n'
  output.write(fns)
  output.close()
  
  print("after run of files")
  wt_pose=pose_from_pdb(pdbs[0])
  print("got pose")


  print("calling getscore")
  chain_wt=wt_pose.split_by_chain()
  print("split 1 worked")
  print(chain_wt[1].pdb_info().name())
  print("got 1")
  print(chain_wt[2].pdb_info().name())
  print("got 2")
  wt_a=getscore(chain_wt[1]) 
  print("finished wt a")
  print("split 2 worked")
  wt_b=getscore(chain_wt[2])  
  print("finished wt b")
  
  sf = get_fa_scorefxn()
  interface_mover = InterfaceAnalyzerMover(1, False, sf, False, False, False, False )
  interface_mover.apply(wt_pose)
  wt_interface=interface_mover.get_interface_dG()

  indiv_energy2 = 'Indiv_energies_AnalyseComplex_' + name[:-8] + '.wt.fxout'

  #wild type first
  org_file=open(indiv_energy2,'w')
  line0="Pdb\tchian\tenergy\n"
  line1=name + '\t'+'A'+'\t'+ str(wt_a)+'\n'
  line2=name + '\t'+'B'+'\t'+ str(wt_b)+'\n'
  org_file.write(line0)
  org_file.write(line1)
  org_file.write(line2)
  org_file.close()

  inter_energy2 = 'Interaction_AnalyseComplex_' + name[:-8] + '.wt.fxout'

  org_inter=open(inter_energy2,"w")
  line1=name+'\t'+'A'+'\t'+'B'+'\t'+str(wt_interface)+'\n'
  org_inter.write(line1)
  org_inter.close()
  
def makeFoldxRepair(name):

  print("cleaning:", name)
  cleanATOM(name)

def makeInitRepair(name):
  name =  name + '.pdb' 
  initial_pose = pose_from_pdb(name)

  # Set up ScoreFunction.
  sf = get_fa_scorefxn()


  # Pack
  pre_pre_packing_score = sf(initial_pose)

  task = standard_packer_task(initial_pose)
  task.restrict_to_repacking()
  task.or_include_current(True)
  pack_rotamers_mover = RotamerTrialsMover(sf, task)
  pack_rotamers_mover.apply(initial_pose)

  #minimize till it doesnt get any better
  init_min=minimize(initial_pose)
  init_score=getscore(init_min)
  min_rep=minimize(init_min)
  min_rep_score=getscore(min_rep)
  threshold=.0001

  while abs(init_score - min_rep_score) > threshold :
    init_score=min_rep_score
    min_rep=minimize(min_rep)
    min_rep_score=getscore(min_rep)
  
  min_rep.dump_pdb(str(name))

def runFoldxRepair(name):
 
  name =  name + '.pdb'    
  makeFoldxRepair(name)
 
  #print('\n\nRepair exit status = ')


def runFoldxRepair_list(name, pdbs):
  output = open('list.txt', 'w')
  fns = ''
  for pdb in pdbs:
    fns += pdb + '\n'
    makeFoldxRepair(pdb)
  output.write(fns)
  output.close()
  



def checkOutputMutator(prefix):
  files = glob.glob('*_' + prefix + '.pdb')
  if len(files) == 2:
    return(True)
  else:
    print('\n\nWrong number of output mutator files.\n')
    return(False)

def checkOutputAnalyzeComplex(prefix):
  indiv_energy1 = 'Indiv_energies_AnalyseComplex_' + prefix + '.fxout'
  indiv_energy2 = 'Indiv_energies_AnalyseComplex_' + prefix + '.wt.fxout'
  inter_energy1 = 'Interaction_AnalyseComplex_' + prefix + '.fxout'
  inter_energy2 = 'Interaction_AnalyseComplex_' + prefix + '.wt.fxout'

  if os.path.isfile(indiv_energy1) and os.path.isfile(indiv_energy2) and os.path.isfile(inter_energy1) and os.path.isfile(inter_energy2):
    print('\n\nThese files are present: ' + indiv_energy1 + ' ' + indiv_energy2 + ' ' + inter_energy1 + ' ' + inter_energy2 + '\n')
    return(True)
  else:
    print('\n\nSome files are missing from Analyze Complex\n')
    return(False)

class Scores:
  def __init__(self):    
    self.files = []
    self.total_energies = []
    self.ids = []
    self.wts = []
    self.best_mutant = []
    self.group1 = []
    self.group2 = []
    self.clashes1 = []
    self.clashes2 = []
    self.stability = []
    self.interaction_energies = []
    
  def parseAnalyzeComplex(self):
    self.parseFiles('Interaction_AnalyseComplex_')
    if 'wt' not in self.files[0] and len(self.files) > 1:
      new_files = [self.files[1], self.files[0]]
      self.files = new_files
      
    for a_file in self.files:
      lines = open(a_file, 'r').readlines()
      start = True
      for a_line in lines:
        energy_list = a_line.split('\t')
        print(energy_list)
        #if not start and energy_list[0] == 'Pdb':
        #start = True
        if start and len(energy_list) > 1:
          print("setting id:",self.removeWhiteSpace(energy_list[0]))
          (self.ids).append(self.removeWhiteSpace(energy_list[0]))
          (self.group1).append(self.removeWhiteSpace(energy_list[1]))
          (self.group2).append(self.removeWhiteSpace(energy_list[2]))
          #(self.clashes1).append(self.removeWhiteSpace(energy_list[3]))
          #(self.clashes2).append(self.removeWhiteSpace(energy_list[4]))
          (self.interaction_energies).append(self.removeWhiteSpace(energy_list[3]))

    self.parseFiles('Indiv_energies_')
    if 'wt' not in self.files[0] and len(self.files) > 1:
      new_files = [self.files[1], self.files[0]]
      self.files = new_files
          
    count = 0
    for a_file in self.files:
      lines = open(a_file, 'r').readlines()
      start = False
      stability = []
      for a_line in lines:
        print(a_line)
        energy_list = a_line.split('\t')
        if not start and energy_list[0] == 'Pdb':
          start = True
        elif start and len(energy_list) > 1 and self.ids[count] == energy_list[0]:
          stability.append(self.removeWhiteSpace(energy_list[2]))
          print("stability in loop:",stability)
      (self.stability).append(stability)
      count += 1
           
  def parseStability(self):
    self.parseFiles('Stability.txt')
    for a_file in self.files:
      lines = open(a_file, 'r').readlines()
      start = False
      for a_line in lines:
        energy_list = a_line.split('\t')
        if not start and energy_list[0] == '':
          start = True
        elif start and len(energy_list) > 1:
          (self.ids).append(energy_list[0])
          (self.total_energies).append(energy_list[1])
      
  def parseFiles(self, prefix):
    import os
    self.files = []
    raw_files = raw_files = os.listdir('.')
    
    for i in raw_files:
      if i.startswith(prefix):
        (self.files).append(i)
   
  def removeWhiteSpace(self, string):
    return("".join(string.split()))
    
  def findBestMutant(self):
    best_score = np.array(self.energies).argsort()[:1][::-1]
    self.best_mutant = self.ids[best_score]

  def getTotalEnergies(self):
    return(self.total_energies)
  def getWTenergies(self):
    return(self.wts)
  def getIds(self):
    return(self.ids)
  def getBestMutant(self):
    return(self.best_mutant)
  def getGroup1(self):
    return(self.group1)
  def getGroup2(self):
    return(self.group2)
  def getClashes1(self):
    return(self.clashes1)
  def getClashes2(self):
    return(self.clashes2)
  def getInteractionEnergies(self):
    return(self.interaction_energies)
  def getStability1(self):
    print("stability is 1,", self.stability)
    return(self.stability[0])
  def getStability2(self):
    print("stability is 2,", self.stability)
    return(self.stability[1])

  def cleanUp(self, delete_files):
    import os, glob
    for name in delete_files:
      for files in glob.glob(name):
        if os.path.isfile(files):
          os.remove(files)
    for files in glob.glob('run_*.foldx'):
      if os.path.isfile(files):
        os.remove(files)
    for files in glob.glob('*out'):
      if os.path.isfile(files):
        os.remove(files)
    if os.path.isfile('missing.txt'): 
      os.remove('missing.txt')
    if os.path.isfile('Stability.txt'):
      os.remove('Stability.txt')
    if os.path.isfile('list.txt'):
      os.remove('list.txt')
    if os.path.isfile('individual_list.txt'):
      os.remove('individual_list.txt')
    if os.path.isfile('scanning_output.txt'):
      os.remove('scanning_output.txt')
