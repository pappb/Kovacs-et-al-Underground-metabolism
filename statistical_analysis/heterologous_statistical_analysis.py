# Script to statistically asess the industrial potential of heterologous reactions
# Author: Szabolcs Cselgo Kovacs
# Contact: kovszasz@gmail.com


import sys,cobra, cobra.io, cobra.test, csv,pandas,codecs, email, smtplib,sys,random,pickle
#from cameo import *
#from cameo import load_model
from cobra import Reaction, Metabolite, Model
import cobra, cobra.io, cobra.test, csv,pandas,codecs,pickle
from cobra.util.solver import linear_reaction_coefficients
from random import randint
#import cameo
import numpy as np
#model=load_model("heterologous.xml")
model="../heterologous.xml"
model=cobra.io.read_sbml_model(model)
model.solver='gurobi'

def white_space_remover(string):#while reading the '.txt' files, the lines are separated wit the '\n' character. This function removes '\n' from the end of the strings.
    return string.replace('\n','')


def binary_saver(obj,name):#We store the results and other important variables binary. It is easer and faster to read and save.
    f=''
    outfile=open(f+name,'wb')
    pickle.dump(obj,outfile)
    outfile.close()

def binary_loader(name):#We store the results and other important variables binary. It is easer and faster to read and save.
    f=''
    infile=open(f+name,'rb')
    return pickle.load(infile)

heterologous_reactions={}#the key is the ID of the heterologous reaction, the item is information about its reversibility.
for i in model.reactions:
    if i.id.find('Meta')>-1:
        heterologous_reactions[i.id]=i.reversibility
        i.knock_out() #all heterologous reaction will be knocked out, because the randomly generated reaction sets are loaded by the path_random input

path_random='RandomHeterologousReactionSet' #input("Reaction set:") #a path for the randomly generated heterologous reaction sets
random_reaction_set=binary_loader(str(path_random)) #loading the randomly generated heterologous reaction sets
class ProductionPrediction():
    def __init__(self,model,target,C_source,random_set,start_seed,end_seed,heterologous_reactions):
        """
        model: heterologous metabolic network
        target: the ID of the target compound
        random_set: the randomly generated heterologous reaction set
        start_seed: the random seeds are between 8011-9011 (1000 randomly generated reaction set)
        end_seed: the random seeds are between 8011-9011 (1000 randomly generated reaction set)
        heterologous_reactions: dict of the heterologous reactions, containing information about their reversibility
        """
        self.model=model
        self.model.reactions.get_by_id("Ec_biomass_iJO1366_core_53p95M").lower_bound=0.09915628719962316 #10% of the maximal biomass production value of native metabolic network. This value based on the compound_production.py results
        self.model.reactions.get_by_id(C_source).lower_bound=-10
        self.model.reactions.get_by_id(C_source).upper_bound=1000
        self.heterologous_reactions=heterologous_reactions
        self.target=target
        self.C_source=C_source
        self.seed_list=[i for i in range(int(start_seed),int(end_seed))] #generating the random seed range
        self.reaction_set=random_set #a dictionary, where the key is a random seed number, and the item is a list containing the randomly generated heterologous reactions
        self.result=[]
        self.create_target_demand_reaction(self.target)
        for i in self.seed_list: #iterating the among the random seed values
            self.result.append(self.run_prediction(self.target,self.C_source,self.reaction_set[i]))
        self.model.remove_reactions('DEMAND_target')

    def c_atom_number(self,a):#For the yield calculation it is required to get information about the carbon atom numbe of the given compund. This function calculates it from the formula of the compound.
        C={}
        for i in ['Cl','H','O','S','P','C','N']:
            if a.find(i)>-1:
                a=a.replace(i,','+i+';')
        l1=a.split(',')
        for i in l1:
            p=i.split(';')
            if len(p)==2:
                if p[1]=='':
                    C[p[0]]=1
                else:
                    C[p[0]]=int(p[1])
        return C['C']

    def create_target_demand_reaction(self,target):
        target_metabol=self.model.metabolites.get_by_id(target)
        target_demand=Reaction(id="DEMAND_target")
        target_demand.add_metabolites({target_metabol:-1.0})
        self.model.add_reaction(target_demand)
        self.model.reactions.get_by_id("DEMAND_target").lower_bound=0
        self.model.reactions.get_by_id("DEMAND_target").upper_bound=1000
        self.model.objective={self.model.reactions.DEMAND_target:1}

    def run_prediction(self,target,C_source,heterologous_reaction_list):
        """
        This function simulates the production of the target compound, while only the randomly generated heterologous reaction set is turned on
        """
        for i in self.heterologous_reactions.keys(): # if any undesired heterologous reaction might be turned on (undesired: not part of the randomly generated heterologous reaction list)
            self.model.reactions.get_by_id(i).knock_out()
        for i in heterologous_reaction_list:
            if self.heterologous_reactions[i]:
                self.model.reactions.get_by_id(i).lower_bound=-1000
                self.model.reactions.get_by_id(i).upper_bound=1000
            else:
                self.model.reactions.get_by_id(i).lower_bound=0
                self.model.reactions.get_by_id(i).upper_bound=1000
        try:
            p=cobra.flux_analysis.pfba(self.model)

            print(p.fluxes['DEMAND_target'])
            compound_yield=abs((p.fluxes["DEMAND_target"]/p.fluxes[C_source])*(self.c_atom_number(self.model.metabolites.get_by_id(target).formula)/6)) # 6 represents the carbon atom number of glucose (default carbon soruce)
        except:
            compound_yield=0
        return compound_yield
start_seed=sys.argv[1]
end_seed=sys.argv[2]
name=sys.argv[3]


targets=[]
with open("value_added_compound_list_for_heterologous_model.txt") as target_compounds:
    for i in target_compounds:
        targets.append(white_space_remover(i))

result={}
#targets=['ptrc_c']
carbon_source="EX_glc_lp_e_rp_"
for target in targets:
    print(target)
    result[target]=ProductionPrediction(model,target,carbon_source,random_reaction_set,start_seed,end_seed,heterologous_reactions).result

binary_saver(result,'results/'+str(name))
