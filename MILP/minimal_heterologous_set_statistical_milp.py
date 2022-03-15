# Script to predict minimal heterologous reaction sets for value-added compound production
# Created: 23 June 2019
# Author:Szabolcs Cselgő Kovács
#
# Citation: This script is based on a knock-in OptCouple example MILP script (Jensen et al., 2019)
#       Kristian Jensen, Valentijn Broeken, Anne Sofie Lærke Hansen, Nikolaus Sonnenschein, Markus J. Herrgård,
#            OptCouple: Joint simulation of gene knockouts, insertions and medium modifications for prediction of growth-coupled strain designs,
#            Metabolic Engineering Communications, Volume 8, 2019,e00087, ISSN 2214-0301, https://doi.org/10.1016/j.mec.2019.e00087.
#
#
# Contact: kovszasz@gmail.com
# Loading dependencies
import os,pickle,sys
from multiprocessing import Pool
os.environ["OPTLANG_USE_SYMENGINE"] = "False"
import matplotlib
matplotlib.use('agg')
#import cameo
import cobra
from cobra.test import create_test_model
from optlang.duality import convert_linear_problem_to_dual
import logging
import optlang
#from cameo.models import universal
#from cameo.data import metanetx
logger=logging.getLogger()
model = cobra.io.read_sbml_model("../heterologous.xml")
model.reactions.ATPM.lower_bound=0
model.reactions.ATPM.upper_bound=0
model.reactions.Ec_biomass_iJO1366_core_53p95M.lower_bound=0
model.reactions.EX_glc_lp_e_rp_.lower_bound=-10
model.solver='cplex'

def binary_loader(name):#We store the results and other important variables binary. It is easer and faster to read and save.
    f=''
    infile=open(f+name,'rb')
    return pickle.load(infile)


start_seed=sys.argv[1]
end_seed=sys.argv[2]
name=sys.argv[3]
#path_of_random_set = sys.argv[4]
heterologous_reactions=binary_loader('../statistical_analysis/MetaCycAll')
random_reaction_set=binary_loader('../statistical_analysis/RandomHeterologousReactionSet')

def c_atom_number(s):
    e=model.reactions.get_by_id(s).metabolites
    for i in e.keys():
        a=model.metabolites.get_by_id(i.id).formula
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

def setup_milp(model, target, remove_blocked=False, exclude_reaction_ids=set()):
    """
    This function constructs the MILP.
    exclude_reaction_ids takes a list of reaction ids that shouldn't be considered for heterologous addition
    (i.e. spontaneous reactions and exchange reactions). These reactions are thus always allowed to have flux within their model bounds.
    """
    original_model = model
    model = model.copy()
    model.objective=target
    for i in heterologous_reactions.keys(): # turns off each heterologous reaction in order to get the only the native metabolic network.
        model.reactions.get_by_id(i).lower_bound=0
        model.reactions.get_by_id(i).upper_bound=0

    for r in model.reactions:
        if r.id.find('MetaCyc')>-1:
            r.type='heterologous'
        else:
            r.type='native'

    # Set the solver to Gurobi for the fastest result. Set to CPLEX if Gurobi is not available.
    if "gurobi" in cobra.util.solver.solvers.keys():
        logger.info("Changing solver to Gurobi and tweaking some parameters.")
        if "gurobi_interface" not in model.solver.interface.__name__:
            model.solver = "gurobi"
        # The tolerances are set to the minimum value. This gives maximum precision.
        problem = model.solver.problem
        problem.params.NodeMethod = 1 # primal simplex node relaxation
        problem.params.FeasibilityTol = 1e-9 #If a flux limited to 0 by a constraint, which range around it is still considered the same as 0 > set smallest possible
        problem.params.OptimalityTol = 1e-3 #how sure the solver has to be about this optimum being really the best it has.
        problem.params.IntFeasTol = 1e-9 #If a value is set to an integer, how much may it still vary? > set smallest possible
        problem.params.MIPgapAbs = 1e-9
        problem.params.MIPgap = 1e-9
        problem.params.Threads=1 #In order to reduce memory usage (increased calculation time)
        problem.params.TimeLimit = 200 # Use max 200 seconds when called, return best solution after that
        problem.params.PoolSearchMode = 1 #0 for only finding the optimum, 1 for finding more solutions (but no quality guaranteed), 2 for finding the n best possible solutions
        problem.params.PoolSolutions = 10 # Number of solutions kept when finding the optimal solution
        problem.params.PoolGap = 0.9 # only store solutions within 90% of the optimal objective value

    elif "cplex" in cobra.util.solver.solvers.keys():
        logger.warning("Changing solver to CPLEX, as Gurobi is not available. This may cause a big slowdown and limit options afterwards.")
        if "cplex_interface" not in model.solver.interface.__name__:
            model.solver = "cplex"
        # The tolerances are set to the minimum value. This gives maximum precision.
        problem = model.solver.problem
        problem.parameters.mip.strategy.startalgorithm.set(1) # primal simplex node relaxation
        problem.parameters.simplex.tolerances.feasibility.set(1e-9) #If a flux limited to 0 by a constraint, which range around it is still considered the same as 0 > set smallest possible
        problem.parameters.simplex.tolerances.optimality.set(1e-3) #possibly fine with 1e-3, try if allowed. Is how sure the solver has to be about this optimum being really the best it has.
        problem.parameters.mip.tolerances.integrality.set(1e-9) #If a value is set to an integer, how much may it still vary? > set smallest possible
        problem.parameters.mip.tolerances.absmipgap.set(1e-9)
        problem.parameters.mip.tolerances.mipgap.set(1e-9)
        problem.parameters.mip.pool.relgap.set(0.9) # For populate: find all solutions within 10% of the optimum for relgap = 0.1
        problem.parameters.timelimit.set(200) # Use max 200 seconds for solving
        problem.parameters.mip.limits.populate.set(20) # Find max 20 solutions (=default)

    else:
        logger.warning("You are trying to run 'Hamlet Hot Rig' with %s. This might not end well." %
                      model.solver.interface.__name__.split(".")[-1])
        pass
    # Remove reactions that are blocked: no flux through these reactions possible. This will reduce the search space for the solver, if not done already.
    if remove_blocked:
        blocked_reactions = cameo.flux_analysis.analysis.find_blocked_reactions(model)
        model.remove_reactions(blocked_reactions)

    # Make dual
    model_with = model.copy() # This variable looks unnecessary, but is kept out of fear of messing stuff up
    model_with.optimize()
    dual_problem = convert_linear_problem_to_dual(model_with.solver)
    logger.debug("Dual problem successfully created")

    # Combine primal and dual
    primal_problem = model.solver

    for var in dual_problem.variables:  # All variables in the dual are copied to the primal
        var = primal_problem.interface.Variable.clone(var)
        primal_problem.add(var)
    for const in dual_problem.constraints:  # All constraints in the dual are copied to the primal
        const = primal_problem.interface.Constraint.clone(const, model=primal_problem)
        primal_problem.add(const)
    logger.debug("Dual and primal combined")

    dual_problem.optimize()

    # Dictionaries to hold the binary control variables:
    heterologous_y_vars = {} # 1 for 'knockin', 0 for inactive
    medium_y_vars = {}       # 1 for medium addition (up to -10), 0 for no addition

    # Now the fun stuff
    constrained_dual_vars = set()

    # For the knockins and medium additions:
    for reaction in [r for r in model.reactions if r.type == "heterologous"]:
        # Add constraint variables
        interface = model.solver.interface
        y_var = interface.Variable("y_" + reaction.id, type="binary")

        # Constrain the primal: flux through reactions maximum within (-1000, 1000), or smaller boundaries defined before
        model.solver.add(interface.Constraint(reaction.flux_expression - 1000 * y_var, ub=0, name="primal_y_const_"+reaction.id+"_ub"))
        model.solver.add(interface.Constraint(reaction.flux_expression + 1000 * y_var, lb=0, name="primal_y_const_"+reaction.id+"_lb"))

        # Constrain the dual
        constrained_vars = []

        if reaction.upper_bound != 0:
            dual_forward_ub = model.solver.variables["dual_" + reaction.forward_variable.name + "_ub"]
            model.solver.add(interface.Constraint(dual_forward_ub - 1000 * (1 - y_var), ub=0))
            constrained_vars.append(dual_forward_ub)
        if reaction.lower_bound != 0:
            dual_reverse_ub = model.solver.variables["dual_" + reaction.reverse_variable.name + "_ub"]
            model.solver.add(interface.Constraint(dual_reverse_ub - 1000 * (1 - y_var), ub=0))
            constrained_vars.append(dual_reverse_ub)
        constrained_dual_vars.update(constrained_vars)

        # Add y variable to the corresponding modifications dictionary
        heterologous_y_vars[y_var] = reaction

    logger.debug("Control variables created")

    # Add number of heterologous switch contraint constraint
    heterologous_turn_on = model.solver.interface.Constraint(
                            optlang.symbolics.Add(*heterologous_y_vars), lb=0, ub=0, name="heterologous_reaction_constraint"
                        )
    model.solver.add(heterologous_turn_on)

    # Set the objective
    primal_objective = model.solver.objective
    dual_objective = interface.Objective.clone(
        dual_problem.objective, model=model.solver
    )

    switch_objective=interface.Objective(heterologous_turn_on.expression, direction='min')

    full_objective = interface.Objective(primal_objective.expression-dual_objective.expression, direction="max")
    model.objective = full_objective

    return model,primal_objective,dual_objective,full_objective,switch_objective


def run(target_metabol,seed_list):
    target_id = "EX_"+target_metabol
    try:
        target=model.reactions.get_by_id(target_id)
    except:
        # Add the target reaction and a needed exporter for the metabolite created.
        EX_ppa_c = cobra.Reaction(target_id, name=target_id)
        model.add_reactions([EX_ppa_c])
        model.reactions.get_by_id(EX_ppa_c.id).reaction = target_metabol+" -->"
        model.reactions.get_by_id(EX_ppa_c.id).gene_reaction_rule = "targetgene"

    biomass_id = "Ec_biomass_iJO1366_core_53p95M"
    milp,primal_objective,dual_objective,full_objective,switch_objective=setup_milp(model,target_id,remove_blocked=False)
    milp.reactions.Ec_biomass_iJO1366_core_53p95M.lower_bound=0.09915628719962316
    milp.reactions.ATPM.upper_bound=1000
    milp.reactions.ATPM.lower_bound=3.15
    milp.objective=primal_objective
    milp.optimize()
    partial_result=[]
    new_yield_demand=milp.reactions.get_by_id(target_id).flux*1.05
    if abs(new_yield_demand)<1e-4:
        #new_yield_demand=0.94 ---old solution 2020.10.28
        new_yield_demand=0.05
    print('New yield demand:\t',new_yield_demand)

    for seed in seed_list:
        print('Seed:\t',seed)
        partial_milp=milp.copy()
        for i in heterologous_reactions.keys(): # turns off each heterologous reaction in order to get the only the native metabolic network.
            partial_milp.reactions.get_by_id(i).lower_bound=0
            partial_milp.reactions.get_by_id(i).upper_bound=0
        for i in random_reaction_set[seed]: # turns on each desired heterologous reaction, which is included in the turned_on_heterologous_reactions list
            if heterologous_reactions[i]:
                partial_milp.reactions.get_by_id(i).lower_bound=-1000
                partial_milp.reactions.get_by_id(i).upper_bound=1000
            else:
                partial_milp.reactions.get_by_id(i).lower_bound=0
                partial_milp.reactions.get_by_id(i).upper_bound=1000

        partial_milp.solver.constraints.heterologous_reaction_constraint.ub=563
        partial_milp.objective=switch_objective
        partial_milp.reactions.get_by_id(target_id).lower_bound=new_yield_demand
        try:
            partial_milp.optimize()
            #partial_result.append({partial_milp.solver.objective.value:len([ i.id for i in partial_milp.reactions if i.id.find('MetaCyc')>-1 and abs(i.flux)>1e-4])})
            partial_result.append(round(partial_milp.solver.objective.value))
        except:
            print('Infeasible at '+target_metabol+' at seed '+str(seed))
        del partial_milp
    #print(target_id,'\t',milp.reactions.get_by_id(target_id).flux)
    return partial_result

def call_thread(payload):
    #f=open('partial_result'+str(start_seed)+'_'+str(end_seed)+'.txt','w')
    target=payload['target']
    start_seed=payload['start_seed']
    end_seed=payload['end_seed']
    #for target in targets:
    partial_result=[]
    seed_list=[i for i in range(int(start_seed),int(end_seed))] 
    partial_result=run(target,seed_list)
    return partial_result

def binary_saver(obj,name):#We store the results and other important variables binary. It is easer and faster to read and save.
    f=''
    outfile=open(f+name,'wb')
    pickle.dump(obj,outfile)
    outfile.close()


def white_space_remover(string):#while reading the '.txt' files, the lines are separated wit the '\n' character. This function removes '\n' from the end of the strings.
    return string.replace('\n','')

targets=[]
with open("value_added_compound_list_for_milp.txt") as target_compounds:
    for i in target_compounds:
        targets.append(white_space_remover(i))
results={}
#targets=['BUTANOL_c']
p = Pool(8)
csv = open('results/result'+str(name)+'.csv','w')
counter=0
all_=len(targets)
for target in targets:
    print(target,'\t',str((counter/all_)*100)+"%")
    #results[target]=call_thread({'target':target,'start_seed':20000,'end_seed':20002})
    pr=call_thread({'target':target,'start_seed':start_seed,'end_seed':end_seed})
    line=target
    for l in pr:
        line=line+','+str(l)
    csv.write(line+'\n')
    results[target]=pr
    #results[target] = p.map(call_thread,[
        #{'target':target,'start_seed':20000,'end_seed':20002},
        #{'target':target,'start_seed':20000,'end_seed':20125},
        #{'target':target,'start_seed':20125,'end_seed':20250},
        #{'target':target,'start_seed':20250,'end_seed':20375},
        #{'target':target,'start_seed':20375,'end_seed':20500},
        #{'target':target,'start_seed':20500,'end_seed':20625},
        #{'target':target,'start_seed':20625,'end_seed':20750},
        #{'target':target,'start_seed':20750,'end_seed':20875},
        #{'target':target,'start_seed':20875,'end_seed':21000}
       # ])
    counter=counter+1
print(results)
csv.close()
binary_saver(results,'results/MilpHetStat'+str(name))
