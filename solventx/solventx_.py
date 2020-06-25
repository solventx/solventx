"""
This describes the solvent extraction model for simulating REE components
separation in a multi-stage process. This contains the solvent extraction object 
definition and methods (including objective function/reward and constraints) 
that can be passed to a RL Agent or an external optimizer
- Research supported by the US Department of Energy 

This Implementation embeds the RBFOPT Black-box function Implementation adapted
for optimizing this solvent extraction design. 
RBFOPT details:
-Licensed under Revised BSD license, see LICENSE.
-(C) Copyright Singapore University of Technology and Design 2014.
-(C) Copyright International Business Machines Corporation 2016.
-Research partially supported by SUTD-MIT International Design Center.


"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

from scipy.optimize import root
import numpy as np
import cantera as ct
import pandas as pd
import rbfopt

from solventx import result_struct as rs
from solventx import config

import operator
import os



class solventx(rbfopt.RbfoptBlackBox):

    """
     Attributes
    ----------

    dimension : int
        Dimension of the problem.
        
    var_lower : 1D numpy.ndarray[float]
        Lower bounds of the decision variables.

    var_upper : 1D numpy.ndarray[float]
        Upper bounds of the decision variables.

    var_type : 1D numpy.ndarray[char]
        An array of length equal to dimension, specifying the type of
        each variable. Possible types are 'R' for real (continuous)
        variables, and 'I' for integer (discrete) variables.

    integer_vars : 1D numpy.ndarray[int]
        A list of indices of the variables that must assume integer
        values.
    
    """          
    coltypes                = ['Extraction','Scrub','Strip']
    ml2l                    = 0.001 # mililiter to liter
    scale                   = 1 # volumetric scaling from unit ml/sec

    g_p_kg                  = 1000 # grams to kg
    s_p_h                   = 3600 # seconds to hours
    target_conc             =  45 #g/L


    def __init__(self,confDict, confEnvDict, ree_mass):

                 
        """Constructor.
        """

     
        self.confDict = confDict 
        solventxHome = self.confDict["solventxHome"]
#        reeComps = self.confDict['compositions']
        self.xmlData = self.confDict['xmlData']
        self.modulesData = self.confDict['modules']      
               
        self.xml = os.path.join(solventxHome,self.xmlData['xml'],self.xmlData['phase']+'_'+''.join(self.modulesData["input"])+'.xml')
        
        # Set required data
                 
        self.phase_names    = self.confDict["phasenames"] # from xml input file
        self.phase          = ct.import_phases(self.xml,self.phase_names) 

        # Derived and/or reusable system parameters        
        self.column         = self.coltypes   # Column name
        self.solv           = self.confDict["solvents"] # solvent list -.i.e., electrolyte, extractant and organic diluent
        
        # ree by modules
        self.ree            = self.modulesData["input"] #self.REEs[self.moduleID] # (rare earth) metal list
                   
        # Cantera indices
        self.mix            = ct.Mixture(self.phase)
        self.aq             = self.mix.phase_index(self.phase_names[0])
        self.org            = self.mix.phase_index(self.phase_names[1])
        self.ns             = self.mix.n_species
        self.naq            = self.mix.phase(self.aq).n_species
        self.norg           = self.mix.phase(self.org).n_species
        
        self.HA_Index       = self.mix.species_index(self.org,'(HA)2(org)') # index of extractant in canera species list
        self.Hp_Index       = self.mix.species_index(self.aq,'H+') # index of H+ in cantera species list
        self.Cl_Index       = self.mix.species_index(self.aq,'Cl-') # index of Cl in cantera species list        
       
        
        self.canteranames   = self.mix.species_names
        self.fixed_species  = ['H2O(L)','OH-', 'Cl-', 'dodecane']
        self.canteravars    = [ij for ij in self.canteranames if ij not in self.fixed_species] # 'Cl-',

        self.nsy                     = len(self.canteravars)
        self.naqy                    = len([ij for ij in self.mix.species_names[:self.naq] if ij not in self.fixed_species])
        self.norgy                   = len([ij for ij in self.mix.species_names[self.naq:] if ij not in self.fixed_species])
        
        self.mwre,\
        self.mwslv          = self.get_mw() # g/mol
        self.rhoslv         = [1000, 960, 750] # [g/L]
                        
        self.purity_spec    = .985 # 
        self.recov_spec     = .80 # 
        
        self.Ns             = [0,0,0]

        self.nsp            = pd.DataFrame() # feed streams (aq and org) for each column       
        self.nsp0           = pd.DataFrame() # feed streams (aq and org) for each column       

        self.y              = {} # all compositions
        self.Ns             = {}
        self.constraint    = {}
        
  
        self.get_conc(ree_mass)

        self.get_process() 
        
        """ create variable space parameters and update optim attributes"""        
        var_upper, var_lower, var_type = self.create_var_space(confEnvDict)
        
        self.dimension      = len(self.var_space['mutable'])

        self.var_lower      = np.array([item for item in var_lower.values()])
        self.var_upper      = np.array([item for item in var_upper.values()])             
        self.var_type       = np.array([item for item in var_type.values()])
        
        
        # Define penalty by column type
        p0 = [5,0,10]
        pen = {}
        for item,jtem in zip(self.column, p0):
            pen[item] = jtem

        self.penalty = pen
            
        

#%%
    def get_conc(self,ree_mass):
        self.ree_mass       = ree_mass
        self.vol            = sum(self.ree_mass) / self.g_p_kg / self.target_conc   # [kg/hr]*[g/kg] /[g/L] = [L/hr] 
        self.ree_conc       = [(ij/sum(self.ree_mass)) * self.target_conc for ij in self.ree_mass] #([kg/hr]/[kg/hr] )* [g/L] = [g/L]


    # -- end function
  
#%%
    def get_mw(self, conv=ml2l):
        """ Initialize parameters for cantera simulation. init() calls this function"""
        
        mx                  = ct.Mixture(self.phase)
        aq                  = mx.phase_index(self.phase_names[0])
        org                 = mx.phase_index(self.phase_names[1])

        mwre                = np.zeros(len(self.ree)) # molecular weight of rees
        mwslv               = np.zeros(len(self.solv)) # mw & densities for 'solvents'         
        
        for re in self.ree:
            mwre[self.ree.index(re)]         = mx.phase(aq).molecular_weights[mx.phase(aq).species_index(re+'+++')]
            
        for so in self.solv:   
            if so == 'H2O(L)':
                mwslv[self.solv.index(so)]   = mx.phase(aq).molecular_weights[mx.phase(aq).species_index(so)]
            else:
                mwslv[self.solv.index(so)]   = mx.phase(org).molecular_weights[mx.phase(org).species_index(so)]

        return mwre, mwslv



#%%
    def get_process(self):
        """Get products."""
        
        
        input_components = self.confDict['modules']['input']
        strip_components = self.confDict['modules']['output']['strip']
    
    
        n_components = len(input_components)
        config_key = ''
        
        print(f'Looping through following modules config:{list(config.valid_processes.keys())}')
        for key,config_dict in config.valid_processes.items():
            
           
            if set(input_components) == set(config_dict['input']):
                if set(strip_components) ==  set(config_dict['strip']):
                    
                    config_key = key
        
        if config_key:
            print(f'Found the following process config:{config_key}')            
        else:
            raise ValueError(f'No valid configuration found for input:{input_components},strip:{strip_components}!')            
       
        modules = config.valid_processes[config_key]['modules']
        x = []
       
        print(f'Process config {config_key}:Input:{input_components},Number of modules:{len(modules)}')
        print('Modules info:')
        for key,module in modules.items():
            x.extend(module['x'])
            print(f'Module {key}:{module["strip_group"]}')
        print(f'x0:{x}')
        
        self.design_variables = x
        self.modules = modules
        self.num_input = n_components
        self.config_key = config_key
        
    


#%%
    def create_var_space(self, confEnvDict, input_feeds=1,): # 


        var_space = {
            'immutable': {},
            'mutable':   {}, #var, index in obj.variables
        }
        
        mod_space = {}
        x_space = {}
        var_upper = {}
        var_lower = {}
        var_type = {}
        
        immutable_var_names = ['mol_frac']
        feed_var_names = ['(HA)2(org)']

        index = 0
        for feed_num in range(input_feeds):
            var_space['mutable'][f'{feed_var_names[0]}-{feed_num}'] = index
            index += 1 


        for key, value in self.modules.items():
            
            mod_space[f'module-{key}'] = key
            x_space[f'module-{key}'] = value['x']

            for i, var in enumerate(value['mvn']):
                var_space['mutable'][f'{var}-{key}'] = index
                index += 1


        for comp_num in range(len(self.confDict['modules']['input'])):
            var_space['immutable'][f'{immutable_var_names[0]}-{self.confDict["modules"]["input"][comp_num]}'] = index
            index += 1

            
        self.var_space = var_space
        mutable = self.var_space['mutable']
        immutable = self.var_space['immutable']
        self.combined_var_space = combine_dict(mutable, immutable)
        
        self.mod_space = mod_space
        self.x_space = x_space
        
        var_config= confEnvDict['variable_config']
        
        for item in mutable.keys():
            name,no = item.split('-')
            var_upper[item] = var_config[name]['upper']
            var_lower[item] = var_config[name]['lower']
            var_type[item] = var_config[name]['type']
            
        return var_upper, var_lower, var_type
        
      
#%%		
    def flow_path_exist(self, var): # Not used in current implementation
        '''
            Does a proper module exist to connect to var. 
            var: Extraction-0
            returns True if Extraction-1 exists
        '''
        try:
            name, module = var.split('-')
        except ValueError: #variable doesnt exist in dictionary
            return False
        
        #get next module location
        if name == 'Extraction':
            next_ = get_next(module, 'left')
        elif name == 'Strip':
            next_ = get_next(module, 'right')
        
        index = self.combined_var_space.get(f'Extraction-{next_}')
        return (False, self.variables[index] > 0)[index != None]




#%%
    def create_nsp_open(self, name, num ):# g_p_kg=g_p_kg ): #, h0, target, xml, cantera_data, experiments):
        
        """ Create mole flows for column feed streams in a given module. Arranges
            them in the dataframe, nsp, in the form that Cantera expects """
               
        nre                 = np.zeros(len(self.ree)) # REE specie moles
        salts               = np.zeros(len(self.ree)) # neutral complexe moles
        
        strip_ree           = config.valid_processes[self.config_key]['modules'][num]['strip_group']
        is_scrub            = [1 if re in strip_ree else 0 for re in self.ree]

        
        # Determine if there is a parent column. If so, match the organic flow into current extraction stage
        # with flow from previous module strip column (which may or may not be the parent module)
        parent_col = get_parent('Extraction-'+num)
        if parent_col:
            pname, pnum = parent_col.split('-')
        else:
            pnum = ''
            
        if pnum in self.modules.keys(): # if the parent is a module in this configuration
#        if int(num) > 0: # 
            nnum = str(int(num)-1) # since modules are numbered sequentially following organic flow
            salts = np.array(self.y['Strip-'+nnum][self.canteravars.index('(HA)2(org)')+1:self.nsy]) # organic exit rees
            n_HA = np.array(self.y['Strip-'+nnum][self.canteravars.index('(HA)2(org)')]) # organic exit rees
            vol_HA = n_HA / (self.rhoslv[self.solv.index('(HA)2(org)')]/self.mwslv[self.solv.index('(HA)2(org)')])
            n_dodec = self.nsp0['Strip-'+nnum][self.canteranames.index('dodecane')]
            vol_dodec = n_dodec/(self.rhoslv[self.solv.index('dodecane')]/self.mwslv[self.solv.index('dodecane')])
            orgvol = vol_HA + vol_dodec
            

        else:            
        #        Compositions   
#            orgvol              = self.orgvol[int(num)]
            orgvol              = self.orgvol[0]
            vol_HA              = orgvol *  (self.x[self.combined_var_space['(HA)2(org)-0']])
            vol_dodec           = orgvol - vol_HA                
            n_HA                = vol_HA * self.rhoslv[self.solv.index('(HA)2(org)')]/self.mwslv[self.solv.index('(HA)2(org)')] # [L/hr]*[g/L]/[g/mol] = [mol/hr]
            n_dodec             = vol_dodec * self.rhoslv[self.solv.index('dodecane')]/self.mwslv[self.solv.index('dodecane')] # [L/hr]*[g/L]/[g/mol] = [mol/hr]            
 
        aqvols              = [orgvol/(self.x[self.combined_var_space['OA Extraction-'+num]]),orgvol/(self.x[self.combined_var_space['OA Scrub-'+num]]), orgvol/(self.x[self.combined_var_space['OA Strip-'+num]]) ]
            
        for k in range(len(self.column)): # k represents different columns - extraction, scrub or strip
    
            n_H2O           = aqvols[k] * self.rhoslv[self.solv.index('H2O(L)')]/self.mwslv[self.solv.index('H2O(L)')]  # [L/hr]*[g/L]/[g/mol] = [mol/hr]
            n_Hp            = aqvols[k] * (self.x[self.combined_var_space['H+ '+self.column[k]+'-'+num]])   # [L/hr]*[g/L]/[g/mol] = [mol/hr]
    
          
            if k==0: # extraction column: 
                #check if there's a parent column. if there isn't, use primary feed, otherwise, 
                # take the corresponding parent aqueous composition data

                
                if parent_col and pnum in self.modules.keys(): # current column has parent column, which actually exists in this configuration
                    myree =  np.array(self.y[parent_col][-self.nsy+self.canteravars.index('H+')+1:-self.norgy])  
                    n_H2O = self.nsp[parent_col][self.canteranames.index('H2O(L)')]   
                    n_Hp  = self.x[self.combined_var_space['H+ '+self.column[k]+'-'+num]] * n_H2O / (self.rhoslv[self.solv.index('H2O(L)')]/self.mwslv[self.solv.index('H2O(L)')])
                    for re in self.ree:                     
                        nre[self.ree.index(re)] = myree[self.ree.index(re)] #[mol/hr]

                else:
                    for re in self.ree:                     
                        nre[self.ree.index(re)] = self.ree_conc[self.ree.index(re)] * aqvols[k] / self.mwre[self.ree.index(re)] # [g/L]/[L/hr]/[g/mol] = [mol/hr]
    
            elif k==1:
                for re in self.ree:                     
                    nre[self.ree.index(re)] =  is_scrub[self.ree.index(re)] * self.x[self.combined_var_space['Recycle-'+num]] * self.ree_conc[self.ree.index(re)] * aqvols[k] / self.mwre[self.ree.index(re)] # [1]*[g/L]/[L/hr]/[g/mol] = [mol/hr]
                
                                                                                                                                       # 0.05 makes it a small value
            else:
                for re in self.ree:
                    nre[self.ree.index(re)] = 0.0              

            n_Cl            = 3*(sum(nre)) + n_Hp # Cl- mole balance, all REEs come in as chlorides from leaching
    
            n_specs         = [n_H2O,n_Hp,0,n_Cl]+[ii for ii in nre]+[n_HA,n_dodec] +[ij for ij in salts]
    
     
            # store in pandas dataframe
            self.nsp[self.column[k]+'-'+num]       = n_specs 
            self.nsp0[self.column[k]+'-'+num]      = n_specs

                
#%%
    def eval_column(self, num, col):
        
        """ This function evaluates the column to compute stream compositions
            for all stages """
       
        Ns                  = int(self.x[self.combined_var_space[col+'-'+num]]) # Ns (number of stages per column)

        # if Number of stages is zero, populate with default values        
        if Ns == 0:
            resy = rs.result_struct([0]*len(self.canteravars), None,'No stages', 10000)
            self.nsp[col+'-'+num]       = [0 for ii in range(len(self.nsp0[col+'-'+num]))]                     
                    
        else:
            ycol                = self.inity(col,num, Ns) # initialize y (stream vector)      
            resy               = root(eColOne, ycol, args=(self, num, col, Ns), method='df-sane', options=None) # method='hybr', options=None) #options={'disp':True, 'maxfev':15}    

        return resy


#%%
    def update_nsp(self, resy, prev_col, num):
        
        col_index = self.column.index(prev_col)+1
        
        if col_index <= 2:            
            col = self.column[col_index] #self.column.index(col_index)]
                        
            self.nsp[col+'-'+num][self.naq:]  =   [ resy.x[self.canteravars.index('(HA)2(org)')] ] + [self.nsp[col+'-'+num][self.canteranames.index('dodecane')] ]+\
                    [jk for jk in resy.x[self.canteravars.index('(HA)2(org)')+1:self.nsy]]  # org exit from previous column

#%%
    def objective(self, size=1):
        
        """ 
        Compute objective function
        in this case, make it recority
        sum(recovery[i] * purity[i])
        """
        recority = {}
        for key in self.recovery.keys():
            recority[key] = np.array(self.purity[key])*np.array(self.recovery[key])
        
        self.recority = recority
        fun = sum(sum(recority.values()))
        return fun/size
            
#%% 
    def constraintfun(self):
        
        """ 
        Compute constraint function
        use weighted penalty based on column type
        c_i = max(0, spec-value_i)
        """
        c_purity = {}
        c_recovery = {}

        for key in self.purity.keys():
            name, no = key.split('-')
            c_purity[key] = max(0, ((self.purity_spec -self.purity[key])*self.penalty[name]))
            c_recovery[key] = max(0, sum([self.recov_spec-ktem for ktem in self.recovery[key]])*self.penalty[name])

        self.constraint = {"purity":c_purity, "recovery":c_recovery}
        
        return sum(c_purity.values()),sum(c_recovery.values())


#%%      
    def evaluate(self, x,): #
        
        """ This is the simpler implementation of the process column design
            it avoids the need to converge recycle streams. For now, it is not
            adequate for multi-module processes (I don't yet have a reliable way of defining
            the equivalent recycle amount)"""

        self.x                  = [i for i in x] 
        
        
        self.status         = {} # all status
        self.msg            = {} # all msgs
        self.fun            = {} # all fun vals
        self.recycle        = {}

       # Assuming org volumetric flow is the same at every extraction stage - This is for initialization
       
        self.orgvol              = [self.vol * (self.x[self.combined_var_space['OA Extraction-'+self.mod_space[next(iter(self.mod_space))]]]) #---------------> should this go into loop;
                                        for ij in range(len(self.mod_space))]
        
        
        # Store all numbers of stages
        for name, num in self.mod_space.items():

#            self.orgvol              = [self.vol * (self.x[self.combined_var_space['OA Extraction-'+num]]) 
#                                for ij in range(len(self.mod_space))]

            ###########construct nsp - dataframe of feed species molar amounts ######################

            self.create_nsp_open(name, num) 


            for col in self.column:
                
                self.Ns.update( {col+'-'+num: int(self.x[self.combined_var_space[col+'-'+num]]) } )
                
#            ########################## Evaluate extraction ########################
            
                resy               = self.eval_column(num,col)        
        
                self.y.update({col+'-'+num:[ij for ij in resy.x]})
                self.status.update({col+'-'+num:resy.success})
                self.msg.update({col+'-'+num:resy.message})
                self.fun.update({col+'-'+num:resy.fun})

                self.update_nsp(resy, col, num)

#       # estimate reward         
        self.recovery_open()#(resy.x)     
        
        # estimate objective
        fun = self.objective()
        
        # estimate constraints
        c_pur, c_rec = self.constraintfun()
        
        return fun + c_pur + c_rec
      

############################################################################### 
            
 
#%%
    def inity(self, col, num, Ns): # Initialize y

        y_i = []

        for m in range(Ns): 

            y_i.append([self.nsp[col+'-'+num][self.canteranames.index(jk)] for jk in self.canteravars]) # pick values from nsp that corespond to the variable entries
        
        y = np.array([ii for ij in y_i for ii in ij]) 
                      
        return y
            


#%%
    def recovery_open(self):
        
        target = {}
        recovery = {}
        streams = {}
        feeds = {}
        purity = {}
    
        mod_nums = [key for key in config.valid_processes[self.config_key]['modules']]
    
        for key in config.valid_processes[self.config_key]['modules']:
        
            for col in self.column:
                
                if get_child(col+'-'+key, mod_nums) == None:
                
                    if col == 'Strip':
                        target[col+'-'+key] = config.valid_processes[self.config_key]['modules'][key]['strip_group'] # strip target
                    else:
                        target[col+'-'+key] = config.valid_processes[self.config_key]['modules'][key]['ext_group'] # non-strip target                                     
                    
                    streams[col+'-'+key] = self.y[col+'-'+key][-self.nsy+self.canteravars.index('H+')+1:-self.norgy]
                    feeds[col+'-'+key] = self.nsp[col+'-'+key][self.canteranames.index('Cl-')+1:self.naq].values
                    
                    pure                 = [ij/sum(streams[col+'-'+key]) for ij in streams[col+'-'+key]]
                    purity[col+'-'+key]  = sum([pure[self.ree.index(ik)] for ik in target[col+'-'+key] ] )
                    
    
        feed_input = [0 for item in feeds[col+'-'+key]]            
        for key, value in feeds.items():
            feed_input = [ij + jk for ij, jk in zip(feed_input, value)] 
    
        for key, value in streams.items():
            recover =  [ij/jk for ij, jk in zip(value, feed_input)]       
            recovery[key] = [recover[self.ree.index(ik)] for ik in target[key]]
                    
                    
        self.streams = streams
        self.feeds = feeds
        self.total_feeds = {ij:jk for ij, jk in zip(self.ree, feed_input)}
        self.purity = purity
        self.recovery = recovery
        self.target_rees = target
    

#%%
    def get_dimension(self):
        """Return the dimension of the problem.

        Returns
        -------
        int
            The dimension of the problem.
        """
        return self.dimension
    # -- end function
    
    def get_var_lower(self):        
        """Return the array of lower bounds on the variables.

        Returns
        -------
        List[float]
            Lower bounds of the decision variables.
        """
        return self.var_lower
    # -- end function
        
    def get_var_upper(self):
        """Return the array of upper bounds on the variables.

        Returns
        -------
        List[float]
            Upper bounds of the decision variables.
        """
        return self.var_upper
    # -- end function

    def get_var_type(self):
        """Return the type of each variable.
        
        Returns
        -------
        1D numpy.ndarray[char]
            An array of length equal to dimension, specifying the type
            of each variable. Possible types are 'R' for real
            (continuous) variables, and 'I' for integer (discrete)
            variables.
        """
        return self.var_type        
# -- end class

#####################################################################################################       
"""------------------------------HELPER FUNCTIONS------------------------------"""


#reverse key value pairs 
#%%
def reverse_dict(D):
	return {v: k for k, v in D.items()}



#combine a and b and return the new dictionary
def combine_dict(a, b):
	c = {}
	c.update(a)
	c.update(b)
	return c


def get_parent(var):
  '''
    takes variable name, determines module number, then returns the column_names
    and module number that feeds into that module_num

    Extraction-1 -> Extraction-0
	Extraction-4 -> Strip-1
    Strip-4 -> None
    Extraction-0 -> None
  '''

  name, num = var.split('-')
  num = int(num)
  if  name != 'Extraction' or num < 1:
    return None
  parentNum = (num-1)//2
#  print('child', var, '   parent', (f'Extraction-{parentNum}',f'Strip-{parentNum}')[num % 2 == 0])
  return (f'Strip-{parentNum}',f'Extraction-{parentNum}')[num % 2 == 0]



def get_child(var, mod_nums):
  '''
    takes current module number (var), get child number from 
    column (coltype), check if child number in mod_nums,
    return child otherwise return None

    Extraction-1 -> Extraction-4 if it exists
	Extraction-0 -> Extraction-2 if it exists
    Strip-0 -> Extraction-1 if it exists

  '''
 
  name, num = var.split('-')
  if name == 'Scrub':
      return None
  else:
      child = get_next(num,name)
      if str(child) in mod_nums:
          return 'Extraction-'+str(child)
      else:
          return None
      
        
        
def get_next(item_num, dir):
    '''
      returns next node number in tree in dir specified
          0      level: 0
        / \    
        1   2           1
      / \ / \
      3  4 5  6         2
  
      get_next(2, 'left') => 5
      get_next(1, 'right') => 4 
    '''
    item_num   = int(item_num)
    level, pos = get_level(item_num)
    next_level = pow(2, level+1) #items in next level
    add        = next_level-(pow(2, level)-pos)
    left       = item_num + add
    
    return (left, left+1)[dir=='right']


        
def get_level(curNum):
    '''
      helper function to get level and pos of a number in a binary tree
          0      level: 0
        / \    
        1   2           1
      / \ / \
      3  4 5  6         2
  
    get_level(2) => 1, 1
    get_level(4) => 2, 1 
    '''
    curNum = int(curNum)
    level = 0
    while pow(2, level+1) <= curNum+1:
      level+=1
    pos = curNum - pow(2, level) + 1 #0 indexed
    return level, pos


   

#%%
def eColOne (yI, obj, num, column, Ns, ree=[]) : # This should be equivalent to actual Lyon et al 4 recycling

    y                       = np.array([i for i in yI])
    naq                     = obj.naq
    nsy                     = obj.nsy
    norgy                   = obj.norgy
    stream                  = np.zeros(obj.ns)
    
###############################################################################
    
    count = 0 # adjust reference to column inlet stream  

    if Ns == 1: # single stage column
        stream[:naq]        = obj.nsp[column+'-'+num][:naq] # aqueous feed 
        stream[naq:]        = obj.nsp[column+'-'+num][naq:]
        obj.mix.species_moles    = stream.copy() 
        obj.mix.equilibrate('TP',log_level=0) # default  # maxsteps=100, maxiter=20
        y[count:count+nsy] = [obj.mix.species_moles[obj.canteranames.index(ij)] for ij in obj.canteravars]    #  # update aq and org variable vector 

        
    else:
    
        for i in range(1,Ns+1): # for each stage in column extraction   
            
            if i == 1: # feed aqueous stage 1  
                stream[:naq]        = obj.nsp[column+'-'+num][:naq] # aqueous feed 
                stream[naq:]        = [y[count+nsy+obj.canteravars.index('(HA)2(org)')]] +\
                                        [obj.nsp[column+'-'+num][obj.canteranames.index('dodecane')] ]+\
                                        [jk for jk in y[count+(nsy)+obj.canteravars.index('(HA)2(org)')+1:count+(2*nsy)] ]
                                        
            elif i == Ns : #feed organic stage N  
                stream[naq:]        = obj.nsp[column+'-'+num][naq:]
                
                cl_minus            = y[count-nsy+obj.canteravars.index('H+')] + 3*sum([jk for jk in y[count-nsy+obj.canteravars.index('H+')+1:count-norgy]])            
                stream[:naq]        = [obj.nsp[column+'-'+num][obj.canteranames.index('H2O(L)')]] +\
                                        [y[count-nsy+obj.canteravars.index('H+')]]+ [obj.nsp[column+'-'+num][obj.canteranames.index('OH-')] ]+\
                                        [cl_minus]+\
                                        [jk for jk in y[count-nsy+obj.canteravars.index('H+')+1:count-norgy]] # counts from after Cl up until the last ree
                                    # The use of count-nsy follows the logic that at the beginning of the second stage (N=2), count has been updated by nsy, but
                                    # this indexing for y is just starting (at zero); so to correct for this lag, we have to subtract nsy from count each time                                
                           
            else:
                cl_minus            = y[count-nsy+obj.canteravars.index('H+')] + 3*sum([jk for jk in y[count-nsy+obj.canteravars.index('H+')+1:count-norgy]])            
                stream[:naq]        = [obj.nsp[column+'-'+num][obj.canteranames.index('H2O(L)')]] +\
                                        [y[count-nsy+obj.canteravars.index('H+')]]+ [obj.nsp[column+'-'+num][obj.canteranames.index('OH-')] ]+\
                                        [cl_minus]+\
                                        [jk for jk in y[count-nsy+obj.canteravars.index('H+')+1:count-norgy]] # equivalent to count-nsy+naqy
    
                
                stream[naq:]        = [y[count+nsy+obj.canteravars.index('(HA)2(org)')]]+\
                                        [obj.nsp[column+'-'+num][obj.canteranames.index('dodecane')] ]+\
                                        [jk for jk in y[count+(nsy)+obj.canteravars.index('(HA)2(org)')+1:count+(2*nsy)] ]


            obj.mix.species_moles    = stream.copy() 
    
            obj.mix.equilibrate('TP',log_level=0) # default  
    
                        
            y[count:count+nsy] = [obj.mix.species_moles[obj.canteranames.index(ij)] for ij in obj.canteravars]    #  # update aq and org variable vector
            
            count += nsy

    return yI-y  


#%%
def flatten(iter):
  if type(iter) == dict:
    return np.array(list(iter.values())).reshape(1, -1)[0]
  elif type(iter) == list:
    return np.array(iter).reshape(1, -1)[0]


def get_sums(d,key):
  i = len(d) 	 #matrix rows
  j = len(d[key])  #matrix cols
  s = flatten(d) #1d list of dictionary `d`
  return list(np.sum(s.reshape(i, j), axis=0)) #sum cols of s





def highest_remaining(values, remove_n_lowest):
  '''
    values: numpy array containing values of components
  '''
  if type(values).__module__ == np.__name__: #values = numpy array
    vals = values.tolist()
  else:
    vals = values
  inds = [i for i in range(len(vals))]
  if remove_n_lowest > len(vals):
    indices = []
  elif remove_n_lowest >= 0:
    pairs = [(i, val) for i, val in zip(inds, vals)]

    remaining = sorted(pairs, key=operator.itemgetter(1), reverse=True)[:-remove_n_lowest]
    indices = [i for i, val in sorted(remaining, key=operator.itemgetter(0))]
  else:
    indices = inds
  return indices

