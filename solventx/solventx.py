

from scipy.optimize import root
import numpy as np
import cantera as ct
import pandas as pd
from solventx import result_struct
from solventx import utilities
import operator
import os


class solventx:
    
    immutable_var_names     = ['mol_frac']
    mutable_var_names       = ['H+ Extraction', 'H+ Scrub', 'H+ Strip', 'OA Extraction', 'OA Scrub', 'OA Strip', 'Recycle', 'Extraction', 'Scrub', 'Strip']
    feed_var_names          = ['(HA)2(org)']

#    var_names               = ['(HA)2(org)',	'H+ Extraction', 'H+ Scrub',	'H+ Strip',	'OA Extraction',	'OA Scrub',	'OA Strip', 'Recycle','Extraction','Scrub', 'Strip']#,	'Nd Scrub',	'Pr Scrub',	'Ce Scrub',	'La Scrub']#,	'Nd',	'Pr',	'Ce',	'La','Factor']
   
    coltypes                = ['Extraction','Scrub','Strip']
    #REEs                    = {1:['Nd','Pr','Ce','La'],2:['Nd','Pr'],3:['Ce','La'], 4:['Nd'],5:['Ce']}
    ml2l                    = 0.001 # mililiter to liter
    scale                   = 1 # volumetric scaling from unit ml/sec

    g_2_kg         = 1 #0.001 # grams to kg
    s_2_h          = 1.0/3600 # seconds to hours
    moduleID = '' #dummy variable
    target_conc =  1.0  #dummy variable
    
    def __init__(self,config_file,prep_capex=0, prep_opex=0, prep_revenue=0, prep_npv=0):
                 
        """Constructor.
        """
        
        self.confDict = utilities.read_config(config_file)
        solventxHome = self.confDict["solventxHome"]
        self.csvData = self.confDict['csvData']
        self.xmlData = self.confDict['xmlData']
        self.modulesData = self.confDict['modules']      
        
        self.dfree = pd.read_csv(os.path.join(solventxHome, self.csvData['dfree'])) # ree feed compositions (g/L)
        self.main_sx_feed = pd.read_csv(os.path.join(solventxHome, self.csvData['main_sx_feed'])) 
        self.mw = pd.read_csv(os.path.join(solventxHome, self.csvData['mw'])) 
        self.reecols = self.dfree.columns        
       
        
        for ktem in self.reecols:  #### Create custom data columns   #########################################
            self.dfree[ktem+' '+str(self.moduleID)] = self.dfree[ktem] * self.dfree['H2O Volume[L/hr]'] # g/L * L/hr = [g/hr]
        
        """ Load and join data from file """ #not be relevant to current problem
        self.inputs = pd.read_csv(os.path.join(solventxHome, self.csvData['ndprcela']))
        for item in self.inputs.columns:
            self.dfree[item] = self.inputs[item][0]
        
       
        self.xml = os.path.join(solventxHome,self.xmlData['xml'],self.xmlData['phase']+'_'+''.join(self.modulesData["input"])+'.xml')
        
        # Set required data
        self.df             = self.dfree              
        self.mainsxflow     = self.main_sx_feed   # kg/hr
                 
        self.phase_names    = self.confDict["phasenames"] # from xml input file
        self.phase          = ct.import_phases(self.xml,self.phase_names) 

        # Derived and/or reusable system parameters        
        self.column         = self.coltypes  #list(self.mod_input['Section']) # Column name
        self.solv           = self.confDict["solvents"] # solvent list -.i.e., electrolyte, extractant and organic diluent
        
        # ree by modules
        self.ree            = self.modulesData["input"] #self.REEs[self.moduleID] # (rare earth) metal list
        self.ree_strip      = self.modulesData["1"]["strip_group"] #self.REEs[self.stripGroup] # Strip target
        self.is_ree         = [1 if re in self.ree_strip else 0 for re in self.ree ]
        self.MID            = self.moduleID
                   
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

        # Feed volume
        self.ree_mass       = [self.df[ij+'Cl3 '+str(self.MID)][0]*(self.mw[ij][0]/self.mw[ij][0]) for ij in self.ree] # kg/hr of ree chloride
        mree                = sum(self.ree_mass) # kg/hr
        self.vol        = mree / self.g_2_kg / self.target_conc   # [kg/hr]*[g/kg] /[g/L] = [L/hr] 
        self.df['H2O Volume[L/hr] '+str(self.MID)] = self.vol
        
        self.purity_spec    = .99 # not needed?
        self.recov_spec     = .99 # not needed?
        
        self.revenue        = [0,0,0]
        self.Ns             = [0,0,0]

        self.nsp            = pd.DataFrame() # feed streams (aq and org) for each column       
        self.nsp0           = pd.DataFrame() # feed streams (aq and org) for each column       

        self.y              = {} # all compositions
        self.Ns             = {}

    # -- end function
  

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



    def create_var_space(self, x, components=4, input_feeds=1,): # Creates containers for process variables

        self.x = x # do you prefer list or nparray?

        var_space = {
            'immutable': {},
            'mutable':   {}, #var, index in obj.variables
        }
        
        mod_space = {}
        x_space = {}
        
        immutable_var_names = ['mol_frac']
        mutable_var_names   = ['H+ Extraction', 'H+ Scrub', 'H+ Strip', 'OA Extraction', 'OA Scrub', 'OA Strip', 'Recycle', 'Extraction', 'Scrub', 'Strip']
        feed_var_names = ['(HA)2(org)']
        lenx = len(mutable_var_names + feed_var_names)

        index = 0
        for feed_num in range(input_feeds):
            var_space['mutable'][f'{feed_var_names[0]}-{feed_num}'] = index
            index += 1 
            
        count = 0
        for module_num in range(components-1):
            
            mod_space[f'module-{module_num}'] = module_num
            x_space[f'module-{module_num}'] = self.x[count:count+lenx]
            count += lenx
            
            for i, var in enumerate(mutable_var_names):
                var_space['mutable'][f'{var}-{module_num}'] = index
                index += 1

        for comp_num in range(components):
            var_space['immutable'][f'{immutable_var_names[0]}-{comp_num}'] = index
            index += 1
        self.var_space = var_space
        mutable = self.var_space['mutable']
        immutable = self.var_space['immutable']
        self.combined_var_space = combine_dict(mutable, immutable)
        
        self.mod_space = mod_space
        self.x_space = x_space
        self.num_feeds = len(mod_space) * len(self.column)



		
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


    

    def create_nsp_open(self, module, lim, init_recyc=0.0, g_2_kg=g_2_kg ): #, h0, target, xml, cantera_data, experiments):
        
        """ Create mole flows for column feed streams in a given module. Arranges
            them in the dataframe, nsp, in the form that Cantera expects """
               
        nre                 = np.zeros(len(self.ree)) # REE specie moles
        salts               = np.zeros(len(self.ree)) # neutral complexe moles
        name, num           = module.split('-')        

                       
#        Flows
        orgvol              = self.orgvol[self.mod_space[module]]
        aqvols              = [orgvol/(self.x[self.combined_var_space['OA Extraction-'+num]]),orgvol/(self.x[self.combined_var_space['OA Scrub-'+num]]), orgvol/(self.x[self.combined_var_space['OA Strip-'+num]]) ]
        

#        Compositions   
        vol_HA              = orgvol *  (self.x[self.combined_var_space['(HA)2(org)-0']])
        vol_dodec           = orgvol - vol_HA
        n_HA                = vol_HA * self.rhoslv[self.solv.index('(HA)2(org)')]/self.mwslv[self.solv.index('(HA)2(org)')] # [L/hr]*[g/L]/[g/mol] = [mol/hr]
        n_dodec             = vol_dodec * self.rhoslv[self.solv.index('dodecane')]/self.mwslv[self.solv.index('dodecane')] # [L/hr]*[g/L]/[g/mol] = [mol/hr]            

        for k in range(len(self.column)): # k represents different columns - extraction, scrub or strip

            # check for module ID, then determine upstream source
            parent_col = get_parent(self.column[k]+'-'+num) 
            
            n_H2O           = aqvols[k] * self.rhoslv[self.solv.index('H2O(L)')]/self.mwslv[self.solv.index('H2O(L)')]  # [L/hr]*[g/L]/[g/mol] = [mol/hr]
            n_Hp            = aqvols[k] * (self.x[self.combined_var_space['H+ '+self.column[k]+'-'+num]])   # [L/hr]*[g/L]/[g/mol] = [mol/hr]

          
            if k==0:
                for re in self.ree:                     
                    nre[self.ree.index(re)] = self.ree_mass[self.ree.index(re)]/g_2_kg / self.mwre[self.ree.index(re)] # [kg/hr]/[kg/g]/[g/mol] = [mol/hr]

            elif k==1:
                uncle = get_parent(self.column[k-1]+'-'+num) # if a parent module exists for ext column  
                if uncle:                       
                    myree =  np.array(self.y[uncle][-self.nsy+self.canteravars.index('H+')+1:-self.norgy])

                    target_indices = highest_remaining(myree, lim)
                    print('target_indices', target_indices)                     
                    target_ree = myree[target_indices[:round( (len(myree)-lim)/2 )]] # replace 2 with something automatic later
                    print('target_ree', target_ree)
                    my_is_ree = [1 if ij in target_ree else 0 for ij in myree]
                    
                    for re in self.ree: 
                        nre[self.ree.index(re)] = my_is_ree[self.ree.index(re)] *self.x[self.combined_var_space['Recycle-'+num]] *init_recyc* myree[self.ree.index(re)] # = [mol/hr]                   
                else:
                    
                    for re in self.ree:                     
                        nre[self.ree.index(re)] = self.is_ree[self.ree.index(re)] *self.x[self.combined_var_space['Recycle-'+num]] *init_recyc* self.ree_mass[self.ree.index(re)]/g_2_kg / self.mwre[self.ree.index(re)] # [kg/hr]/[kg/g]/[g/mol] = [mol/hr]
                                                                                                                        # 0.05 makes it a small value
            else:
                for re in self.ree:
                    nre[self.ree.index(re)] = 0.0                
                            
            n_Cl            = 3*(sum(nre)) + n_Hp # Cl- mole balance, all REEs come in as chlorides from leaching


            # check for module ID, then determine upstream source
            if parent_col: 
                afa, okwa           = parent_col.split('-')        
                
                # mass balance for chloride ion
                if afa == 'Strip': # if strip column, subtract the recycled portion to estimate equivalent recycle ratio 

                    strips = np.array(self.y[parent_col][-self.nsy+self.canteravars.index('H+')+1:-self.norgy])
                    targ_ind = highest_remaining(strips, lim)
                    strip_ree = strips[targ_ind[0]]
                    scrub_ree = self.nsp['Scrub-'+okwa][self.canteranames.index('Cl-')+1+targ_ind[0]] # corresponding index in scrub-in

                    recycle = scrub_ree/strip_ree # equivalent recycle ratio
                    self.recycle.update({module:recycle})
                    
                    nre = [ij*(1-recycle) for ij in self.y[parent_col][-self.nsy+self.canteravars.index('H+')+1:-self.norgy]] # no of moles of rare earth
                    n_H2O = self.nsp['Strip-'+okwa][self.canteranames.index('H2O(L)')]*(1-recycle)
                else:
                    nre = self.y[parent_col][-self.nsy+self.canteravars.index('H+')+1:-self.norgy]
                    n_H2O =self.nsp[parent_col][self.canteranames.index('H2O(L)')]
                  
                n_Cl = n_Hp + 3*sum(nre)

            n_specs         = [n_H2O,n_Hp,0,n_Cl]+[ii for ii in nre]+[n_HA,n_dodec] +[ij for ij in salts]

 
            # store in pandas dataframe
            self.nsp[self.column[k]+'-'+num]       = n_specs 
            self.nsp0[self.column[k]+'-'+num]      = n_specs 





    def create_nsp_loop(self, module, lim, init_recyc=0.0, g_2_kg=g_2_kg ): #, h0, target, xml, cantera_data, experiments):
        
        """ Create mole flows for column feed streams in a given module. Arranges
            them in the dataframe, nsp, in the form that Cantera expects """
               
        nre                 = np.zeros(len(self.ree)) # REE specie moles
        salts               = np.zeros(len(self.ree)) # neutral complexe moles
        name, num           = module.split('-')        

                       
#        Flows
        orgvol              = self.orgvol[self.mod_space[module]]
        aqvols              = [orgvol/(self.x[self.combined_var_space['OA Extraction-'+num]]),orgvol/(self.x[self.combined_var_space['OA Scrub-'+num]]), orgvol/(self.x[self.combined_var_space['OA Strip-'+num]]) ]
        

#        Compositions   
        vol_HA              = orgvol *  (self.x[self.combined_var_space['(HA)2(org)-0']])
        vol_dodec           = orgvol - vol_HA
        n_HA                = vol_HA * self.rhoslv[self.solv.index('(HA)2(org)')]/self.mwslv[self.solv.index('(HA)2(org)')] # [L/hr]*[g/L]/[g/mol] = [mol/hr]
        n_dodec             = vol_dodec * self.rhoslv[self.solv.index('dodecane')]/self.mwslv[self.solv.index('dodecane')] # [L/hr]*[g/L]/[g/mol] = [mol/hr]            

        for k in range(len(self.column)):

            # check for module ID, then determine upstream source
            parent_col = get_parent(self.column[k]+'-'+num) 
            
            n_H2O           = aqvols[k] * self.rhoslv[self.solv.index('H2O(L)')]/self.mwslv[self.solv.index('H2O(L)')]  # [L/hr]*[g/L]/[g/mol] = [mol/hr]
            n_Hp            = aqvols[k] * (self.x[self.combined_var_space['H+ '+self.column[k]+'-'+num]])   # [L/hr]*[g/L]/[g/mol] = [mol/hr]

          
            if k==0:
                for re in self.ree:                     
                    nre[self.ree.index(re)] = self.ree_mass[self.ree.index(re)]/g_2_kg / self.mwre[self.ree.index(re)] # [kg/hr]/[kg/g]/[g/mol] = [mol/hr]


            elif k==1:
                uncle = get_parent(self.column[k-1]+'-'+num) # if a parent module exists for ext column  
                if uncle:                       
                    myree =  np.array(self.y[uncle][-self.nsy+self.canteravars.index('H+')+1:-self.norgy]) #     maxids = np.array(myree).argsort()[-2:][::-1]

                    # Replace section with something more automatated
                    target_indices = highest_remaining(myree, lim) 
                    target_ree = myree[target_indices[:round( (len(myree)-lim)/2 )]] # replace 2 with something automatic later
                    my_is_ree = [1 if ij in target_ree else 0 for ij in myree]
                    self.tindex.update({module:target_indices})
                    for re in self.ree: 
                        nre[self.ree.index(re)] = my_is_ree[self.ree.index(re)] *self.x[self.combined_var_space['Recycle-'+num]] *init_recyc* myree[self.ree.index(re)] # = [mol/hr]                   
                else:
                    
                    for re in self.ree:                     
                        nre[self.ree.index(re)] = self.is_ree[self.ree.index(re)] *self.x[self.combined_var_space['Recycle-'+num]] *init_recyc* self.ree_mass[self.ree.index(re)]/g_2_kg / self.mwre[self.ree.index(re)] # [kg/hr]/[kg/g]/[g/mol] = [mol/hr]
                                                                                                                        # 0.05 makes it a small value
            else:
                for re in self.ree:
                    nre[self.ree.index(re)] = 0.0                
                            
            n_Cl            = 3*(sum(nre)) + n_Hp # Cl- mole balance, all REEs come in as chlorides from leaching


           
            if parent_col: 
                afa, okwa           = parent_col.split('-')        
                
                if afa == 'Strip':
                    nre = [ij*(1-self.x[self.combined_var_space['Recycle-'+num]]) for ij in self.y[parent_col][-self.nsy+self.canteravars.index('H+')+1:-self.norgy]]
                    n_H2O = self.nsp['Strip-'+okwa][self.canteranames.index('H2O(L)')]*(1-self.x[self.combined_var_space['Recycle-'+num]])
                    self.recycle.update({module:self.x[self.combined_var_space['Recycle-'+num]]})

                else:
                    nre = self.y[parent_col][-self.nsy+self.canteravars.index('H+')+1:-self.norgy]
                    n_H2O = self.nsp[parent_col][self.canteranames.index('H2O(L)')]
                   
                n_Cl = n_Hp + 3*sum(nre)

            n_specs         = [n_H2O,n_Hp,0,n_Cl]+[ii for ii in nre]+[n_HA,n_dodec] +[ij for ij in salts]

            # store in pandas dataframe
            self.nsp[self.column[k]+'-'+num]       = n_specs 
            self.nsp0[self.column[k]+'-'+num]      = n_specs 


                
                

    def eval_column(self, module, col):
        
        """ This function evaluates the column to compute stream compositions
            for all stages """

        name, num           = module.split('-')        
        Ns                  = int(self.x[self.combined_var_space[col+'-'+num]]) # Ns (number of stages per column)

        # if Number of stages is zero, populate with default values        
        if Ns == 0:
            resy = rs.result_struct([0]*len(self.canteravars), None,'No stages', 0)
            
        else:
            ycol                = self.inity(col,num, Ns) # initialize y (stream vector)      
            resy               = root(eColOne, ycol, args=(self, num, col, Ns), method='df-sane', options=None) # method='hybr', options=None) #options={'disp':True, 'maxfev':15}    

        return resy




    def evaluate_open(self, x, lim=2, g_2_kg=g_2_kg ): #
        
        """ This is the simpler implementation of the process column design
            it avoids the need to converge recycle streams. For now, it is not
            adequate for multi-module processes (I don't yet have a reliable way of defining
            the equivalent recycle amount)"""

        self.x                  = [i for i in x] 
        #        MID = ' '+str(self.MID)
        
        self.status         = {} # all status
        self.msg            = {} # all msgs
        self.fun            = {} # all fun vals
        self.recycle        = {}

        # Store all numbers of stages
        for key in self.mod_space:
            name, num           = key.split('-')
            for item in self.column:
                self.Ns.update( {item+'-'+num: int(self.x[self.combined_var_space[item+'-'+num]]) } )
        
        # Assuming org volumetric flow is the same at every extraction stage - This is for initialization
        self.orgvol              = [self.vol * (self.x[self.combined_var_space['OA Extraction-0']]) 
                                        for ij in range(len(self.mod_space))]
        
########################construct nsp - dataframe of feed species molar amounts ######################

        
        for module in self.mod_space:

            name, num           = module.split('-')

            self.create_nsp_open(module, lim, init_recyc=1)

            # if No extraction column #########################################
            if self.x[self.combined_var_space['Extraction-'+num]] == 0:
                
                rez = rs.result_struct([0]*len(self.canteravars), True,'No stages', 0)
            
                for item in self.column:
                    
                    self.y.update({item+'-'+num:rez.x})
                    self.status.update({item+'-'+num:rez.status})
                    self.msg.update({item+'-'+num:rez.message})
                    self.fun.update({item+'-'+num:rez.fun})

            else:    
        ########################## Evaluate extraction ########################
        
                resye               = self.eval_column(module,'Extraction')        
        
                self.y.update({'Extraction-'+num:[ij for ij in resye.x]})
                self.status.update({'Extraction-'+num:resye.success})
                self.msg.update({'Extraction-'+num:resye.message})
                self.fun.update({'Extraction-'+num:resye.fun})
    
        
        ########################## Evaluate scrub  ##############################
                
                self.nsp['Scrub-'+num][self.naq:]  =   [ resye.x[self.canteravars.index('(HA)2(org)')] ] + [self.nsp['Scrub-'+num][self.canteranames.index('dodecane')] ]+\
                                    [jk for jk in resye.x[self.canteravars.index('(HA)2(org)')+1:self.nsy]]  # org exit from feed stage 1
        
                    
                resysc               = self.eval_column(module, 'Scrub')        
    
                self.y.update({'Scrub-'+num:[ij for ij in resysc.x]})
                self.status.update({'Scrub-'+num:resysc.success})
                self.msg.update({'Scrub-'+num:resysc.message})
                self.fun.update({'Scrub-'+num:resysc.fun})
        
        ########################## Evaluate strip ##############################
        
                self.nsp['Strip-'+num][self.naq:]  =   [ resysc.x[self.canteravars.index('(HA)2(org)')] ] + [self.nsp['Strip-'+num][self.canteranames.index('dodecane')] ]+\
                                        [jk for jk in resysc.x[self.canteravars.index('(HA)2(org)')+1:self.nsy]]  # org exit from feed 'Strip'age 1
        
        
                resyst               = self.eval_column(module,'Strip')        
    
                self.y.update({'Strip-'+num:[ij for ij in resyst.x]})
                self.status.update({'Strip-'+num:resyst.success})
                self.msg.update({'Strip-'+num:resyst.message})
                self.fun.update({'Strip-'+num:resyst.fun})


                self.recovery_open()#(resy.x)        
                self.recovery_loop()#(resy.x)        

############################################################################### 
            
        


    def evaluate_loop(self, x, lim=0, g_2_kg=g_2_kg ): #, h0, target, xml, cantera_data, experiments):
        """ This is a more representative implementation of the process column design
            it explicitly converges recycle streams - but this makes it very expensive."""

        self.x                  = [i for i in x] 

        self.status         = {} # all status
        self.msg            = {} # all msgs
        self.fun            = {} # all fun vals
        self.tindex         = {} # target indices for nre recycle convergence for downstream modules
        self.recycle        = {}


        # Store all stage numbers

        for key in self.mod_space:
            name, num           = key.split('-')
            for item in self.column:
                self.Ns.update( {item+'-'+num: int(self.x[self.combined_var_space[item+'-'+num]]) } )

            
        # Assuming org is the same at every extraction stage - This is for initialization
        self.orgvol              = [self.vol * (self.x[self.combined_var_space['OA Extraction-0']]) 
                                        for ij in range(len(self.mod_space))]
        
########################construct nsp#  dataframe of feed species molar amounts ######################

        for module in self.mod_space:

            name, num           = module.split('-')

            self.create_nsp_loop(module, lim, init_recyc=0.05)
            

            # if No extraction column #########################################
            if self.x[self.combined_var_space['Extraction-'+num]] == 0:
                
                rez = rs.result_struct([0]*len(self.canteravars), True,'No stages', 0)
            
                for item in self.column:
                    
                    self.y.update({item+'-'+num:rez.x})
                    self.status.update({item+'-'+num:rez.status})
                    self.msg.update({item+'-'+num:rez.message})
                    self.fun.update({item+'-'+num:rez.fun})                   

            else:    
                    
    ########################## Evaluate extraction ############################
    
                resye               = self.eval_column(module, 'Extraction')        
        
                self.y.update({'Extraction-'+num:[ij for ij in resye.x]})
                self.status.update({'Extraction-'+num:resye.success})
                self.msg.update({'Extraction-'+num:resye.message})
                self.fun.update({'Extraction-'+num:resye.fun})
        
        ########################## Evaluate scrub & strip #####################
                if num == '0' or lim==0:
                    nre                 = np.array([ij for ij in self.nsp['Scrub-'+num][self.canteranames.index('Cl-')+1:self.naq-lim]]) # only rare earths                    
                    converge            = root(tearcolumns, nre, args=(self, resye.x, lim, num, module, 'Scrub','Strip'), method='excitingmixing', options=None)

                else:
                    # replace with something that scales with no. of components
                    the_rees            = np.array([ij for ij in self.nsp['Scrub-'+num][self.canteranames.index('Cl-')+1:self.naq]])#-lim]])
                    nre                 = the_rees[self.tindex[module][:round( (len(the_rees)-lim)/2 )]]
                    converge            = root(tearcol, nre, args=(self, resye.x, num, module, 'Scrub','Strip', lim), method='excitingmixing', options=None)

        #############################update y vector of stream compositions####      
                count = 1
                for item in self.resy:
    
                    self.y.update({self.column[count]+'-'+num:[ij for ij in item.x]})
                    self.status.update({self.column[count]+'-'+num:item.success})
                    self.msg.update({self.column[count]+'-'+num:item.message})
                    self.fun.update({self.column[count]+'-'+num:item.fun})
                    count += 1
                    
                self.converge = converge                
                
                self.recovery_open()#(resy.x)        
                self.recovery_loop()#(resy.x)        
        



    def inity(self, col, num, Ns): # Initialize y

        y_i = []

        for m in range(Ns): 

            y_i.append([self.nsp[col+'-'+num][self.canteranames.index(jk)] for jk in self.canteravars]) # pick values from nsp that corespond to the variable entries
        
        y = np.array([ii for ij in y_i for ii in ij]) 
                      
        return y
            


    def  recovery_open(self):
        
        col_out = {} 
        max_pur, purity  = {}, {}
        max_recov, recovery = {}, {}
        argmax = {}
        reemax = {}
        feed_in   = {} #self.nsp['Extraction-0'][self.canteranames.index('Cl-')+1:self.naq]        
        scrub_in  = {} # for closed, do not use scrub in

        parents = [get_parent(item) for item in self.y.keys() if get_parent(item) != None]            
        
        for key, value in self.y.items():
#            parent_col = get_parent(key)            
            if key not in parents:               
                feed_in.update({key: [ij for ij in self.nsp['Extraction-0'][self.canteranames.index('Cl-')+1:self.naq] ]})
                stream = value[-self.nsy+self.canteravars.index('H+')+1:-self.norgy]
                col_out.update({key:stream})

                purity.update({key: [i/sum(stream) for i in stream]})
                max_pur.update({key:max(purity[key])})
                argmax.update({key:np.argmax(purity[key])})
                reemax.update({key:self.ree[argmax[key]]})
                
                name, num           = key.split('-')        
                if name == 'Scrub':
                    scrub_in.update({key:[ij for ij in self.nsp['Scrub-'+num][self.canteranames.index('Cl-')+1:self.naq] ]})
                else:
                    scrub_in.update({key:[0 for i in self.nsp['Scrub-'+num][self.canteranames.index('Cl-')+1:self.naq]]})
                   
        scrubs = get_sums(scrub_in,'Scrub-0') 
        feeds  = get_sums(feed_in,'Scrub-0')
        feeds = [ij/len(feed_in) for ij in feeds]
        
        ins = [ij+jk for ij,jk in zip(feeds, scrubs)]
        
        for key, value in col_out.items():
            recovery.update({key:[ij/jk for ij,jk in zip(col_out[key],ins) ]})
            
            max_recov.update({key: recovery[key][argmax[key]]})

         
        self.recovery = recovery
        self.max_recov = max_recov
        self.ree_max   = reemax
        
        self.purity  = purity
        self.max_pur = max_pur
        self.raffinates = col_out
        
        self.feed_in    = feed_in
        self.parent_cols = parents
        
        self.total_feed = feeds
        self.total_scrubs = scrubs




    def  recovery_loop(self):
        
        col_out = {} 
        max_pur, purity  = {}, {}
        max_recov, recovery = {}, {}
        argmax = {}  
        reemax = {}
        feed_in   = {} #self.nsp['Extraction-0'][self.canteranames.index('Cl-')+1:self.naq]        
        scrub_in  = {} # for closed, do not use scrub in

        parents = [get_parent(item) for item in self.y.keys() if get_parent(item) != None]            
        
        for key, value in self.y.items():

            if key not in parents:                
                feed_in.update({key: [ij for ij in self.nsp['Extraction-0'][self.canteranames.index('Cl-')+1:self.naq]]})
                stream_0 = value[-self.nsy+self.canteravars.index('H+')+1:-self.norgy]

                name, num           = key.split('-')        
                if name == 'Strip':
                    stream = [ij-jk for ij,jk in zip(value[-self.nsy+self.canteravars.index('H+')+1:-self.norgy], self.nsp['Scrub-'+num][self.canteranames.index('Cl-')+1:self.naq])]
                else:
                    stream = stream_0.copy() #value[-self.nsy+self.canteravars.index('H+')+1:-self.norgy]
                    
                col_out.update({key:stream})
                purity.update({key: [i/sum(stream_0) for i in stream_0]})
                max_pur.update({key:max(purity[key])})
                argmax.update({key:np.argmax(purity[key])})
                reemax.update({key:self.ree[argmax[key]]})

                name, num           = key.split('-')        
                if name == 'Scrub':
                    scrub_in.update({key:[ij for ij in self.nsp['Scrub-'+num][self.canteranames.index('Cl-')+1:self.naq] ]})
                else:
                    scrub_in.update({key:[0 for ij in self.nsp['Scrub-'+num][self.canteranames.index('Cl-')+1:self.naq]]})


        scrubs = get_sums(scrub_in,'Scrub-0')             
        feeds  = get_sums(feed_in,'Scrub-0')
        feeds = [ij/len(feed_in) for ij in feeds]
        
                
        for key, value in col_out.items():
            recovery.update({key:[ij/jk for ij,jk in zip(col_out[key],feed_in[key])]})
            
            max_recov.update({key: recovery[key][argmax[key]]})
         
        self.recovery_ = recovery
        self.max_recov_ = max_recov
        self.ree_max_   = reemax
        
        self.purity_  = purity
        self.max_pur_ = max_pur
        self.raffinates_ = col_out
        
        self.feed_in_    = feed_in
        self.parent_cols_ = parents
        
        self.total_feed = feeds
        self.total_scrubs = scrubs

      
        
# -- end class

#####################################################################################################       
#------------------------------HELPER FUNCTIONS------------------------------


#reverse key value pairs 
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
  return (f'Extraction-{parentNum}',f'Strip-{parentNum}')[num % 2 == 0]




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




def tearcolumns(nre, obj, y, lim, num, module, sc, st): #ree
   
#    print('len y in tear', len(y))
    
    resy                = []
    recycle             = obj.x[obj.combined_var_space['Recycle-'+num]]
    
    # update scrub feed with recycled stream
    obj.nsp[sc+'-'+num] = [ij for ij in obj.nsp0[sc+'-'+num]] # update with primary scrub feed - no address sharing
    

    obj.nsp[sc+'-'+num][obj.canteranames.index('Cl-')+1:obj.naq-lim] = [max(0,ij) for ij in nre] # update ree !!!!!!!!!!REMOVED +=!!!!!!!!!!
    obj.nsp[sc+'-'+num][obj.canteranames.index('Cl-')] =  obj.nsp[sc+'-'+num][obj.canteranames.index('H+')] + 3*sum(obj.nsp[sc+'-'+num][obj.canteranames.index('Cl-')+1:obj.naq] )# update cl
    
    
    obj.nsp[sc+'-'+num][obj.naq:]  =   [ y[obj.canteravars.index('(HA)2(org)')] ] + [obj.nsp[sc+'-'+num][obj.canteranames.index('dodecane')] ]+\
                            [jk for jk in y[obj.canteravars.index('(HA)2(org)')+1:obj.nsy]]  # org exit from feed stage 1

        
    resy.append( obj.eval_column(module, sc) )
    y1                  = resy[0].x

    ####################################################################

    obj.nsp[st+'-'+num][obj.naq:]  =   [ y1[obj.canteravars.index('(HA)2(org)')] ] + [obj.nsp[st+'-'+num][obj.canteranames.index('dodecane')] ]+\
                            [jk for jk in y1[obj.canteravars.index('(HA)2(org)')+1:obj.nsy]]  # org exit from feed stage 1

    resy.append( obj.eval_column(module, st) )
    
    obj.resy = resy
    
    # compute strip exit ree molar amount
    rhs = resy[1].x[-obj.nsy+obj.canteravars.index('H+')+1:-obj.norgy-lim] # 
    
#    print (sum(recycle * rhs - nre))

    return (recycle * rhs) - nre



def tearcol(nre, obj, y, num, module, sc, st, lim=0): #ree
   
#    print('len y in tear', len(y))
#    index               = obj.tindex[module]
    resy                = []
    recycle             = obj.x[obj.combined_var_space['Recycle-'+num]]
    
    # update scrub feed with recycled stream
    obj.nsp[sc+'-'+num] = [ij for ij in obj.nsp0[sc+'-'+num]] # update with primary scrub feed - no address sharing
    
    nindices = obj.canteranames.index('Cl-')+1 + np.array(obj.tindex[module][:round( (len(obj.ree)-lim)/2 )])[0] # first ree is the target one from recycle options
    
    obj.nsp[sc+'-'+num][nindices] = np.array([max(0,ij) for ij in nre])[0]  # update ree 
    obj.nsp[sc+'-'+num][obj.canteranames.index('Cl-')] =  obj.nsp[sc+'-'+num][obj.canteranames.index('H+')] + 3*sum(obj.nsp[sc+'-'+num][obj.canteranames.index('Cl-')+1:obj.naq] )# update cl

    
    obj.nsp[sc+'-'+num][obj.naq:]  =   [ y[obj.canteravars.index('(HA)2(org)')] ] + [obj.nsp[sc+'-'+num][obj.canteranames.index('dodecane')] ]+\
                            [jk for jk in y[obj.canteravars.index('(HA)2(org)')+1:obj.nsy]]  # org exit from feed stage 1

        
    resy.append( obj.eval_column(module, sc) )
    y1                  = resy[0].x

    ####################################################################

    obj.nsp[st+'-'+num][obj.naq:]  =   [ y1[obj.canteravars.index('(HA)2(org)')] ] + [obj.nsp[st+'-'+num][obj.canteranames.index('dodecane')] ]+\
                            [jk for jk in y1[obj.canteravars.index('(HA)2(org)')+1:obj.nsy]]  # org exit from feed stage 1

    resy.append( obj.eval_column(module, st) )
    
    obj.resy = resy
    
    # compute strip exit ree molar amount
    yindices = -obj.nsy+obj.canteravars.index('H+')+1 + np.array(obj.tindex[module][:round( (len(obj.ree)-lim)/2)])[0]
    rhs = resy[1].x[yindices] # 
    
    return (recycle * rhs) - nre


    
    


def eColOne (yI, obj, num, column, Ns, ree=[]) : # This should be equivalent to actual Lyon et al 4 recycling

    y                       = np.array([i for i in yI])
    naq                     = obj.naq
    nsy                     = obj.nsy
    norgy                   = obj.norgy
    stream                  = np.zeros(obj.ns)
    
###############################################################################
#    print(Ns, len(y))
    
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
    
    #        try:
            obj.mix.equilibrate('TP',log_level=0) # default  
    
                        
    #        obj.mix.equilibrate('TP',log_level=0, maxsteps=100, maxiter=20) # default                      
            y[count:count+nsy] = [obj.mix.species_moles[obj.canteranames.index(ij)] for ij in obj.canteravars]    #  # update aq and org variable vector
    #        except:
            
            count += nsy

    return yI-y  



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

""" root(tear) timings (seconds):
    
lm: 39
broyden1: 11
broyden2: no convergence
anderson: 11
linearmixing: no convergence
diagbroyden:10
excitingmixing:4
krylov:8
hybr: 13
df-sane: 12

"""