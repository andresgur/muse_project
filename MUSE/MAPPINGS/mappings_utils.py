# @Author: Maxime Parra
# @Date:   07-07-2021
# @Email:  maxime.parra@irap.omp.eu
# @Last modified by:   mparra
# @Last modified time: 19-07-2021

import os
import sys
import numpy as np
sys.path.append('/home/mparra/PARRA/Scripts/Python/MUSE/BPT')
import pandas as pd
import glob
import matplotlib as mpl

"""
To do : implement a better plotting/single_grid for cut models (i.e. Alarie and Allen2008_cut) 
in order to see correctly the variations of time/distance instead of directly plotting everything 
(which might lead to wrong segmentation in the meshes)
"""

#change this to your own database path OR the github url for the tables
shock_database_path='/home/mparra/PARRA/prog/MAPPINGS/V/shock_database/tables/'


def plot_dgrid(dgrid,axis,BPT,color,label):
    
    '''
    function to plot the 2d grid of a shock_model in a specific BPT diagram
    
    alpha decreases (transparency increases) with increasing shock velocities
    
    we use two orders of the tables to draw separately (and with different alphas) the segments of the 2d surface 
    between each shck_vel and each mag_fld parameter
    
    A single label is displayed for an alpha of 1
    '''
    
    #sorting the grid in the two directions
    dgrid_v=dgrid.sort_values(['mag_fld','shck_vel'])
    dgrid_h=dgrid.sort_values(['shck_vel','mag_fld'])
    
    #getting the range of parameters
    vel_range=np.unique(dgrid.shck_vel)
    b_range=np.unique(dgrid.mag_fld)
    
    #and our alpha array
    alpha_arr=np.logspace(0,-1,len(vel_range))
    
    #first direction of the mesh, here we can do a single loop since our line are grouped by shock velocity and thus the 
    #grid_h segments have a single alpha
    for i in range(len(vel_range)):
        grid_h=dgrid_h[dgrid_h['shck_vel']==vel_range[i]]
        grid_h.plot(ax=axis,x=BPT.x_column,y=BPT.y_column,color=color,alpha=alpha_arr[i],legend=False,label=label if i==0 else '')
    
    #second direction of the mesh, here we need to do two loops since every segment in a line has a different alpha
    for i in range(len(b_range)):
        grid_v=dgrid_v[dgrid_v['mag_fld']==b_range[i]]
        for j in range(len(vel_range)):
            segment_v=grid_v.iloc[j:j+2]
            segment_v.plot(ax=axis,x=BPT.x_column,y=BPT.y_column,color=color,alpha=alpha_arr[j],legend=False,label='')
            
def plot_shock_mod(models,axis,BPT,labels):
    
    '''
    function to plot a list of models in a specific BPT diagram
    
    "models" should be either a dataframe or a list of dataframes
    In order to plot 2d surfaces correctly, keep one density/abundance/carbon to dust ratio for all the dataframes
    Also preferably have only one parameter varying in the array (since there is only one color scale)
    (e.g. densities or abundances)
    
    we use the plot_dgrid function to plot with variable transparency (here increasing with the shock speed values)
    
    viridis color map by default
    
    labels is assumed to be the same length as the models 
    '''
    viridis=mpl.cm.get_cmap('viridis',len(models))
    col_vars=viridis(range(len(models)))
    
    for i in range(len(models)):    
        curr_map_model=models[i]
        plot_dgrid(curr_map_model,axis,BPT,col_vars[i],labels[i])
            
def shock_model(modelname):
      
    """
    General function for the MAPPINGSV models, which uses the tables from the 3mdb database downloaded with gridmaker.
    Contains subclasses of different models depending on the abundance
    
    Parameter
    ----------
    model: str
        Gives which type of model is used here. Use :
            -'Alarie' for Alarie19s
            -'Allen' or '2008' for Allen2008
            -'cut' for Allen2008-cut
            -'Gutkin' for Gutkin16
        
    Arguments
    ---------
    csv: csv array
        Returns the list of csvs for every file associated to this model in the tables directory
    
    abund: str/float array
        Returns the list of abundances available for each model. 
            -If the abundances are floats, you can access the corresponding models with the .abund(value) attribute, 
            which returns a class containing the model
            -If the abundances are str, you can directly access the the models with the .yourabundancestrings attributes
        
    see shock_model_single for a description of the arguments for a single model
    """
    
    #tables directory
    currdir=os.getcwd()
    os.chdir(shock_database_path)
    
    #getting all the model files corresponding to this model name
    if 'Allen' in modelname or '2008' in modelname and 'cut' not in modelname:
        modelname='Allen2008'
        modelfiles=glob.glob(modelname+"_*.php")
    elif 'cut' in modelname:
        modelfiles=glob.glob("*cut*.php")
    else:
        modelfiles=glob.glob(modelname+"*.php")
    
    #and storing the csvs and model parameters in an array
    csv_list=np.array([None]*len(modelfiles))
    abund_list=np.array([None]*len(modelfiles))
    abund_strlist=np.array([None]*len(modelfiles))
    
    for i in range(len(modelfiles)):
        
        filename=modelfiles[i]
        csv_list[i]=pd.read_csv(filename,sep='\t')
        
        #searching for the abundances
        
        if 'Alarie' in modelname:
            #for the PNe model the only abundance is solar
            abund_list[i]='Solar'
        
        elif 'Gutkin' in modelname:
            abund_strlist[i]=filename[filename.find('d')+1:filename.find('C')-1]

        else:
            #for Allen and allen cut
            abund_list[i]=filename[filename.find('_')+1:filename.find('.php')]

     
    #for the Gutkin model we can get rid of the abundances repetitions using floats for the abundances
    if 'Gutkin' in modelname:
        #float abundances
        abund_float=np.core.defchararray.add('0.',abund_strlist)
        abund_float=np.sort(np.unique(abund_float.astype(float)))
        #string abundances used for other purposes
        abund_list=abund_float.astype(str).tolist()
        abund_list=np.array([elem.replace('.','_') for elem in abund_list])

    class shock_model_single:
        """
       subclass for a single model for a specific abundance (and dust ratio if it is variable)
       Parameters
       ----------
       csv: 
           Gives the csv to use
       dens,mag_fld,shck_vel:
           shortcuts for the lists of values
       grids:
           grid: 1-dimensionnal grid (the only remaining variable is the shock velocity 
                                            AND time/distance for the cut and Alarie models)
           lgrid gives logarithmic values for the ratios
           
           d_grid/d_lgrid: 2-dimensionnal grids (remaining variables : mag_fld/shck_vel
                                                 AND time/distance for the cut and Alarie models)
       
       Those arguments are used in the plotting functions above to superpose grids with bpts
       """
        
        def __init__(self, csv):
            if len(csv)==1:
                self.csv=csv[0]
            else:
                self.csv=csv
            self.dens=np.unique(self.csv['dens'])
            self.mag_fld=np.unique(self.csv['mag_fld'])
            self.shck_vel=np.unique(self.csv['shck_vel'])
            
        #one dimensionnal grid (lines) functions
        def grid(self, dens,mag_fld):
            return self.csv[self.csv.dens==dens][self.csv.mag_fld==mag_fld][self.csv.columns[2:]]
    
        def lgrid(self, dens,mag_fld):
            #we keep the velocity in linear and check for the presnece of infs. to transform into nans.
            log_csv=self.csv[self.csv.dens==dens][self.csv.mag_fld==mag_fld][self.csv.columns[2:]]
            log_csv=log_csv.copy()
            log_csv=log_csv.replace([np.inf,-np.inf],np.nan)
            log_csv[log_csv.columns[-4:]]=np.log10(log_csv[log_csv.columns[-4:]])
            return log_csv
            
        #two dimensionnal surfaces grids functions
        def d_grid(self,dens):
            return self.csv[self.csv.dens==dens][self.csv.columns[1:]],

        def d_lgrid(self,dens):
            #we keep the velocity in linear and check for the presnece of infs. to transform into nans.
            log_csv=self.csv[self.csv.dens==dens][self.csv.columns[1:]]
            log_csv=log_csv.copy()
            log_csv=log_csv.replace([np.inf,-np.inf],np.nan)
            log_csv[log_csv.columns[-4:]]=np.log10(log_csv[log_csv.columns[-4:]])
            return log_csv
        
    class shock_mod:
        
        """
        Global model class associated with the function shock_model
        contains information on the entire model and enables access to the more specific elements
        """
        
        def __init__(self):
            
            #list of the csvs
            self.csv=csv_list
            
            #list of STRING abundances
            self.abund=abund_list 
            if 'Gutkin' in modelname:
                #list of float abundances specifically for the Gutkin models
                self.abund_float=abund_float
                for i in range(len(abund_list)):
                    abund_var=abund_float[i]
                    
                    #dynamically adding attributes for the float abundances values, they can't be called except with getattr() 
                    #We'll rather directly abund_single here but we still create it just in case
                    setattr(self,str(abund_list[i]),self.abund_single(abund_var))
            else:
                for i in range(len(abund_list)):
                    #dynamically adding attributes for the str abundance values
                    setattr(self,abund_list[i],self.abund_single(abund_list[i]))
                    
        def abund_single(self,abundance):
            
            """
            Subwrapper for carbon to dust proportion specification for Gutkin
            Parameters
            ----------
            
            abund: str
                Set of abundances corresponding to the model
            """
            
            class shock_model_with_abundances:
                def __init__(self):
                    #dichotomy since Gutkin models will have float abundances
                    if type(abundance) in [float,np.float64,int,np.int64]:
                        #here we recheck in the file list to get the csvs corresponding to the different carbon to dust fraction
                        abund_allfloat=np.core.defchararray.add('0.',abund_strlist)
                        abund_allfloat=abund_allfloat.astype(float)
                        positions=np.where(abund_allfloat==abundance)[0]
                        
                        #now we can determine who is where and attribute the corresponding attributes
                        for i in positions:
                            if '0d26' in modelfiles[i]:
                                self.ctd_low=shock_model_single(csv_list[i])
                            if '1d00' in modelfiles[i]:
                                self.ctd_high=shock_model_single(csv_list[i])
                                
                    elif type(abundance) in [str,np.str_]:
                        position=np.where(abund_list==abundance)
                        #here since we have one less "level" of classes due to the absence of carbon to dust ratio variability
                        #we directly replace the class parameters by the single shock model class
                        self.__dict__=shock_model_single(csv_list[position]).__dict__
                        self.__class__=shock_model_single(csv_list[position]).__class__
                        
            return shock_model_with_abundances()
    
    #returning to the initial directory
    os.chdir(currdir)
    return shock_mod()
