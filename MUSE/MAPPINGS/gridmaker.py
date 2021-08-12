# @Author: Maxime Parra
# @Date:   29-06-2021
# @Email:  maxime.parra@irap.omp.eu
# @Last modified by:   mparra
# @Last modified time: 19-07-2021

'''
This script downloads MAPPINGSV shock only models and stores them in the current directory

Database paper : https://doi.org/10.22201/ia.01851101p.2019.55.02.21

The database is on SQL i.e. blocked on IRAP network

The resulting csv files contain the LINEAR BPT ratios (BPT uses log10)

pymysql needs to be installed

The model names and abundance are manually inputed in the code.
As such, it can be wise to check if they are new models/parameters on the dtb, at :
http://3mdb.astro.unam.mx:3686/paramexplorer (blocked on IRAP network due to SQL port)

model parameter units :

density                  cm⁻³
magnetic field           microGauss
shock velocity           km.s⁻¹
time                     s
distance                 cm
'''


import os
import numpy as np
import pymysql
import pandas as pd

os.system('mkdir -p tables')
os.chdir('tables')

''' DATABASE CONNECTION'''

# Get credentials stored from the environment variables
host = os.environ['MdB_HOST']
user = os.environ['MdB_USER']
passwd = os.environ['MdB_PASSWD']
port = os.environ['MdB_PORT']
                                
# Connect to the database
db = pymysql.connect(host=host, user=user, passwd=passwd, port=int(port), db='3MdBs')

''' QUERIES AND SAVING'''

#models list
models=np.array(['Allen08','Gutkin16','Alarie19s','Allen08-cut'])

#abundances list (2D array with abundances list for each model)
abunds=np.array([None]*4)
abunds[0]=np.core.defchararray.add('Allen2008_',np.array(['Dopita2005','LMC','SMC','Solar','TwiceSolar']))
abunds[2]=np.array(['3MdB-PNe2014-solar'])
abunds[3]=abunds[0]

abund_gutkin=[]
for id_1 in ['0001','0002','0005','001','002','004','006','008','010','014','01524','017','01','020','02','030','03','040','04']:
    if id_1 not in ['010','01','020','02','030','03','040','04']:
        for id_2 in ['0d26','1d00']:
            abund_gutkin.append('Gutkin16_ISM0d'+id_1+'_C'+id_2)
    else:
        if len(id_1)==2:
            abund_gutkin.append('Gutkin16_ISM0d'+id_1+'_C'+'0d26')
        elif len(id_1)==3:
            abund_gutkin.append('Gutkin16_ISM0d'+id_1+'_C'+'1d00')
            
abunds[1]=np.array(abund_gutkin)
            
abunds_str=abunds.copy()
abunds_str[3]=np.core.defchararray.add('Allen2008-cut_',np.array(['Dopita2005','LMC','SMC','Solar','TwiceSolar']))
abunds_str[2]=np.array(['Alarie19s_3MdB-PNe2014-solar'])

model_types=['shock','precursor','shock_plus_precursor']

for k in range(3):
    os.mkdir(model_types[k])
    os.chdir(model_types[k])
    for i in range(len(models)):
        if i<2:
            for j in range(len(abunds[i])):
                
                # SQL request with Pandas
                curr_csv = pd.read_sql("""SELECT   shock_params.preshck_dens AS dens,
                                                 shock_params.mag_fld AS mag_fld,
                                                 shock_params.shck_vel AS shck_vel, 
                                                 emis_VI.OIII_5007/emis_VI.HI_4861 AS OIII_Hb,
                                                 emis_VI.NII_6583/emis_VI.HI_6563 AS NII_Ha,
                                                (emis_VI.SII_6716+emis_VI.SII_6716)/emis_VI.HI_6563 AS SII_Ha,
                                                 emis_VI.OI_6300/emis_VI.HI_6563 AS OI_Ha
                                        FROM shock_params 
                                        INNER JOIN emis_VI ON emis_VI.ModelID=shock_params.ModelID
                                        INNER JOIN abundances ON abundances.AbundID=shock_params.AbundID
                                        WHERE emis_VI.model_type='"""+model_types[k]+"""'
                                        AND shock_params.ref='"""+models[i]+"""'
                                        AND abundances.name='"""+abunds[i][j]+"""'
                                        ORDER BY preshck_dens, mag_fld, shck_vel;""", con=db)
                                        
                curr_csv.to_csv(r'./'+abunds_str[i][j]+'.php',sep='\t',mode='w+',index=False)
        else:
            for j in range(len(abunds[i])):
                
                # SQL request with Pandas
                curr_csv = pd.read_sql("""SELECT   shock_params.preshck_dens AS dens,
                                                 shock_params.mag_fld AS mag_fld,
                                                 shock_params.shck_vel AS shck_vel,
                                                 shock_params.time AS time, 
                                                 shock_params.distance AS distance,
                                                 emis_VI.OIII_5007/emis_VI.HI_4861 AS OIII_Hb,
                                                 emis_VI.NII_6583/emis_VI.HI_6563 AS NII_Ha,
                                                (emis_VI.SII_6716+emis_VI.SII_6716)/emis_VI.HI_6563 AS SII_Ha,
                                                 emis_VI.OI_6300/emis_VI.HI_6563 AS OI_Ha
                                        FROM shock_params 
                                        INNER JOIN emis_VI ON emis_VI.ModelID=shock_params.ModelID
                                        INNER JOIN abundances ON abundances.AbundID=shock_params.AbundID
                                        WHERE emis_VI.model_type='"""+model_types[k]+"""'
                                        AND shock_params.ref='"""+models[i]+"""'
                                        AND abundances.name='"""+abunds[i][j]+"""'
                                        ORDER BY preshck_dens, mag_fld, shck_vel, time, distance;""", con=db)
                                        
                curr_csv.to_csv(r'./'+abunds_str[i][j]+'.php',sep='\t',mode='w+',index=False)
    os.chdir('..')
