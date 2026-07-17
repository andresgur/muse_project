#!/usr/bin/env python3
"""
Download MAPPINGS V shock models from 3MdBs and save as tab-separated files.
Database: http://3mdb.astro.unam.mx:3686/
"""

import os
import pymysql
import pandas as pd


outdir = "tables"
os.makedirs(outdir, exist_ok=True)
os.chdir(outdir)

# --- CONNECTION ---
db = pymysql.connect(host='3mdb.astro.unam.mx', user='OVN_user',
                     password='oiii5007', port=3306, database='3MdBs')

# --- MODEL / ABUNDANCE DEFINITIONS ---
# Each entry: (ref name in DB, abundance name in DB, output file prefix)
models = []

# Allen08 and Gutkin16 (no time/distance columns)
allen_abunds = ['Dopita2005', 'LMC', 'SMC', 'Solar', 'TwiceSolar']
for ab in allen_abunds:
    models.append(('Allen08',   f'Allen2008_{ab}',   f'Allen2008_{ab}'))
    models.append(('Allen08',   f'Allen2008_{ab}',   f'Allen2008-cut_{ab}'))  # same query, different label

gutkin_abunds = []
for id_1 in ['0001','0002','0005','001','002','004','006','008','010','014','01524','017','01','020','02','030','03','040','04']:
    if id_1 not in ['010','01','020','02','030','03','040','04']:
        for id_2 in ['0d26','1d00']:
            gutkin_abunds.append(f'Gutkin16_ISM0d{id_1}_C{id_2}')
    else:
        id_2 = '0d26' if len(id_1) == 2 else '1d00'
        gutkin_abunds.append(f'Gutkin16_ISM0d{id_1}_C{id_2}')

for ab in gutkin_abunds:
    models.append(('Gutkin16', ab, ab))

# Alarie19s (has time/distance columns)
models_td = [('Alarie19s', '3MdB-PNe2014-solar', 'Alarie19s_3MdB-PNe2014-solar')]

# --- SQL TEMPLATES ---
SELECT_BASE = """
    SELECT shock_params.preshck_dens  AS dens,
           shock_params.mag_fld       AS mag_fld,
           shock_params.shck_vel      AS shck_vel,
           {extra_cols}
           emis_VI.OIII_5007/emis_VI.HI_4861              AS OIII_Hb,
           emis_VI.NII_6583/emis_VI.HI_6563               AS NII_Ha,
           (emis_VI.SII_6716+emis_VI.SII_6731)/emis_VI.HI_6563 AS SII_Ha,
           emis_VI.OI_6300/emis_VI.HI_6563                AS OI_Ha,
           emis_VI.HeII_4686/emis_VI.HI_4861              AS HeII_Hb,
           emis_VI.HI_4861              AS Hb
    FROM shock_params
    INNER JOIN emis_VI      ON emis_VI.ModelID   = shock_params.ModelID
    INNER JOIN abundances   ON abundances.AbundID = shock_params.AbundID
    WHERE emis_VI.model_type = '{model_type}'
      AND shock_params.ref   = '{ref}'
      AND abundances.name    = '{abund}'
    ORDER BY preshck_dens, mag_fld, shck_vel{extra_order};
"""

EXTRA_COLS_TD    = "shock_params.time AS time, shock_params.distance AS distance,"
EXTRA_ORDER_TD   = ", time, distance"

# --- DOWNLOAD LOOP ---
for model_type in ['shock', 'precursor', 'shock_plus_precursor']:
    os.makedirs(model_type, exist_ok=True)
    os.chdir(model_type)

    for ref, abund, outname in models:
        query = SELECT_BASE.format(
            extra_cols='', extra_order='',
            model_type=model_type, ref=ref, abund=abund
        )
        pd.read_sql(query, con=db).to_csv(f'{outname}.php', sep='\t', index=False)

    for ref, abund, outname in models_td:
        query = SELECT_BASE.format(
            extra_cols=EXTRA_COLS_TD, extra_order=EXTRA_ORDER_TD,
            model_type=model_type, ref=ref, abund=abund
        )
        pd.read_sql(query, con=db).to_csv(f'{outname}.php', sep='\t', index=False)

    print(f"Downloaded {model_type} models")
    os.chdir('..')

db.close()

print(f"Tables stored in {outdir}")