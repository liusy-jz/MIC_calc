import numpy as np
import pandas as pd
from minepy import MINE

# read water content class classification results
salt_phe=pd.read_csv('phe_heatmap.data.csv')
salt_label=pd.DataFrame(salt_phe.loc[:,'label'])
salt_label.index=salt_phe.loc[:,'Taxa'].values

# read metabolome data 
mt=pd.read_csv('Met-QI-all_compounds.csv')
mt.index=mt.Compound
mt=mt.drop('Compound',axis=1)
mt=mt.T
mt.columns.name=None
con=['salt']*(mt.shape[0]/2)+['water']*(mt.shape[0]/2)
mt.index=pd.MultiIndex.from_tuples([(y,x.split('.')[0])for x,y in zip(mt.index,con)])

# extract metabolites under salt condition
salt=mt.loc['salt',:]
salt=pd.merge(salt,salt_label.loc[salt_label.label!=1,:],how='inner',left_index=True,right_index=True)

# extract metabolites under ck condition
water=mt.loc['water',:]
water=pd.merge(water,salt_label.loc[salt_label.label!=1,:],how='inner',left_index=True,right_index=True)

# MINE initialization
mine=MINE()

# metabolite of salt condition MIC calculate
salt['label']=salt['label'].replace(2,1)
salt_X = salt.values[:,:-1]
salt_y = salt['label'].values

mic_scores=list()
for i in range(salt_X.shape[1]):
    mine.compute_score(salt_X[:, i], salt_y)
    mic_scores.append(mine.mic())

salt.iloc[:,np.append(np.argsort(mic_scores)[::-1][:200],[-1])].to_csv('MIC_top200_salt.csv')

# metabolite of ck condition MIC calculate
water['label']=water['label'].replace(2,1)
water_X = water.values[:,:-1]
water_y = water['label'].values

mic_scores=list()
for i in range(water_X.shape[1]):
    mine.compute_score(water_X[:, i], water_y)
    mic_scores.append(mine.mic())

water.iloc[:,np.append(np.argsort(mic_scores)[::-1][:200],[-1])].to_csv('MIC_top200_water.csv')
