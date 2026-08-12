# TCAM microbiota data
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

from mprod.dimensionality_reduction import TCAM
from mprod import table2tensor
import random
from itertools import combinations
from multiprocessing import Pool
import numpy as np
import pandas as pd
import re
from sklearn.metrics import pairwise_distances
from scipy.stats import f
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy.spatial.distance import cdist
from scipy.special import betainc
from scipy.stats import f as fdist
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os

# Set Seaborn style
sns.set_style('ticks')

# Permanova functions
def pairwise_melt_meta(meta, data):

    _dm = pairwise_distances(data.loc[meta['SampleID']])
    for i in range(_dm.shape[0]):
        _dm[i,i:] = None
    _dm_df = pd.DataFrame(_dm, index = meta['SampleID'].rename('sample1'), columns=meta['SampleID'].rename('sample2'))

    _dm_melt = pd.melt(_dm_df.reset_index(), id_vars=['sample1'], value_name='d').dropna()

    metacols = ['rGroup','Participant', 'SampleID']
    _dm_melt2 = meta.loc[:, metacols].merge(_dm_melt, right_on = 'sample2', left_on = 'SampleID', how = 'right').drop('SampleID', axis=1)
    _dm_melt_meta = meta.loc[:, metacols].merge(_dm_melt2, right_on = 'sample1', left_on= 'SampleID', how = 'right', suffixes = ("1","2")).drop('SampleID', axis=1)
    return _dm_melt_meta

def random_combination(iterable, r):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n), r))
    return tuple(pool[i] for i in indices)

def ss_total(data):
    return (data['d']**2).sum()

def ss_within(data):
    if 'rGroup11' in data.columns:
        # Group11 Group22 are the permuted labels of samples in Group1, Group2 resp'
        return data.query('rGroup11 == rGroup22').groupby('rGroup11')['d'].apply(lambda x: (x**2).sum()).sum()
    else:
        return data.query('rGroup1 == rGroup2').groupby('rGroup1')['d'].apply(lambda x: (x**2).sum()).sum()

def ss_between(data, ss_t):
    return ss_t - ss_within(data)

def Fstat(data, ss_t = None, N = None):
    # p = 2 as there are always 2 groups in our cases
    if N is None:
        N = len(set(data['sample1'].tolist() + data['sample2'].tolist()))
    if ss_t is None:
        ss_t = ss_total(data)
    ss_w = ss_within(data)
    ss_b = ss_t - ss_w
    return ss_b / (ss_w / (N - 2))

def gen_rperm(meta, r):
    _mice_map = meta.loc[:, ['Participant', 'rGroup']].copy()
    _mice_map.drop_duplicates(inplace=True)
    _mice_map.index = _mice_map['Participant']
    gsizes = _mice_map.groupby('rGroup').size()
    g1size, g1label, g2label = gsizes[0], gsizes.index[0], gsizes.index[-1]

    for combo in  random_combination(combinations(_mice_map.index, g1size), r):
        _mice_map['Group_perm'] = g2label
        _mice_map.loc[combo ,'Group_perm'] = g1label
        yield _mice_map['Group_perm'].to_dict()

def _permfs(args):
    __data, tv, N, permdict = args
    _data = __data.copy()
    _data['rGroup11'] = _data['Participant1'].map(permdict)
    _data['rGroup22'] = _data['Participant2'].map(permdict)
    fstat = Fstat(_data, tv, N)
    return fstat

def run_permanova(meta, data, nperms = 1000, n_processes = 1):
    _dmm = pairwise_melt_meta(meta, data)

    _dmm['rGroup11'] = _dmm['rGroup1'].copy()
    _dmm['rGroup22'] = _dmm['rGroup2'].copy()

    N = len(set(_dmm['sample1'].tolist() + _dmm['sample2'].tolist()))
    total_var = ss_total(_dmm)
    fs_obs = Fstat(_dmm, total_var, N)
    
    if n_processes > 1:
        with Pool(processes=n_processes) as p:
            fs_perms = p.imap(_permfs, ((_dmm, total_var, N, dd) for dd in gen_rperm(meta, nperms)))
            fs_perms = np.array([x for x in fs_perms])
    else:
        fs_perms = np.array([
            _permfs((_dmm, total_var, N, dd))
            for dd in gen_rperm(meta, nperms)
        ])

    return ((fs_perms > fs_obs).sum() ) / (nperms ) 

# Data import
file_path = "new_tcam/imputed_microbiome_data_2.csv"
os.makedirs("pythonoutput_2", exist_ok=True)
data_raw = pd.read_csv(file_path, index_col=[0, 1], sep=",")

# Reset MultiIndex correctly
df1 = data_raw.reset_index()  # Converts MultiIndex into columns
df1.rename(columns={df1.columns[0]: "ID"}, inplace=True)  # Rename first index level
print(df1.head())
print(df1.columns)

# filepath: /omics/groups/OE0554/internal_temp/barbara/projects/als/scripts/2_6_tcam.py
# Sort
df2 = df1.sort_values(by=['MouseID', 'Age_ints'], ascending=[True, True])
print("Sorted DataFrame:")
print(df2.head())
print(df2.columns)

# Reshape and normalize
pseudocount = 0.1
meta_cols = ['ID', 'MouseID', 'Age_weeks', 'Age_ints', 'Genotype', 'Sex', 'GenotypePerSex']
feature_cols = [col for col in df2.columns if col not in meta_cols]
cols = ['MouseID', 'Age_ints'] + feature_cols
df3 = df2[cols]\
    .set_index(['MouseID', 'Age_ints'])\
    .groupby(level='MouseID')\
    .apply(lambda x: np.log2((x + pseudocount) / (x.loc[x.index.get_level_values('Age_ints') == 6].mean() + pseudocount)))
print("Reordered DataFrame with new index:")
print(df3.head())
print(df3.shape)
df3.index = df3.index.droplevel(0)
df3.index

# Check for missing data and display only columns with missing values
missing_data = df3.isnull().sum()
missing_data = missing_data[missing_data > 0]
print(missing_data) # no missing data

# TCAM
# Allow missing MouseID x Age_ints combinations; returned tensor will contain masked entries.
tensor_data, mode1_map, mode3_map = table2tensor(df3, missing_flag=True)
mode1_reverse_map = {val:k for k,val in mode1_map.items()}
print(tensor_data.shape)
print(mode1_map)
print(mode1_reverse_map)

tcam = TCAM(n_components=None)    # will produce minimal number of components such that
                                # total explained var > 99%

transformed_data = tcam.fit_transform(tensor_data)
df_tca = pd.DataFrame(transformed_data).rename(index = mode1_reverse_map)
rounded_expvar = np.round(100*tcam.explained_variance_ratio_, 2)
df_tca.columns = [f'F{i}:{val}%' for i,val in enumerate(rounded_expvar, start = 1)]
df_plot = df1.iloc[:,1:7].merge(df_tca, left_on = 'MouseID', right_index=True)
f1 = df_tca.columns[0]
f2 = df_tca.columns[1]

fig,ax  = plt.subplots(1,1, figsize=[4,4], dpi = 500)
sns.scatterplot(data = df_plot, x = f1, y = f2 , hue='GenotypePerSex', ax = ax, edgecolor = 'k', linewidths=.1)
ax.legend(loc='center left', fontsize = 6, fancybox = False, framealpha = 1, edgecolor = 'k', bbox_to_anchor=(1, 0.5))
plt.savefig('new_tcam/scatterplot.png', format='png', dpi=300, bbox_inches='tight')

df_plot.to_csv('new_tcam/pythonoutput_2/df_plot.csv', index=False)

# Permanova
n_factors = (np.cumsum(tcam.explained_variance_ratio_) < .2).sum() + 1  # PC for 20% variance
perm_res = pd.DataFrame(columns=['g1', 'g2', 'p'])
# Ensure meta_perm is correctly defined
meta_perm = df1.iloc[:,1:7]
meta_perm['Participant'] = meta_perm['MouseID']
print("meta_perm head:")
print(meta_perm.head())

# Check unique values in 'GenotypePerSex'
unique_values = meta_perm['GenotypePerSex'].unique()
print("Unique values in 'GenotypePerSex':")
print(unique_values)

# Check if there are enough unique values to generate combinations
for g1, g2 in combinations(unique_values, 2):
    print(g1)
    print(g2)
    meta_comp = meta_perm.loc[meta_perm['GenotypePerSex'].isin([g1, g2]) & (meta_perm['Age_ints'] == 6)].copy()
    meta_comp.rename(columns={'GenotypePerSex': 'rGroup'}, inplace=True)
    
    # Ensure 'SampleID' column is present
    if 'SampleID' not in meta_comp.columns:
        meta_comp['SampleID'] = meta_comp['MouseID']
    
    print("meta_comp head:")
    print(meta_comp.head())
    
    table_comp = df_tca.loc[meta_comp['MouseID']].copy()
    
    n_all = meta_comp.loc[meta_comp['rGroup'].isin([g1, g2])].shape[0]
    n_1 = meta_comp.loc[meta_comp['rGroup'].isin([g1])].shape[0]
    n_2 = meta_comp.loc[meta_comp['rGroup'].isin([g2])].shape[0]
    n_perm = min(1000, n_all)  # Ensure n_perm is not larger than the number of available samples
    
    # Ensure run_permanova is correctly defined and used
    try:
        p = run_permanova(meta_comp, table_comp.iloc[:, :n_factors].copy(), nperms=n_perm, n_processes=1)
        print("PERMANOVA p-value for {} vs {}: {}".format(g1, g2, p))
        perm_res.loc[perm_res.shape[0], :] = [g1, g2, p]
    except Exception as e:
        print("Error running PERMANOVA for {} vs {}: {}".format(g1, g2, e))

print(perm_res)
perm_res.to_csv("new_tcam/pythonoutput_2/perm_res.csv")

# Loadings
loadings = tcam.mode2_loadings
df_loadings = pd.DataFrame(loadings, index =  df3.iloc[:,0:].columns)
df_loadings.columns = [f'F{i}:{val}%' for i,val in enumerate(rounded_expvar, start = 1)]
df_loadings.to_csv('new_tcam/pythonoutput_2/df_loadings.csv', index=True)

# PC1
loadings_f1 = df_loadings.iloc[:,0].sort_values().copy()
top_bottom = [ loadings_f1.nlargest(10), loadings_f1.nsmallest(10) ]
loadings_f1 = pd.concat(top_bottom)
loadings_f1 = loadings_f1.sort_values()
loadings_f1 = loadings_f1.rename('F1')
loadings_f1 = loadings_f1.reset_index()

fig,ax  = plt.subplots(1,1, figsize=[1.2,2], dpi = 500)
ax.barh(y = loadings_f1.query('F1 < 0').index, width = loadings_f1.query('F1 < 0')['F1'] 
       , linewidth = .25, edgecolor = 'k')
ax.barh(y = loadings_f1.query('F1 > 0').index.astype(int), width = loadings_f1.query('F1 > 0')['F1']
        , linewidth = .25, edgecolor = 'k')
ax.set_yticks([-.7,loadings_f1.index.max() + 3])
ax.set_yticklabels([])

loadings_f1.index.names
for i, row in loadings_f1.iterrows():
    sname = str(row['index'])
    sname = sname.replace(";__", "")
    sname = sname.replace(";s__", "_")
    sname = re.sub(";s__$", "", sname)
    sname = re.sub(r'\w__[;]*','', sname)
    sname = re.sub("[;_]*$","",  sname)
    sname = sname.split(";")[-1]
    ax.text(x=0.7, y=i, s = sname, fontdict={'size':3})
    
ax.yaxis.set_tick_params(size = 3, length = 400,direction = 'inout')
sns.despine(top=True, bottom=True, right=True, trim=True)
ax.spines['left'].set_position('zero')
plt.savefig('new_tcam/pythonoutput_2/barplot_pc1.png', format='png', dpi=300, bbox_inches='tight')

# PC 2
loadings_f2 = df_loadings.iloc[:,1].sort_values().copy()
top_bottom = [ loadings_f2.nlargest(10), loadings_f2.nsmallest(10) ]
loadings_f2 = pd.concat(top_bottom)
loadings_f2 = loadings_f2.sort_values()
loadings_f2 = loadings_f2.rename('F2')
loadings_f2 = loadings_f2.reset_index()

fig,ax  = plt.subplots(1,1, figsize=[1.2,2], dpi = 500)
ax.barh(y = loadings_f2.query('F2 < 0').index, width = loadings_f2.query('F2 < 0')['F2'] 
       , linewidth = .25, edgecolor = 'k')
ax.barh(y = loadings_f2.query('F2 > 0').index.astype(int), width = loadings_f2.query('F2 > 0')['F2']
        , linewidth = .25, edgecolor = 'k')
ax.set_yticks([-.7,loadings_f2.index.max() + 3])
ax.set_yticklabels([])

loadings_f2.index.names
for i, row in loadings_f2.iterrows():
    sname = str(row['index'])
    sname = sname.replace(";__", "")
    sname = sname.replace(";s__", "_")
    sname = re.sub(";s__$", "", sname)
    sname = re.sub(r'\w__[;]*','', sname)
    sname = re.sub("[;_]*$","",  sname)
    sname = sname.split(";")[-1]
    ax.text(x=0.7, y=i, s = sname, fontdict={'size':3})
    
ax.yaxis.set_tick_params(size = 3, length = 400, direction = 'inout')
sns.despine(top=True, bottom=True, right=True, trim=True)
ax.spines['left'].set_position('zero')
plt.savefig('new_tcam/pythonoutput_2/barplot_pc2.png', format='png', dpi=300, bbox_inches='tight')