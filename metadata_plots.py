#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go


# Load and concat datasets
luad = pd.read_csv('luad.csv', index_col=0)
luad['cancer_type'] = 'luad'

lusc = pd.read_csv('lusc.csv', index_col=0)
lusc['cancer_type'] = 'lusc'

data = pd.concat([luad, lusc])

data.age = ['wiek >= 65 lat' if i == 1 else 'wiek <65 lat' for i in data.age]
data.sex = ['mężczyzna' if i == 1 else 'kobieta' for i in data.sex]

# Translate english statements to polish
smoking_translator = {
    'light_smoker': 'palacz sporadyczny', #  (<10 paczek papierosów rocznie)
    'heavy_smoker': 'palacz nałogowy', # (>=10 paczek papierosów rocznie)
    'non_smoker': 'niepalący'
}

disease_stage_translator = {
    'early': 'wczesne', # stage I
    'middle': 'lokalnie-zaawans.', # stage II + IIIa
    'late': 'zaawansowane' # stage IIIb + IV
}

data.smoking_status = [smoking_translator.get(i) for i in data.smoking_status]
data.disease_stage = [disease_stage_translator.get(i) for i in data.disease_stage]


# Gender distribution

fig = go.Figure()
current_data_luad = data.loc[data.cancer_type.str.fullmatch('luad'), 'sex'].value_counts()
fig.add_trace(go.Bar(
    x=current_data_luad.index,
    y=current_data_luad.values.tolist(),
    name='LUAD',
    marker_color='#330C73',
    opacity=0.75,
    text=current_data_luad.values.tolist(),
    textposition='auto'
))

current_data_lusc = data.loc[data.cancer_type.str.fullmatch('lusc'), 'sex'].value_counts()
fig.add_trace(go.Bar(
    x=current_data_lusc.index,
    y=current_data_lusc.values.tolist(),
    name='LUSC',
    marker_color='#EB89B5',
    opacity=0.75,
    text=current_data_lusc.values.tolist(),
    textposition='auto'
))

fig.update_layout(
    title_text='Rozkład płci w zależności od typu raka',
    xaxis_title_text='Płeć',
    yaxis_title_text='Liczność',
    bargap=0.2,
    bargroupgap=0.1,
    font=dict(
        size=20,
    )
)
fig.show()

# Age distribution

fig = go.Figure()
current_data_luad = data.loc[data.cancer_type.str.fullmatch('luad'), 'age'].value_counts()
fig.add_trace(go.Bar(
    x=current_data_luad.index,
    y=current_data_luad.values.tolist(),
    name='LUAD',
    marker_color='#330C73',
    opacity=0.75,
    text=current_data_luad.values.tolist(),
    textposition='auto'
))

current_data_lusc = data.loc[data.cancer_type.str.fullmatch('lusc'), 'age'].value_counts()
fig.add_trace(go.Bar(
    x=current_data_lusc.index,
    y=current_data_lusc.values.tolist(),
    name='LUSC',
    marker_color='#EB89B5',
    opacity=0.75,
    text=current_data_lusc.values.tolist(),
    textposition='auto'
))

fig.update_layout(
    title_text='Rozkład wieku w zależności od typu raka',
    xaxis_title_text='Kategoria wiekowa',
    yaxis_title_text='Liczność',
    bargap=0.2,
    bargroupgap=0.1,
    font=dict(
        size=20,
    )
)
fig.show()

# Smoking status

fig = go.Figure()

current_data_luad = data.loc[data.cancer_type.str.fullmatch('luad'), 'smoking_status'].value_counts()
fig.add_trace(go.Bar(
    x=current_data_luad.index,
    y=current_data_luad.values.tolist(),
    name='LUAD',
    marker_color='#330C73',
    opacity=0.75,
    text=current_data_luad.values.tolist(),
    textposition='auto'
))

current_data_lusc = data.loc[data.cancer_type.str.fullmatch('lusc'), 'smoking_status'].value_counts()
fig.add_trace(go.Bar(
    x=current_data_lusc.index,
    y=current_data_lusc.values.tolist(),
    name='LUSC',
    marker_color='#EB89B5',
    opacity=0.75,
    text=current_data_lusc.values.tolist(),
    textposition='auto'
))

fig.update_layout(
    title_text='Rozkład profili palących w zależności od typu raka', # title of plot
    xaxis_title_text='Profile palących', # xaxis label
    yaxis_title_text='Liczność', # yaxis label
    bargap=0.2, # gap between bars of adjacent location coordinates
    bargroupgap=0.1, # gap between bars of the same location coordinates
    font=dict(
        size=20,
    )
)
fig.show()

# Smoking history

fig = go.Figure()

current_data_luad = data.loc[data.cancer_type.str.fullmatch('luad'), 'smoking_status'].value_counts()
fig.add_trace(go.Bar(
    x=current_data_luad.index,
    y=current_data_luad.values.tolist(),
    name='LUAD',
    marker_color='#330C73',
    opacity=0.75,
    text=current_data_luad.values.tolist(),
    textposition='auto'
))

current_data_lusc = data.loc[data.cancer_type.str.fullmatch('lusc'), 'smoking_status'].value_counts()
fig.add_trace(go.Bar(
    x=current_data_lusc.index,
    y=current_data_lusc.values.tolist(),
    name='LUSC',
    marker_color='#EB89B5',
    opacity=0.75,
    text=current_data_lusc.values.tolist(),
    textposition='auto'
))

fig.update_layout(
    title_text='Rozkład profili palących w zależności od typu raka',
    xaxis_title_text='Profile palących',
    yaxis_title_text='Liczność',
    bargap=0.2,
    bargroupgap=0.1,
    font=dict(
        size=20,
    )
)
fig.show()

# Disease stage

fig = go.Figure()

current_data_luad = data.loc[data.cancer_type.str.fullmatch('luad'), 'disease_stage'].value_counts()
fig.add_trace(go.Bar(
    x=current_data_luad.index,
    y=current_data_luad.values.tolist(),
    name='LUAD',
    marker_color='#330C73',
    opacity=0.75,
    text=current_data_luad.values.tolist(),
    textposition='auto'
))

current_data_lusc = data.loc[data.cancer_type.str.fullmatch('lusc'), 'disease_stage'].value_counts()
fig.add_trace(go.Bar(
    x=current_data_lusc.index,
    y=current_data_lusc.values.tolist(),
    name='LUSC',
    marker_color='#EB89B5',
    opacity=0.75,
    text=current_data_lusc.values.tolist(),
    textposition='auto'
))

fig.update_layout(
    title_text='Rozkład zaawansowania w zależności od typu raka',
    xaxis_title_text='Zaawansowanie raka',
    yaxis_title_text='Liczność',
    bargap=0.2,
    bargroupgap=0.1,
    font=dict(
        size=20,
    )
)

fig.show()


# Mutations

mutations_all = ['AKT1', 'ALK', 'BRAF', 'CD274', 'CDKN2A', 'DDR2', 'EGFR', 
                 'ERBB2', 'FGFR1', 'FGFR3', 'KEAP1', 'KRAS', 'MAP2K1', 'MET',
                 'NF1', 'NFE2L2', 'NRAS', 'NTRK1', 'PIK3CA', 'PTEN', 'RB1',
                 'RET', 'ROS1', 'SMARCA4', 'STK11', 'TP53', 'U2AF1L4']

mutations = ['TP53', 'KRAS', 'EGFR', 'NFE2L2', 'CDKN2A', 'RET', 'STK11', 'RB1']

fig = go.Figure()

current_data_luad = data.loc[data.cancer_type.str.fullmatch('luad'), mutations].sum()
fig.add_trace(go.Bar(
    x=current_data_luad.index,
    y=current_data_luad.values.tolist(),
    name='LUAD',
    marker_color='#330C73',
    opacity=0.75,
    text=current_data_luad.values.tolist(),
    textposition='auto'
))

current_data_lusc = data.loc[data.cancer_type.str.fullmatch('lusc'), mutations].sum()
fig.add_trace(go.Bar(
    x=current_data_lusc.index,
    y=current_data_lusc.values.tolist(),
    name='LUSC',
    marker_color='#EB89B5',
    opacity=0.75,
    text=current_data_lusc.values.tolist(),
    textposition='auto'
))

fig.update_layout(
    title_text='Liczność mutacji w zależności od typu raka',
    xaxis_title_text='Mutacje',
    yaxis_title_text='Liczność',
    bargap=0.2,
    bargroupgap=0.1,
    font=dict(
        size=20,
    )
)
fig.show()


# Mutations occurrences table
current_data_luad = data.loc[data.cancer_type.str.fullmatch('luad'), mutations_all].sum()
current_data_lusc = data.loc[data.cancer_type.str.fullmatch('lusc'), mutations_all].sum()

mutations_counts = pd.DataFrame(current_data_luad, columns=['luad_count']).join(
    pd.DataFrame(current_data_lusc, columns=['lusc_count']))
mutations_counts['total'] = mutations_counts.luad_count + mutations_counts.lusc_count

mutations_counts.sort_values('total', ascending=False, inplace=True)
print(mutations_counts.to_latex())

# transform english labels to polish
cindex = pd.read_csv('c-index.tsv', sep='\t')
cindex.fillna("brak", inplace=True)
cindex.drop_duplicates(inplace=True)

dataset_translator = {
    'p': 's', # pathologic stage -> statdium
    's': 'p', # smoking status -> palenia
    'ps': 'sp'
}
cindex['dataset'] = [dataset_translator.get(i, i) for i in cindex['dataset']]
cindex.columns = ['typ raka', 'statystyki tkanek', 'mean_cindex', 'sd_cindex', 'best_cindex',
                  'zestaw mutacji', 'wariant danych klinicznych']
tissues_subsets = {
    'clinical + mutations + prevalence (TUMOR, VESSEL, IMMUNE, NECROSIS, STROMA, MIXED)': 'zliczenia tkanek II',
    'clinical + mutations + prevalence (TUMOR, VESSEL, IMMUNE, NECROSIS, STROMA)': 'zliczenia tkanek I',
    'clinical + mutations + microenvironment (VESSEL, IMMUNE, NECROSIS, STROMA, BRONCHI)': 'mikrośrodowisko I',
    'clinical + mutations + microenvironment (VESSEL, IMMUNE, NECROSIS, STROMA, LUNG)': 'mikrośrodowisko II'
}

cindex['statystyki tkanek'] = [tissues_subsets.get(i, i) for i in cindex['statystyki tkanek']]
cindex['typ raka'] = cindex['typ raka'].str.upper()
cindex['wariant danych klinicznych'] = [f" wariant {i if i != 'baza' else 'bazowy'}"
                                       for i in cindex['wariant danych klinicznych']]
cindex['statystyki tkanek'] = [i.replace('clinical + mutations + ', "") for i in cindex['statystyki tkanek']]
tissues_stats_translator = {
    'ITLR': 'indeks ITLR',
    'Morisita_Horn': 'indeks Morisita-Horn',
    'Lymphocyte_Ratio': 'współczynnik limfocytów (LR)',
    'Shannon': 'indeks Shannon',
    'Simpson': 'indeks Simpson',
    'prevalence': 'zliczenia tkanek',
    'microenvironment': 'mikrośrodowisko'
}
cindex['statystyki tkanek'] = [i.split(' + ')[-1] for i in cindex['statystyki tkanek']]
cindex['statystyki tkanek'] = ['brak' if i == 'mutations' or i == 'clinical' else i for i in cindex['statystyki tkanek']]
cindex['statystyki tkanek'] = [tissues_stats_translator.get(i, i)  for i in cindex['statystyki tkanek']]
cindex['statystyki tkanek'] = cindex['statystyki tkanek'].str.replace('VESSEL', 'naczynia krwionośne')
cindex['statystyki tkanek'] = cindex['statystyki tkanek'].str.replace('NECROSIS', 'martwica')
cindex['statystyki tkanek'] = cindex['statystyki tkanek'].str.replace('IMMUNE', 'tkanka immunologiczna')
cindex['statystyki tkanek'] = cindex['statystyki tkanek'].str.replace('TUMOR', 'tkanka nowotworowa')
cindex['statystyki tkanek'] = cindex['statystyki tkanek'].str.replace('MIXED', 'tkanka mieszana')
cindex['statystyki tkanek'] = cindex['statystyki tkanek'].str.replace('BRONCHI', 'oskrzela')
cindex['statystyki tkanek'] = cindex['statystyki tkanek'].str.replace('STROMA', 'tkanka zrębu')
cindex.sort_values(['mean_cindex', 'sd_cindex'], ascending=[False, True], inplace=True)
cindex['cindex'] = [f"{row['mean_cindex']}±{row['sd_cindex']}" for _, row in cindex.iterrows()]

# adjust columns order
luad = cindex.loc[cindex['typ raka'].str.fullmatch('LUAD'),
                  ['wariant danych klinicznych', 'zestaw mutacji', 'statystyki tkanek', 'cindex',]]
luad.columns = ['wariant danych klinicznych', 'zestaw mutacji', 'statystyki tkanek', 'c-indeks']
luad.index = range(luad.shape[0])

lusc = cindex.loc[cindex['typ raka'].str.fullmatch('LUSC'), 
                  ['wariant danych klinicznych', 'zestaw mutacji', 'statystyki tkanek', 'cindex',]]
lusc.columns = ['wariant danych klinicznych', 'zestaw mutacji', 'statystyki tkanek', 'c-indeks']
lusc.index = range(lusc.shape[0])

# convert padnas dataframe to LaTeX notation
luad.to_latex('luad_cindex_table')
lusc.to_latex('lusc_cindex_table.tex')
