from collections import Counter
from scipy.stats import binom, ttest_ind
import numpy as np
import pandas as pd 
import formating 
'''WARNING: pandas dataframe has limitations to how much it can handle, computers with 8-16 GB RAM should be fine unless file has millions of rows. 
More information on limitations and possible solutions:  https://www.geeksforgeeks.org/python/how-many-rows-can-pandas-dataframe-handle/'''

def get_gene_counts_df(df, gene_types):
    '''input: 
        df - df where each row has all information for each barcode 
        gene_types - types of genes to consider (MUST have same column name as excel, ex: ['v_gene', 'j_gene', 'd_gene'])
    output: 
        gene_counts_df - pandas serise of total counts of  heavy and light chains in the dataset'''
    needed_columns = [col for col in df.columns if any(col.endswith(gene) for gene in gene_types)]
    all_genes = pd.concat([df[col] for col in needed_columns], ignore_index=True)
    gene_counts_df = all_genes.value_counts().rename_axis('gene').reset_index(name='count')
    return gene_counts_df.copy()

def get_pair_counts(df, gene_types): 
    '''
    input: 
        df - df where each row has all information for each barcode 
        gene_types - types of genes to consider (MUST have same column name as excel, ex: ['v_gene', 'j_gene', 'd_gene'])
    output: 
        pairs_df - pandas df of all of the heavy and light pairs in each cells 
    '''
    needed_columns = [col for col in df.columns if any(col.endswith(gene) for gene in gene_types)]
    pair_counts = Counter()
    for index, row in df.iterrows():
        cell_genes = row[needed_columns].values.flatten()
        #seperating heavy and light chains
        heavies = []
        lights = []
        for gene in cell_genes: 
            if isinstance(gene, str):
                if gene.startswith('IGH'):
                    heavies.append(gene)
                elif (gene.startswith('IGK') or gene.startswith('IGL')):
                    lights.append(gene)
        for h in heavies: 
            for l in lights: 
                pair_counts[(h, l)] += 1
    #converting counter to dataframe 
    pairs_df = pd.DataFrame([(k[0], k[1], v) for k, v, in pair_counts.items()], 
                            columns=["IGH", "IGKL", "Pair Count"])    
    return pairs_df

def all_heavyLight_genes(dfs, gene_types): 
    '''
    input: 
        - dfs: dfs of all files to consider in experiment, can be dictionary or list. 
        - gene_types: list of gene types to be considered (MUST BE SAME AS NAME IN FILE)
    output: 
        - all_igh: list of all the heavy chain genes in the experiment 
        - all_iglk: list of all light chain genes in the experiment 
    '''
    all_igh = []
    all_iglk = []
    df_list = dfs
    if isinstance(dfs, dict): #so function will work with dictionaries too
        df_list = list(dfs.values())
    for df in df_list:
        needed_cols = [col for col in df.columns if any(col.endswith(gene) for gene in gene_types)]
        all_genes = pd.concat([df[gene] for gene in needed_cols], ignore_index=True).dropna().drop_duplicates()
        for gene in all_genes: 
            if gene.startswith('IGH'):
                all_igh.append(gene)
            elif (gene.startswith('IGK') or gene.startswith('IGL')):
                all_iglk.append(gene)
    all_igh = list(set(all_igh))
    all_iglk = list(set(all_iglk))
    return all_igh, all_iglk

def lightChain_vs_heavyRcount_relFrequency(grouped_df, light_gene_names, ighR_colName): 
    #sorting out data needed 
    heavyR_df = grouped_df[[ighR_colName, 'IGL_v_gene', 'IGK_v_gene']].copy()
    heavyR_df = heavyR_df.dropna(subset=[ighR_colName]) #removing columns with no heavy R count (no heavy chain)
    #reorganizing data so that it is indexed by light genes 
    heavyR_df['Rcount_type'] = heavyR_df[ighR_colName].apply(
        lambda x: 'Rcount_gt_0' if x > 0 else 'Rcount_equalTo_0'
    )
    light_gene_cols = ['IGL_v_gene', 'IGK_v_gene']
    heavyR_df = heavyR_df.melt(id_vars=[ighR_colName, 'Rcount_type'], value_vars=light_gene_cols, value_name='light_gene').dropna(subset=['light_gene'])
    heavyR_df = heavyR_df.drop(columns='variable')
    output_df = (heavyR_df.groupby(['light_gene', 'Rcount_type']).size().reset_index(name='count').pivot(index='light_gene', columns='Rcount_type', values='count')
                .fillna(0).astype(int))
    for col in ['Rcount_equalTo_0', 'Rcount_gt_0']: #so program does not crash if one column has no values
        if col not in output_df.columns:
            output_df[col] = 0
    #make all light genes present 
    output_df = output_df.reindex(light_gene_names, fill_value=0)
    #output_df.index.name = None
    output_df = output_df.reset_index()
    #relative frequency analysis 
    output_df['total'] = output_df['Rcount_equalTo_0']+output_df['Rcount_gt_0']
    output_df['R_equal_0_freq'] = output_df['Rcount_equalTo_0']/output_df['total'] *100
    output_df['R_gt_0_freq'] = output_df['Rcount_gt_0']/output_df['total'] *100
    return output_df 

def relative_abundance_df(df, gene_types, wanted_gene_prefix, all_wanted_genes):
    all_gene_counts_df = get_gene_counts_df(df, gene_types)
    mask = False
    for prefix in wanted_gene_prefix:
        mask |= all_gene_counts_df['gene'].str.startswith(prefix)
    output_df = all_gene_counts_df[mask]
    missing_genes = pd.DataFrame({
        'gene': list(set(all_wanted_genes) - set(output_df['gene'].unique())), 
        'count': 0})
    output_df = pd.concat([output_df, missing_genes], ignore_index=True)
    output_df = formating.sortCols_byGenes(output_df, 6, 'gene')
    output_df['rel_freq'] = output_df['count']/output_df['count'].sum()*100
    return output_df.set_index('gene')

#function for t-test USE FOR COMPAIRING NORMALIZED ABUNDANCE OF LIGHT CHAINS BETWEEN HEALTHY AND SLE
def t_test(healthy_list, sick_list, all_genes_considered):
    equal_pop_var_pvals = []
    equal_pop_var_tstats = []
    welch_p_values = []
    welch_t_stats = []
    for gene in all_genes_considered:
        hcds = [healthy_list[i].loc[gene]['rel_freq'] for i in range(len(healthy_list))]
        sles = [sick_list[i].loc[gene]['rel_freq'] for i in range(len(sick_list))]
        (et_stat, ep_val) = ttest_ind(hcds, sles)
        equal_pop_var_pvals.append(ep_val)
        equal_pop_var_tstats.append(et_stat)
        (wt_stat, wp_val) = ttest_ind(hcds, sles, equal_var=False)
        welch_p_values.append(wp_val) #YOU HAVE 3 DATA POINTS PER ARRAY ASSUME IT IS NORMAL DISTRIBUTION DO NOT NEED THIS!
        welch_t_stats.append(wt_stat)
    return pd.DataFrame({"Gene": all_genes_considered, "Equal Population Variance t_stats":equal_pop_var_tstats, "Equal Population Variance p_values":equal_pop_var_pvals, "Welch t_stats":welch_t_stats, "Welch p_values":welch_p_values})

#t-test USE FOR PAIRS OF LIGHT AND HEAVY CHAINS
def normalize_pairMatrix(pairMatrix_df):
    return pairMatrix_df.div(pairMatrix_df.sum(axis=1), axis=0)

def t_test_matrix(healthyPairM_list, sickPairM_list, all_igh, all_iglk):
    for index in range(len(healthyPairM_list)): 
        if (healthyPairM_list[index].iloc[0].sum() > 1):
            healthyPairM_list[index] = normalize_pairMatrix(healthyPairM_list[index])
    for index in range(len(sickPairM_list)): 
        if (sickPairM_list[index].iloc[0].sum() > 1):
            sickPairM_list[index] = normalize_pairMatrix(sickPairM_list[index])
    p_values = []
    ighs = []
    iglks = []
    t_stats = []
    for iglk in all_iglk:
        for igh in all_igh: 
            hcds = [healthyPairM_list[i][igh][iglk] for i in range(len(healthyPairM_list))]
            sles = [sickPairM_list[i][igh][iglk] for i in range(len(sickPairM_list))]
            (t_stat, p_val) = ttest_ind(hcds, sles)
            t_stats.append(t_stat)
            p_values.append(p_val) 
            ighs.append(igh)
            iglks.append(iglk)
    return pd.DataFrame({"IGH": ighs, "IGKL": iglks, "t_stats":t_stats, "p_values":p_values})

#function for bionomial distribution FOR LIGHT CHAINS AND R STUFF
def binomalDist(iglk_hRcount_df, sucess_colName, total_colName, sucessP_colName):
    p = np.nanmean(iglk_hRcount_df[sucessP_colName].values)/100 #probability of success on a given trial (average of R_equal_0_freq)
    bino_dists = [] #Survival function outputs or 1-cdf 
    for index, row in iglk_hRcount_df.iterrows(): 
        n = row[total_colName] #number of trials (number of observations of a light chain)
        if (n == 0):
            bino_dists.append(np.nan)
            continue
        k = row[sucess_colName] #number of sucesses (num of light chains paired with heavy chains R_equal_0)
        bino_dists.append(binom.sf(k, n, p))
    iglk_hRcount_df["Binomial_Dist_sf"] = bino_dists
    return iglk_hRcount_df

