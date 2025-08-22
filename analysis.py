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

def all_heavyLight_genes(dfs, gene_types, use_index=False): 
    '''
    input: 
        - dfs: dfs of all files to consider in experiment, can be dictionary or list. 
        - gene_types: list of gene types to be considered (MUST BE SAME AS NAME IN FILE)
        - use_index: when true rather than getting genes from columns the index is used for the genes. 
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
        if use_index: 
            all_genes = df.index 
        else: 
            needed_cols = [col for col in df.columns if any(col.endswith(gene) for gene in gene_types)]
            all_genes = pd.concat([df[gene] for gene in needed_cols], ignore_index=True).dropna().drop_duplicates()
        for gene in all_genes: 
            if not isinstance(gene, str):
                gene = str(gene)
            if gene.startswith('IGH'):
                all_igh.append(gene)
            elif (gene.startswith('IGK') or gene.startswith('IGL')):
                all_iglk.append(gene)
    all_igh = list(set(all_igh))
    all_iglk = list(set(all_iglk))
    return all_igh, all_iglk

def lightChain_vs_heavyRcount_relFrequency(grouped_df, light_gene_names, ighR_colName, light_gene_types): 
    '''
    Inputs: 
        grouped_df: df to get the light vs heavy chain R count relative frequency of. This df must be in the format of the outputs of formating.make_combined_rows(). 
        light_gene_names: list of all of the light gene names to group by. 
        ighR_colName: the EXACT name of the column with heavy chain R count in grouped_df. 
        light_gene_types: types of light genes to consider. 
    Outputs: 
        output_df: a new dataframe where there is a row for each light chain gene and the following columns: light_gene, Rcount_equalTo_0, Rcount_gt_0	total, R_equal_0_freq, R_gt_0_freq
    '''
    #sorting out data needed 
    wanted_cols = [f"{chain}_{col}" for chain in ['IGK', 'IGL'] for col in light_gene_types]
    wanted_cols.append(ighR_colName)
    heavyR_df = grouped_df[wanted_cols].copy()
    wanted_cols.remove(ighR_colName)
    heavyR_df = heavyR_df.dropna(subset=[ighR_colName]) #removing columns with no heavy R count (no heavy chain)
    #reorganizing data so that it is indexed by light genes 
    heavyR_df['Rcount_type'] = heavyR_df[ighR_colName].apply(
        lambda x: 'Rcount_gt_0' if x > 0 else 'Rcount_equalTo_0'
    )
    heavyR_df = heavyR_df.melt(id_vars=[ighR_colName, 'Rcount_type'], value_vars=wanted_cols, value_name='light_gene').dropna(subset=['light_gene'])
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

#but this should be for combined heavy chains right....
def lightChain_Hpair_pKa(df, all_iglk, pka_col_name):
    #adding column for combined_cdr3 heavy genes 
    formating.combine_cdr3_genes(df, ['IGH'])
    wanted_col = ['IGK_v_gene', 'IGK_j_gene', 'IGL_v_gene', 'IGL_j_gene', 'IGH_cdr3_genes', 'IGH_cdr3_pKa']
    wanted_df = df[wanted_col]
    #TODO: COME BACK AND FINISH THISSSS
    

def relative_abundance_df(df, gene_types, wanted_gene_prefix, all_wanted_genes):
    '''
    Inputs: 
        df: dataframe to get the counts of gene apperances in. 
        gene_types: types of genes to count. ex: ['v_gene', 'j_gene', 'd_gene'] (MUST MATCH column name endings in df)
        wanted_gene_prefix: prefix of the genes to count. ex: ['IGH'], ['IGK', 'IGL']
        all_wanted_genes: list of all heavy or light chain genes in the experiment (not just file) depending on wanted_gene_prefix value. 
    Output: 
        output_df: new dataframe where in row is a different gene and has the columns: count and rel_freq. The index is the gene names so can search for a specific row by gene name. 
    '''
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

def rel_charge_df(df, col_suffixes, chain_types):
    ''' Function to analyze genes of interest and their pKa/R count
    Input:
        col_suffixes - list of the endings of the columns you want information for on each chain. 
            ex: ['v_gene', 'd_gene', 'j_gene', 'cdr3_pKa', 'R_count']
        chain_types - list of the starts of column names you want. usually the chain types. 
            ex: ['IGH'] or ['IGL', 'IGK']
    Ouputs: 
        output_df: new df with the columns TYPE_cdr3_genes, TYPE_cdr3_pKa_total, TYPE_R_count_total	count, rel_freq, avg_pKa, avg_R_count, pKa_rel_freq, R_count_rel_freq. 
    '''
    if len(chain_types) > 2: 
        print("Function can only accept the following chain type inputs: ['IGH'], ['IGL'], ['IGK'], ['IGL', 'IGK']")
    #making smaller df with only needed information 
    wanted_cols = [f"{chain}_{col}" for chain in chain_types for col in col_suffixes]
    wanted_df = df[wanted_cols]
    #replacing gene columns with one cdr3 column where genes for a cells chain are one string. 
    formating.combine_cdr3_genes(wanted_df, chain_types)
    needed_ends = ['cdr3_genes']
    for end in col_suffixes: 
        if not end.endswith('gene'): needed_ends.append(end)
    needed_cols = [f"{chain}_{end}" for chain in chain_types for end in needed_ends]
    charges_df = wanted_df[needed_cols].copy() #to prevent SettingWithCopyWarning from dropna(inplace=True)
    cdr3_names = [f"{chain}_cdr3_genes" for chain in chain_types]
    prefix = chain_types[0]
    if len(chain_types) == 2: #if want to analyze all light chains need to reformat df so it still has only 3 columns 
        type1_cols = [f"{chain_types[0]}_{end}" for end in needed_ends]
        type1_df = charges_df[type1_cols].rename(columns={col: f"light_{col[4:]}" for col in type1_cols})
        type2_cols = [f"{chain_types[1]}_{end}" for end in needed_ends]
        type2_df = charges_df.drop(type1_cols, axis=1).rename(columns={col: f"light_{col[4:]}" for col in type2_cols})
        charges_df = pd.concat([type1_df, type2_df], ignore_index=True)
        prefix = 'light'
        cdr3_names = 'light_cdr3_genes'
    else: cdr3_names = cdr3_names[0]
    charges_df.dropna(subset=[f"{prefix}_cdr3_pKa"], inplace=True)
    #preparing to make output df 
    cdr3_genes = []
    cdr3_pKa_total = []
    R_count_total = []
    grouped = charges_df.groupby(by=cdr3_names)
    #getting relative abundance to add to output df 
    combined_igh, combined_iglk = all_heavyLight_genes([charges_df], cdr3_names)
    rel_abund_df = relative_abundance_df(charges_df, ['cdr3_genes'], chain_types, combined_igh if 'IGH' in chain_types else combined_iglk)
    for gene_names in rel_abund_df.index.astype(str): #want df to be in the same format so merging the dfs is easy 
        cdr3_genes.append(gene_names)
        cdr3_pKa_total.append(grouped.get_group(gene_names)[f"{prefix}_cdr3_pKa"].sum())
        R_count_total.append(grouped.get_group(gene_names)[f"{prefix}_R_count"].sum())
    total_charge_df = pd.DataFrame({f"{prefix}_cdr3_genes":cdr3_genes, f"{prefix}_cdr3_pKa_total":cdr3_pKa_total, f"{prefix}_R_count_total":R_count_total})
    output_df = pd.merge(total_charge_df, rel_abund_df.reset_index().set_axis([f"{prefix}_cdr3_genes", 'count', 'rel_freq'], axis=1),on=[f"{prefix}_cdr3_genes"], how='outer')
    output_df['avg_pKa'] = output_df[f"{prefix}_cdr3_pKa_total"]/output_df['count']
    output_df['avg_R_count'] = output_df[f"{prefix}_R_count_total"]/output_df['count']
    output_df['pKa_rel_freq'] = output_df['avg_pKa']*output_df['rel_freq']
    output_df['R_count_rel_freq'] = output_df['avg_R_count']*output_df['rel_freq']
    return output_df

#function for t-test USE FOR COMPAIRING NORMALIZED ABUNDANCE OF LIGHT CHAINS BETWEEN HEALTHY AND SLE
def t_test(healthy_list, sick_list, all_genes_considered, rel_freq_col_name):
    '''
    Inputs: 
        healthy_list: list of all the dataframes that contain information about healthy subject files with relative frequency columns. 
                    for example outputs of analysis.rel_charge_df() or analysis.relative_abundance_df(). 
        sick_list: same has healthy_list but for the files of sick subjects. 
        all_genes_considered: list of all genes to consider. usualy all_igh or all_iglk but also can be list of combined_cdr3 genes for heavy or light chains depends on what data you are looking at. 
        rel_freq_col_name: name of the column of the realitivy frequency of the thing to be analyzed. Must match name in dataframes.
    Output: 
        new dataframe with a row for each gene in all_genes_considered and the following columns: Genes, Equal Population Variance t_stats, and Equal Population Variance p_values
    '''
    equal_pop_var_pvals = []
    equal_pop_var_tstats = []
    for gene in all_genes_considered:
        hcds = [healthy_list[i].loc[gene][rel_freq_col_name] for i in range(len(healthy_list))]
        sles = [sick_list[i].loc[gene][rel_freq_col_name] for i in range(len(sick_list))]
        (et_stat, ep_val) = ttest_ind(hcds, sles)
        equal_pop_var_pvals.append(ep_val)
        equal_pop_var_tstats.append(et_stat)
    return pd.DataFrame({"Gene": all_genes_considered, "Equal Population Variance t_stats":equal_pop_var_tstats, "Equal Population Variance p_values":equal_pop_var_pvals})

#t-test USE FOR PAIRS OF LIGHT AND HEAVY CHAINS
def normalize_pairMatrix(pairMatrix_df):
    '''Function to take a given inputed matrix shaped df and output a normalized matrix.'''
    return pairMatrix_df.div(pairMatrix_df.sum(axis=1), axis=0)

def t_test_matrix(healthyPairM_list, sickPairM_list, all_igh, all_iglk):
    '''
    Does everything like t_test except for matrix inputs and outputs a matrix with the t_test values for each heavy and light gene pair. 
    '''
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
    '''
    Inputs: 
        iglk_hRcount_df: dataframe with two different relative populations for each grouping. 
                        For example for each light chain gene the relative freqency or pairs heavy chains that have Rcounts == 0 and Rcounts > 0. 
        sucess_colName: EXACT name of the success column in the dataframe. For example: "Rcount_gt_0"
        total_colName: EXACT name of the column in dataframe which has the total number of observations in any group. 
                      For example: "total" which the the column that has the total number of times a specific light chain is observed. 
        sucessP_colName: EXACT name of the column with the relative frequency of the sucess column. For example: "R_gt_0_freq"
    Ouput: 
        new dataframe with an additional column for the Binomial Distribution survival function value for the information inputed. 
    '''
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