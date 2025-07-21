from collections import Counter
import matplotlib.pyplot as plt
from io import BytesIO
from openpyxl.drawing.image import Image
import re
import seaborn as sns
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.multitest import multipletests
from scipy.stats import hypergeom
import numpy as np
import pandas as pd 
'''WARNING: pandas dataframe has limitations to how much it can handle, computers with 8-16 GB RAM should be fine unless file has millions of rows. 
More information on limitations and possible solutions:  https://www.geeksforgeeks.org/python/how-many-rows-can-pandas-dataframe-handle/'''

''' HELPER FUNCTIONS '''
def remove_rows(df, col_name, good_val):
    return df[df[col_name] == good_val]

def remove_by_threshold(df, gene_counts_df, threshold): #only remove bottom threshold of light chains (remove whole cell row)
    if (threshold >= 1):
        print("The threshold should be a decimal value between 0 to 1 (not inclusive).")
        return
    light_count_df = gene_counts_df[(gene_counts_df['gene'].str.startswith("IGK")) | (gene_counts_df['gene'].str.startswith("IGL"))]
    l_min = light_count_df['count'].quantile(threshold)
    l_to_exc = light_count_df.query("count <= @l_min")["gene"]
    if l_to_exc.empty:
        print("Threshold is too small to exclude any data.")
        return 
    cols_to_check = ["IGK_v_gene", "IGK_d_gene", "IGK_j_gene", "IGL_v_gene", "IGL_d_gene", "IGL_j_gene"]
    exc_bitmask = df[cols_to_check].isin(l_to_exc.to_list()).any(axis=1)
    new_df = df[~exc_bitmask].copy()
    print("The following light chain genes were removed:", l_to_exc.to_list())
    return new_df
    

def get_gene_counts_df(df, gene_types):
    '''input: 
        df - df where each row has all information for each barcode 
        gene_types - types of genes to consider (MUST have same column name as excel, ex: ['v_gene', 'j_gene', 'd_gene'])
    output: 
        gene_counts_df - pandas serise of total counts of  heavy and light chains in the dataset'''
    needed_columns = [col for col in df.columns if any(col.endswith(gene) for gene in gene_types)]
    all_genes = pd.concat([df[col] for col in needed_columns], ignore_index=True)
    gene_counts_df = all_genes.value_counts().rename_axis('gene').reset_index(name='count')
    return gene_counts_df


def make_combined_rows(df, chain_order): 
    ''' Function reorganize data for each cell, all chains together and orded with heavy chains first. 
    Input: 
        df - dataframe where each chain has its own row 
        chain_order - array with order that heavy and light chains are ordered by (ex: ['IGH', 'IGK', 'IGL'])
    Output: 
        New dataframe where each barcode has all chains in one row with an added index to end of each col name
    '''
    df['chain'] = df['chain'].str.strip()
    df = df.sort_values('barcode')    #sort rows 
    rows = []
    all_barcodes = df['barcode'].unique() 
    repeated_cols = [col for col in df.columns if col not in ['barcode', 'raw_clonotype_id']]
    for b in all_barcodes: 
        bchains_df = df[df['barcode'] == b].copy() #df of only chains releated to barcode b 
        bchains_df = bchains_df.sort_values(by='chain', ascending=True)
        combined_row = {}
        combined_row['barcode'] = b
        combined_row['raw_clonotype_id'] = bchains_df.iloc[0]['raw_clonotype_id']
        for chain in chain_order:
            chain_row = bchains_df[bchains_df['chain']==chain]
            if not chain_row.empty:
                row = chain_row.iloc[0]
                for col in repeated_cols: 
                    if col in row: 
                        combined_row[f"{chain}_{col}"] = row[col]
                    else: 
                        combined_row[f"{chain}_{col}"] = None
            else: 
                for col in repeated_cols: 
                    combined_row[f"{chain}_{col}"] = None
        rows.append(combined_row)
    return pd.DataFrame(rows)

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

def amino_acid_count(data, amino_code, trim_heavy=True):
    '''
    input: 
        - data: df with data
        - amino_code: character of amino acid to count 
        - trim_heavy: optional, usually true but if want to include the first 3 amino characters of heavy sequencies then add False to inputs 
    output: no output from this function but does insert column to data 
    '''
    check_for = "_count" if trim_heavy else "_noTrimCount"
    check_cols = [col for col in data.columns if col.endswith(check_for)]
    if len(check_cols) > 0:
        print("ERROR: This column already exists in this df. If you would like to update these values, please drop the related columns and rerun function. Like so:\n")
        print(f"drop_cols = [col for col in data.columns if col.endswith({check_for})]\n")
        print("data = data.drop(columns=drop_cols)\n")
        return 
    needed_cols = [col for col in data.columns if col.endswith("cdr3")].copy()
    for col in needed_cols: 
        col_name = col[:3] + "_" + amino_code + "_count"
        if (trim_heavy == False):
            col_name = col[:3] + "_" + amino_code + "_noTrimCount"
        new_col = []
        if trim_heavy and col.startswith("IGH"):
            new_col = data[col].apply(
                lambda x: x[3:].count(amino_code) if isinstance(x, str) else pd.NA
            )
        else:
            new_col = data[col].apply(
                lambda x: x.count(amino_code) if isinstance(x, str) else pd.NA
            )
        #inserting new col right next to cdr3 col
        col_index = data.columns.get_loc(col) + 1
        data.insert(col_index, col_name, new_col)

def crd3_pKa(data, trim_heavy=True):
    '''
    input: 
        - data: df with data
        - trim_heavy: optional, usually true but if want to include the first 3 amino characters of heavy sequencies then add False to inputs 
    output: no output from this function but does insert column to data 
    '''
    check_for = "_crd3_pKa" if trim_heavy else "_noTrimCrd3_pKa"
    check_cols = [col for col in data.columns if col.endswith(check_for)]
    if len(check_cols) > 0:
        print("ERROR: This column already exists in this df. If you would like to update these values, please drop the related columns and rerun function. Like so:\n")
        print(f"drop_cols = [col for col in data.columns if col.endswith({check_for})]\n")
        print("data = data.drop(columns=drop_cols)\n")
        return 
    needed_cols = [col for col in data.columns if col.endswith("cdr3")]
    new_cols = {col: [] for col in needed_cols} #if initalize dict like dict.fromkeys(needed_cols, []) then all keys share same list object, instead do this!
    for index, row in data.iterrows():
        for col in needed_cols: 
            if isinstance(row[col], str):
                sequence = ""
                if trim_heavy and col.startswith("IGH"):
                    sequence = row[col][3:] 
                else: 
                    sequence = row[col] 
                counts = Counter(sequence)
                pka = counts['R']*12.48+counts['D']*3.65+counts['C']*8.18+counts['H']*6+counts['K']*10.53+counts['E']*4.25+counts['Y']*10.07
                new_cols[col].append(pka)
            else: 
                new_cols[col].append(None)
    for col in needed_cols: 
        col_name = col[:4] + "_crd3_pKa"
        if (trim_heavy == False):
            col_name = col[:4] + "_noTrimCrd3_pKa"
        col_index = data.columns.get_loc(col) + 1
        data.insert(col_index, col_name, new_cols[col])

def all_heavyLight_genes(df_list, gene_types): 
    '''
    input: 
        - df_list: list of dfs of all files to consider in experiment 
        - gene_types: list of gene types to be considered (MUST BE SAME AS NAME IN FILE)
    output: 
        - all_igh: list of all the heavy chain genes in the experiment 
        - all_iglk: list of all light chain genes in the experiment 
    '''
    all_igh = []
    all_iglk = []
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

def pairDf_to_Matrix(pairs_df, all_heavies, all_lights): 
    pairs_df = pairs_df.set_index(["IGH", "IGKL"])  # so index matches full_pairs
    all_pairs = pd.MultiIndex.from_product([all_heavies, all_lights], names=["IGH", "IGKL"])
    pairs_df = pairs_df.reindex(all_pairs, fill_value='')
    pairs_df = pairs_df.reset_index() #so that can pivot 
    matrix_df = pairs_df.pivot(index="IGH", columns="IGKL", values="Pair Count") #matrix is reshaped DataFrame
    return matrix_df

def get_gene_num(gene):
    ''' Converts string number at end of gene name to int/float'''
    match = re.search(r'-(\d+)', gene)
    if match:
        return int(match.group(1))
    else: #if no dash 
        match = re.search(r'(\d+)$', gene)
        return int(match.group(1)) if match else float('inf')
    

def lightChain_vs_heavyRcount_relFrequency(grouped_df, needed_genes, ighR_colName): 
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
    output_df = output_df.reindex(needed_genes, fill_value=0)
    output_df = output_df.reset_index()
    #relative frequency analysis 
    output_df['total'] = output_df['Rcount_equalTo_0']+output_df['Rcount_gt_0']
    output_df['R_gt_0_freq'] = output_df['Rcount_gt_0']/output_df['total'] *100
    output_df['R_equal_0_freq'] = output_df['Rcount_equalTo_0']/output_df['total'] *100
    return output_df 

#TODO function for z-test USE FOR PAIRS OF LIGHT AND HEAVY CHAINS
#TODO function for bionomial distribution FOR LIGHT CHAINS AND R STUFF
#TODO function for t-test USE FOR COMPAIRING NORMALIZED ABUNDANCE OF LIGHT CHAINS BETWEEN HEALTHY AND SLE

def to_excel(grouped_df, specific_count, frequencies_matrix, file_name):
    with pd.ExcelWriter("tai_"+ file_name + ".xlsx", engine='openpyxl') as writer:
        grouped_df.to_excel(writer, sheet_name="Grouped", index=False,  na_rep="")
        start_col = 0
        #writing gene type frequencies to same excel file but different sheet 
        specific_count.to_excel(writer, sheet_name='Analysis', index=True, startcol=start_col, na_rep="")
        start_col += 4
        #writing pairs in form of frequencies matrix to the same excel sheet 
        frequencies_matrix.to_excel(writer, sheet_name='Analysis', index=True, startcol=start_col, na_rep="")

def all_to_excel(organized_data, file_name): 
    ''' TODO: remove tai_ from the naming before giving to other people to use!
    input: 
        Organized data - dictionary of data to export to file broken up by which sheet to put what where and if a graph of the data is needed or not 
            ex:  { sheet_name: {'tables': [{
                                            "data": df,                                #dataframe/matrix (reshaped DataFrame)
                                            "title": table_title,                      # title to put above table 
                                            "graph_type": None,                        #None if no graph needed, bar, sns.barplot, etc
                                            "graph_name(s)":                           #" " if no graph to name
                                            "graph_details": {"x": "Gene", "y": "another df col name"}
                                            "split_by":[],                             #what values to groupby/split by if breaking df into several tables, [] if keeping df together
                                            "sort_by":[],                              #what columns in the df to sort by if need to sort 
                                            "drop_cols":[],                            #what columns to drop from the df 
                                            "include_cols":[]}                         #what columns to indlude in excel file from df
                                    {add as many tables per sheet as needed}, ...]}, 
                    sheet_name: {'tables': [{type: "one_row",                       #if want several different dfs to be on the same row in sheet do this 
                                             'tables': [{'data':df1, ...}, {"data":df2, ...}]}, 
                                            {"data":df3, ... }]} } 
        file_name - name of excel file to make/overwrite 
    output: 
        there is nothing returned from function but file file_name.xlsx is made in this folder
    '''
    with pd.ExcelWriter("tai_"+ file_name + ".xlsx", engine='openpyxl') as writer: #if tai_filename.xlsx already this will exists its contents
        for sheetName, sheet_content in organized_data.items(): 
            print(f"Writing sheet: {sheetName}")
            start_row, start_col = 0, 0
            for item in sheet_content.get("tables"):
                if item.get("type") == "one_row": #want several df on one row 
                    cur_startCol = start_col
                    for table in item.get("tables"):
                        cur_df = table.get('data')
                        if not table.get("split_by"): #Cases: .get() returns None, split_by is explicitly None, or an empty list ([])
                            shape = writer_table_helper(writer, sheetName, cur_df, start_row, cur_startCol, title=table.get("title"), sort_by=table.get("sort_by"), drop_cols=table.get("drop_cols"), include_cols=table.get("include_cols"))
                            num_col = 2 if len(shape) == 1 else shape[1]
                            if g_type:= table.get("graph_type"):
                                insert_graph_helper(writer, sheetName, cur_df, table["graph_type"], table["graph_details"], start_row + shape[0], table.get("graph_name"))
                            cur_startCol += num_col + 2
                        else: #split up df 
                            for value, group in cur_df.groupby(by=table.get("split_by")):
                                shape = writer_table_helper(writer, sheetName, group, start_row, cur_startCol, title=table.get("split_by")[0]+": "+ value, sort_by=table.get("sort_by"), drop_cols=table.get("drop_cols"), include_cols=table.get("include_cols"))
                                cur_startCol += shape[1]+1
                            #ADD GRAPHING PART HERE TOO!
                    start_col = cur_startCol
                else: 
                    cur_df = item.get("data")
                    if not item.get("split_by"): #keep cur_df together 
                        shape = writer_table_helper(writer, sheetName, cur_df, start_row, start_col, title=item.get("title"), sort_by=item.get("sort_by"), drop_cols=item.get("drop_cols"), include_cols=item.get("include_cols"))
                        if g_type:= item.get("graph_type"):
                            insert_graph_helper(writer, sheetName, cur_df, item["graph_type"], item["graph_details"], start_row + shape[0], item.get("graph_name"))
                        start_row += shape[0] + 4
                    else: #split up df
                        cur_startCol = start_col
                        for value, group in cur_df.groupby(by=item.get("split_by")):
                            shape = writer_table_helper(writer, sheetName, group, start_row, cur_startCol, title=item.get("split_by")[0]+": "+ value, sort_by=item.get("sort_by"), drop_cols=item.get("drop_cols"), include_cols=item.get("include_cols"))
                            cur_startCol += shape[1]+1
                        #TODO:ADD IN ALL THE GRAPH STUFF!
        #TODO: add in colorcoding of rows with significance?? 


'''Helper Functions for alternative_to_excel'''            
def writer_table_helper(writer, sheetName, df, start_row, start_col, title=None, sort_by=None, drop_cols=None, include_cols=None):
        if title: 
            pd.DataFrame([title]).to_excel(writer, sheet_name=sheetName, startcol=start_col, startrow=start_row, header=False, index=False)
            start_row += 1
        if sort_by:
            df = df.sort_values(by=sort_by)
        if drop_cols: 
            df = df.drop(columns=drop_cols)
        use_df = df 
        if include_cols: 
            use_df = df[include_cols]
        use_df.to_excel(writer, sheet_name=sheetName, startcol=start_col, startrow=start_row, index=True, na_rep="")
        return use_df.shape

#WAIT MAYBE https://docs.xlwings.org/en/stable/matplotlib.html CAN MAKE THIS MUCH EASIER!! AND THE GRAPHS WILL BE EDITABLE??
def insert_graph_helper(writer, sheet_name, df, graph_type, graph_details, plot_row, graph_name="", graph_width=25): 
    fig, ax = plt.subplots(figsize=(graph_width, 5))
    x = graph_details.get("x")
    y = graph_details.get("y")
    if (graph_type.lower() == "bar"):
        df.plot(x=x, y=y, ax=ax, kind="bar", title=graph_name)
        ax.set_xlabel("Kind of " + x)
        ax.set_ylabel(y)
    elif (graph_type.lower() == "grouped"):
        hue = graph_details.get("hue")
        sns.barplot(data=df, x=x, y=y, hue=hue, ax=ax)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title=hue+" Types")
        plt.title(graph_name)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    imgdata = BytesIO()
    plt.savefig(imgdata, format='png')
    plt.close(fig)
    imgdata.seek(0)
    img = Image(imgdata)
    writer.sheets[sheet_name].add_image(img, f"A{plot_row}")

'''Not Needed Functions (but dont want to delete them yet)'''
def ungrouped_chain_type_counts(df, gene_types): 
    '''
    input: 
        grouped - df where each row has all information for each barcode 
        gene_types - types of genes to consider (MUST have same column name as excel, ex: ['v_gene', 'j_gene', 'd_gene'])
    output: 
        total_heavyLight_counts - pandas serise of total counts of  heavy and light chains in the dataset
        chain_occurance_counts - pandas serise of total heavy and light chains 
        pairs_df - pandas df of all of the heavy and light pairs in each cells 
    '''
    #counting totals for IGH, IGK, IGL chains
    total_heavyLight_counts = df['chain'].value_counts() 
    #counting frequency of all genes in dataset 
    all_genes = pd.concat([df[gene] for gene in gene_types])
    chain_occurance_counts = all_genes.value_counts()
    #counting pairs of heavy and light chains 
    pair_counts = Counter()
    for barcode, chain_infos in df.groupby("barcode"):
        cell_genes = chain_infos[gene_types].values.flatten()
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
    return total_heavyLight_counts, chain_occurance_counts, pairs_df