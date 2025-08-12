'''SCRATCH WORK THAT IS GOOD REFRENCE SOMETIMES BUT OTHERWISE NOT NEEDED HEHE'''
# file_name = 'sle_3_filtered_contig_annotations.csv'
# df = pd.read_csv(file_name)
# remove_unproductive_rows(df)
# needed_genes = ['v_gene', 'j_gene', 'd_gene']
# chain_order = ['IGH', 'IGK', 'IGL']  
# grouped_df = make_combined_rows(df, chain_order)
# amino_acid_count(grouped_df, 'R')
# crd3_pKa(grouped_df)
# #print(grouped_df)
# total_heavyLight_counts, chain_occurance_counts, pairs_df = ungrouped_chain_type_counts(df, needed_genes)
# all_igh, all_idlk = all_heavyLight_genes(grouped_df, needed_genes)
# light_vs_Rcount_df = lightChain_vs_heavyRcount(grouped_df)
# # print(light_vs_Rcount_df)
# # print(light_vs_Rcount_df.shape)
# matrix_df = pairDf_to_Matrix(pairs_df, all_igh, all_idlk)
# organized_data = {"Grouped": {'tables': [{"data": grouped_df}]}, 
#                   "Analysis": {'tables': [{"type": "one_row", "tables": [
#                                                         {"data": chain_occurance_counts, "title": "Total Counts of Genes in Dataset"}, 
#                                                         {"data": matrix_df, "title": "Frequency Matrix of Light and Heavy Chain Gene Pairs"}] }] }
#                   }
# print(chain_occurance_counts)
# print(len(chain_occurance_counts.shape))
# print(len(chain_occurance_counts.columns))
# print(matrix_df)
# print(matrix_df.shape)
# print(len(matrix_df.columns))
#alternative_to_excel(organized_data, "_alt_"+file_name)
#to_excel(grouped_df, chain_occurance_counts, matrix_df, file_name)


#File from Mamas email TODO: MAKE THIS INTO FUNCTIONS! OR ADJUST CURRENT FUNCTIONS TO BE ABLE TO DO THIS!
# file_name = 'Atypical_heavy_light_chains.csv'
# df = pd.read_csv(file_name)
# gene_types = ['v_gene', 'j_gene', 'd_gene']
# needed_cols = [col for col in df.columns if any(col.endswith(gene) for gene in gene_types)]
# all_heavies, all_lights = all_heavyLight_genes(df, gene_types)
# processed_samples = [] #[[name, [cell type, total heavy/light, specific gene count, pair counts], [...]], [...]]
# for name, sample_info in df.groupby("Sample"):
#     sample = [name]
#     for cell_type, cell_info in sample_info.groupby("Celltype (Cluster)"):
#         #counting totals for IGH, IGK, IGL chains
#         total_heavyLight_counts = cell_info['chain'].value_counts() 
#         #counting frequency of all genes in cell type of sample/person  
#         all_genes = pd.concat([cell_info[gene] for gene in needed_cols], ignore_index=True)
#         gene_groups = all_genes.value_counts().rename_axis('Gene').reset_index(name='Count')
#         gene_groups['Prefix'] = gene_groups['Gene'].str[:4] 
#         gene_groups['sortNum'] = gene_groups['Gene'].apply(get_gene_num)
#         gene_groups['Family Total'] = gene_groups.groupby('Prefix')['Count'].transform("sum")
#         gene_groups['Frequency'] = gene_groups['Count']/gene_groups['Family Total']*100
#         #counting pairs of heavy and light chains 
#         pair_counts = Counter()
#         for index, row in cell_info.iterrows():
#             cell_genes = row[needed_cols].values.flatten()
#             #seperating heavy and light chains
#             heavies = []
#             lights = []
#             for gene in cell_genes: 
#                 if isinstance(gene, str):
#                     if gene.startswith('IGH'):
#                         heavies.append(gene)
#                     elif (gene.startswith('IGK') or gene.startswith('IGL')):
#                         lights.append(gene)
#             for h in heavies: 
#                 for l in lights: 
#                     pair_counts[(h, l)] += 1 
#         #converting counter to dataframe 
#         pairs_df = pd.DataFrame([(k[0], k[1], v) for k, v, in pair_counts.items()], 
#                             columns=["IGH", "IGKL", "Pair Count"]) 
#         sample.append([cell_type, total_heavyLight_counts, gene_groups, pairs_df])
#     #making df of gene usage by cell type, want to be able to group by family as well
#     gene_cells_df = df.melt(id_vars="Celltype (Cluster)", value_vars=needed_cols, var_name='Chain Type', value_name='Gene').dropna()
#     gene_cells_df['Family'] = gene_cells_df['Gene'].str[:4]
#     organized_gene_cells_df = gene_cells_df.groupby(['Family', 'Celltype (Cluster)', 'Gene',]).size().reset_index(name='Count')
#     organized_gene_cells_df['Family Total'] = organized_gene_cells_df.groupby(['Family', 'Celltype (Cluster)'])['Count'].transform('sum')
#     organized_gene_cells_df['Frequency'] = organized_gene_cells_df['Count']/organized_gene_cells_df['Family Total'] *100
#     organized_gene_cells_df['sortNum'] = organized_gene_cells_df['Gene'].apply(get_gene_num)
#     processed_samples.append((sample, organized_gene_cells_df))

# #Sending to excel file
# cur_sample = 0
# with pd.ExcelWriter("tA1_"+ file_name + ".xlsx", engine='openpyxl') as writer:
#     for sample, gene_cells in processed_samples: 
#         sample_name = sample[0]
#         start_row = 1
#         for cells in sample[1:]: 
#             start_col = 0
#             #writing cell type label 
#             cell_name_df = pd.DataFrame([f"Cell: Type: {cells[0]}"])
#             cell_name_df.to_excel(writer, sheet_name=sample_name, startrow=start_row, startcol=start_col, header=False, index=False)
#             start_row += 1
#             #writing chain type frequencies to same excel file with 2 blank columns inbetween 
#             cells[1].to_excel(writer, sheet_name=sample_name, index=True, startrow=start_row, startcol=start_col)
#             start_col += 4
#             #writing chain gene frequencies to same excel file 
#             for prefix, gene_group in cells[2].groupby('Prefix'):
#                 sorted_group = gene_group.sort_values(by=['Prefix', 'sortNum']) 
#                 sorted_group = sorted_group.drop(columns=['sortNum'])
#                 sorted_group.to_excel(writer, sheet_name=sample_name, columns=['Gene', 'Count', 'Family Total', 'Frequency'], index=False, startrow=start_row, startcol=start_col)
#                 start_col += 7 #space for 3 rows plus a blank
#                 #ploting the table
#                 plot_row = start_row + sorted_group.shape[0]
#                 fig, ax = plt.subplots()
#                 plot = sorted_group.plot(x='Gene', y='Frequency', ax=ax, kind='bar', xlabel="Gene Type", ylabel="Frequency",  title=f"{prefix} Frequency")
#                 ax.set_ylim(0, 100)
#                 plt.tight_layout()
#                 imgdata = BytesIO() #dont need to save images to disk this way (i think)
#                 plt.savefig(imgdata, format='png')
#                 plt.close(fig)
#                 imgdata.seek(0)
#                 sheet = writer.book[sample_name]
#                 img = Image(imgdata)
#                 img.anchor = f'E{plot_row + 1}'
#                 sheet.add_image(img)
#             start_col += 4
#             #writing pairs in form of frequencies matrix to the same excel sheet
#             frequencies_matrix = pairDf_to_Matrix(cells[3], all_heavies, all_lights)
#             frequencies_matrix.to_excel(writer, sheet_name=sample_name, index=True, startrow=start_row, startcol=start_col, na_rep="")
#             start_row += frequencies_matrix.shape[0] + 3

#         #adding gene type usage by cell type to new sheets
#         #gene_cells.to_excel(writer, sheet_name=sample_name, startrow=start_row, index=False)
#         start_col = cur_sample*7
#         for fam_name, fam_data in gene_cells.groupby('Family'):
#             sorted_fam = fam_data.sort_values(by=['Family', 'sortNum'])
#             sorted_fam = sorted_fam.drop(columns=['sortNum'])
#             pd.DataFrame([[f"Sample: {sample_name}"]]).to_excel(
#                 writer, sheet_name=fam_name, startrow=0, startcol=start_col, index=False, header=False
#             )
#             sorted_fam.to_excel(writer, sheet_name=fam_name, index=False, startrow=1, startcol=start_col)
#             #ploting charts
#             plot_row = start_row + sorted_fam.shape[0]
#             fig, ax = plt.subplots(figsize=(25, 5)) #create new blank figure with enough space for gene names
#             sns.barplot(data=sorted_fam, x='Gene', y='Frequency', hue='Celltype (Cluster)', ax=ax)
#             ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="Cell Type")
#             #ax.set_ylim(0, 68)
#             plt.title(f"{fam_name} Gene Usage By Cell Type in {sample_name}")
#             plt.xticks(rotation=45, ha='right') #rotate x-axis labels and right aligns them so easier to read 
#             plt.tight_layout() #adjest spacing to avoid overlap 
#             imgdata = BytesIO() #creats in-memory bytes stream/virual file so dont need to save to disk
#             plt.savefig(imgdata, format='png')
#             plt.close(fig)
#             imgdata.seek(0) #move pointer to start of bytesIO stream 
#             img = Image(imgdata)
#             img.anchor=f'B{3}'
#             sheet = writer.book[fam_name]
#             sheet.add_image(img)
#         cur_sample += 1 

# '''Not Needed Functions (but dont want to delete them yet)'''
# def ungrouped_chain_type_counts(df, gene_types): 
#     '''
#     input: 
#         grouped - df where each row has all information for each barcode 
#         gene_types - types of genes to consider (MUST have same column name as excel, ex: ['v_gene', 'j_gene', 'd_gene'])
#     output: 
#         total_heavyLight_counts - pandas serise of total counts of  heavy and light chains in the dataset
#         chain_occurance_counts - pandas serise of total heavy and light chains 
#         pairs_df - pandas df of all of the heavy and light pairs in each cells 
#     '''
#     #counting totals for IGH, IGK, IGL chains
#     total_heavyLight_counts = df['chain'].value_counts() 
#     #counting frequency of all genes in dataset 
#     all_genes = pd.concat([df[gene] for gene in gene_types])
#     chain_occurance_counts = all_genes.value_counts()
#     #counting pairs of heavy and light chains 
#     pair_counts = Counter()
#     for barcode, chain_infos in df.groupby("barcode"):
#         cell_genes = chain_infos[gene_types].values.flatten()
#         #seperating heavy and light chains
#         heavies = []
#         lights = []
#         for gene in cell_genes: 
#             if isinstance(gene, str):
#                 if gene.startswith('IGH'):
#                     heavies.append(gene)
#                 elif (gene.startswith('IGK') or gene.startswith('IGL')):
#                     lights.append(gene)
#         for h in heavies: 
#             for l in lights: 
#                 pair_counts[(h, l)] += 1
#     #converting counter to dataframe 
#     pairs_df = pd.DataFrame([(k[0], k[1], v) for k, v, in pair_counts.items()], 
#                             columns=["IGH", "IGKL", "Pair Count"])       
#     return total_heavyLight_counts, chain_occurance_counts, pairs_df
