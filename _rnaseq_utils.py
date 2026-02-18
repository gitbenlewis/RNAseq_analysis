'''
rnaseq_utils
'''
# module at /home/ubuntu/projects/gitbenlewis/RNAseq_analysis/rnaseq_utils.py
# module level imports
import os
import os.path
import pandas as pd
import numpy as np
import re
import shutil
import fnmatch
import re


##### make_samplesheet ##### start
def make_samplesheet(FASTQ_dir, output_file='samplesheet.csv',index_col_name='sample', save_output=True,short_test_list=False):
    '''
    This function takes a Azenta fastaq output directory with fastq files and makes a samplesheet with the following columns:
    sample, fastq_1, fastq_2, strandedness
    This is the format required by nf-core RNAseq pipeline
    Cauation may need to be modified for fastq directory with different naming conventions or contents
    '''
    list_of_fastqs=os.listdir(FASTQ_dir)
    rm_pat='*md5*'
    # Using list comprehension and fnmatch.filter to remove matching strings
    list_of_fastqs = [item for item in list_of_fastqs if not fnmatch.fnmatch(item, rm_pat)]
    pattern = r'_R.*?_001\.fastq\.gz'
    list_of_samples=[re.split(pattern,item)[0] for item in list_of_fastqs]
    # now remove duplicates in sample list
    # list_of_samples=list(set(list_of_samples))
    # now remove duplicates in sample list without changing order
    list_of_samples = list(dict.fromkeys(list_of_samples))
    # now make fastq_1 and fastq_2 lists
    fastq_1_list=[FASTQ_dir+'/'+item+'_R1_001.fastq.gz' for item in list_of_samples]
    fastq_2_list=[FASTQ_dir+'/'+item+'_R2_001.fastq.gz' for item in list_of_samples]

    # now make a samplesheet dataframe
    samplesheet_df=pd.DataFrame({index_col_name:list_of_samples,'fastq_1':fastq_1_list,'fastq_2':fastq_2_list})
    samplesheet_df['strandedness']='auto'
    # now save the samplesheet to csv file with index false and header true
    if save_output:
        samplesheet_df.to_csv(output_file,index=False,header=True)
    if short_test_list:
        # now make a short test list with only 2 rows
        samplesheet_df=samplesheet_df.head(2)
        samplesheet_df.to_csv(output_file,index=False,header=True)
    return samplesheet_df

def make_samplesheet_nonAzenta(FASTQ_dir, output_file='samplesheet.csv',index_col_name='sample', save_output=True,short_test_list=False):
    '''
    This function takes a Non regular output directory (not from Azenta)
      with fastq files and makes a samplesheet with the following columns:
    sample, fastq_1, fastq_2, strandedness
    This is the format required by nf-core RNAseq pipeline
    '''
    ###
    list_of_fastqs=os.listdir(FASTQ_dir)
    rm_pat='*md5*'
    # Using list comprehension and fnmatch.filter to remove matching strings
    list_of_fastqs = [item for item in list_of_fastqs if not fnmatch.fnmatch(item, rm_pat)]
    rm_pat='*.fastq'
    list_of_fastqs = [item for item in list_of_fastqs if not fnmatch.fnmatch(item, rm_pat)]
    # make list of libraries
    pattern = r'_R.*?_001\.fastq\.gz'
    list_of_libraries=[re.split(pattern,item)[0] for item in list_of_fastqs]
    # now remove duplicates in list_of_libraries without changing order
    list_of_libraries = list(dict.fromkeys(list_of_libraries))
    # now make fastq_1 and fastq_2 lists
    fastq_1_list=[FASTQ_dir+'/'+item+'_R1_001.fastq.gz' for item in list_of_libraries]
    fastq_2_list=[FASTQ_dir+'/'+item+'_R2_001.fastq.gz' for item in list_of_libraries]
    # make list of samples
    pattern = r'_L00.*?_R.*?_001\.fastq\.gz'
    list_of_samples=[re.split(pattern,item)[0] for item in list_of_libraries]
    # now make a samplesheet dataframe
    samplesheet_df=pd.DataFrame({index_col_name:list_of_samples,'fastq_1':fastq_1_list,'fastq_2':fastq_2_list})
    samplesheet_df['strandedness']='auto'
    # now order the sample list by the sample column
    samplesheet_df=samplesheet_df.sort_values(by=[index_col_name])
    # now save the samplesheet to csv file with index false and header true
    if save_output:
        samplesheet_df.to_csv(output_file,index=False,header=True)
    if short_test_list:
        # now make a short test list with only 2 rows
        samplesheet_df=samplesheet_df.head(2)
        samplesheet_df.to_csv(output_file,index=False,header=True)
    return samplesheet_df



# add Random_AB column with half of the rows as A and half as B
def add_random_AB_column(df, mask=None, label_A='A', label_B='B', random_seed=None):
    """Add randomized labels for masked rows; preserves original row order."""
    df = df.copy()
    if mask is None:
        mask = pd.Series(True, index=df.index)
    elif isinstance(mask, pd.Series):
        mask = mask.reindex(df.index, fill_value=False)
    else:
        if len(mask) != len(df):
            raise ValueError("mask must be same length as df")
        mask = pd.Series(mask, index=df.index)
    mask = mask.fillna(False)

    col = f'Random_{label_A}_{label_B}'
    df[col] = pd.NA

    row_positions = np.flatnonzero(mask.to_numpy(dtype=bool))
    n_rows = row_positions.size
    if not n_rows:
        return df

    n_A = n_rows // 2
    rng = np.random.default_rng(seed=random_seed)
    rng.shuffle(row_positions)

    col_idx = df.columns.get_loc(col)
    df.iloc[row_positions[:n_A], col_idx] = label_A
    df.iloc[row_positions[n_A:], col_idx] = label_B
    return df


##### make_samplesheet ##### END


#### Start ################################ Editing the NFcore outputs. ####################################


################################ This function adds the mouse homolog to the human gene names gene_id2gene_name
# function to make gene_id2gene_name file
def make_gene_id2gene_name_file(nfcore_output_rnaseq_dir='rnaseq/results/',
                                add_homologs=True,organism='human',
                                h2m_agg_file='/data/h2m_agg.csv',
                                m2h_agg_file='/data/m2h_agg.csv',
                                save_output=True,
                                output_dir=None,
                                also_save_to_src_dir=True,
                                return_df=False,
                                ):
    '''This function reads the tx2gene_file and makes a gene_id2gene_name file
    add_homologs=True,organism='human','mouse'
    parameters:
        nfcore_output_rnaseq_dir (str): Path to the nfcore rnaseq output directory.
        add_homologs (bool): Whether to add homologs to the gene_id2gene_name file.
        organism (str): 'human' or 'mouse' to indicate the reference organism for homolog annotations.
        h2m_agg_file (str): Path to the human to mouse homolog aggregated file.
        m2h_agg_file (str): Path to the mouse to human homolog aggregated file.
        save_output (bool): Whether to save the gene_id2gene_name file.
        output_dir (str): Directory to save the gene_id2gene_name file. If None, saves to nfcore_output_rnaseq_dir/star_salmon.
        return_df (bool): Whether to return the gene_id2gene_name dataframe.
    returns:
        gene_id_df_species2species (pd.DataFrame): Dataframe with gene_id to gene_name and homologs if return_df is True.

    '''
    count_dir=os.path.join(nfcore_output_rnaseq_dir,'star_salmon')
    output_path = os.path.join(output_dir if output_dir is not None else count_dir, 'gene_id2gene_name.csv')
    count_dir_output_path = os.path.join(count_dir, 'gene_id2gene_name.csv')
    tx2gene_file=os.path.join(count_dir,'tx2gene.tsv')
    tx2gene_df=pd.read_csv(tx2gene_file,sep='\t',header=None)
    tx2gene_df.columns=['transcript_id','gene_id','gene_name']
    ### now keep the rows with unique gene_id (keep first) and drop the transcript_id column
    tx2gene_df=tx2gene_df.drop_duplicates(subset='gene_id',keep='first')
    tx2gene_df=tx2gene_df.drop(columns=['transcript_id'])
    # re set index with droping the existing index
    tx2gene_df.reset_index(inplace=True,drop=True)
    #print(tx2gene_df.shape)
    # now make column with unique gene names the second instance of a non unique gene name will have the gene_id appended to the gene_name
    isduplicated=tx2gene_df.duplicated(subset="gene_name",keep='first')
    unique_gene_names=[tx2gene_df['gene_name'][i]+'_'+tx2gene_df['gene_id'][i] if isduplicated[i] else tx2gene_df['gene_name'][i] for i in range(len(tx2gene_df))]
    tx2gene_df['unique_gene_name']=unique_gene_names
    if add_homologs==False:
        gene_id_df_species2species = tx2gene_df
        if save_output:
            tx2gene_df.to_csv(output_path, index=False)
            if also_save_to_src_dir and output_dir is not None:
                tx2gene_df.to_csv(count_dir_output_path, index=False)
    else:
        if organism=='human':
            h2m_aggregated_df=pd.read_csv(h2m_agg_file)
            cols2merge=['gene_id','mouse_gene_id', 'mouse_gene_name', 'mouse_gene_ids_all', 'mouse_gene_names_all']
            species2species_map_df=h2m_aggregated_df.copy()
        elif organism=='mouse':
            m2h_aggregated_df=pd.read_csv(m2h_agg_file)
            cols2merge=['gene_id','human_gene_id', 'human_gene_name', 'human_gene_ids_all', 'human_gene_names_all']
            species2species_map_df=m2h_aggregated_df.copy()
        # only keep cols in cols2merge that are in species2species_map_df
        cols2merge=[col for col in cols2merge if col in species2species_map_df.columns]
        # merge the tx2gene_df with the species2species_map_df on gene_id
        gene_id_df_species2species=pd.merge(tx2gene_df,species2species_map_df[cols2merge],left_on='gene_id',right_on='gene_id',how='left')
        if save_output:
            gene_id_df_species2species.to_csv(output_path, index=False)
            if also_save_to_src_dir and output_dir is not None:
                gene_id_df_species2species.to_csv(count_dir_output_path, index=False)
    if return_df:
        return gene_id_df_species2species
    else:
        return None


################################ functions to filter  the basic NF core RNA seq pipeline out puts prior to downstream analysis

def filter_zeros_n_rowaverage_of_nfcore_output_dir(
        nfcore_output_rnaseq_dir='rnaseq/results/',
        min_rowaverage=10,
        min_non_zeros_perrow=1,
        save_output=True,
        output_dir=None,
        also_save_to_src_dir=True,
        return_dfs=False,): 
    '''
    Filters the nfcore RNAseq output count files based on  minimum row average and minimum number of non-zeros per row .
    opperates on 'salmon.merged.gene_counts.tsv','salmon.merged.gene_counts_length_scaled.tsv','salmon.merged.gene_lengths.tsv'
    parameters:
        nfcore_output_rnaseq_dir (str): Path to the nfcore rnaseq output directory.
        min_rowaverage (int): Minimum row average to filter rows.
        min_non_zeros_perrow (int): Minimum number of non-zeros per row to filter rows.
        save_output (bool): Whether to save the filtered count files.
        output_dir (str): Directory to save the filtered count files. If None, saves to nfcore_output_rnaseq_dir/star_salmon.
        return_df (bool): Whether to return the filtered dataframes.
    returns:
        raw_count_file_df_filtered,length_scaled_count_file_df_filtered, gene_lengths_file_df_filtered: Filtered dataframes.


    This function reads the count files from the nfcore output directory and filters the rows based on the minimum number of non-zeros per row and the row average.
    The function reads the raw count file and the length scaled count file and filters the rows in both files.
    The function then reads the gene lengths file and removes rows that are not in the filtered raw count file.
    The function returns the filtered length scaled count file, the filtered raw count file and the filtered gene lengths file.
    '''
    import os
    ###### make a dictionary of the parameters used in the function and generate above
    parameters = locals().copy()

    def filter_zeros_n_rowaverage(df,  min_rowaverage=10,  min_non_zeros_perrow=1,**kwargs): 
        if min_rowaverage is not None:
            df_filtered = df[df.iloc[:, 2:].mean(axis=1) > min_rowaverage].copy()
            print(f'{df.shape[0] - df_filtered.shape[0]} Rows with row average less than {min_rowaverage} filtered,  {df_filtered.shape[0]} rows remain')
            if min_non_zeros_perrow is not None:
                rows_pre_filter = df_filtered.shape[0]
                df_filtered = df_filtered[(df_filtered.iloc[:, 2:] != 0).sum(axis=1) >= min_non_zeros_perrow].copy()
                print(f' {rows_pre_filter - df_filtered.shape[0]} Rows with fewer than {min_non_zeros_perrow} non-zeros filtered,  {df_filtered.shape[0]} rows remain')
                print(df_filtered.shape, 'is the final shape of the df_filtered')
                return df_filtered
            else:
                print('No filtering based on non-zeros (min_non_zeros_perrow=None) final returned df shape is', df_filtered.shape)
                return df_filtered
        else:
            print('No filtering done  (min_rowaverage=None and min_non_zeros_perrow=None) final returned shape is', df.shape)
            df_filtered = df
            return df_filtered

    # paths to count files
    count_dir = os.path.join(nfcore_output_rnaseq_dir, 'star_salmon')
    raw_count_file = os.path.join(count_dir, 'salmon.merged.gene_counts.tsv')
    length_scaled_count_file = os.path.join(count_dir, 'salmon.merged.gene_counts_length_scaled.tsv')
    gene_lengths_file = os.path.join(count_dir, 'salmon.merged.gene_lengths.tsv')
    # read count files
    raw_count_file_df = pd.read_csv(raw_count_file, sep='\t')
    print('raw_count_file_df shape is', raw_count_file_df.shape)
    length_scaled_count_file_df = pd.read_csv(length_scaled_count_file, sep='\t')
    print('length_scaled_count_file_df shape is', length_scaled_count_file_df.shape)
    # filter count files
    print('Filtering raw count file')
    raw_count_file_df_filtered = filter_zeros_n_rowaverage(raw_count_file_df, 
                                                           min_rowaverage=min_rowaverage, min_non_zeros_perrow=min_non_zeros_perrow,# **parameters
                                                           )
    print('Filtering length scaled count file')
    length_scaled_count_file_df_filtered = filter_zeros_n_rowaverage(length_scaled_count_file_df,
                                                                      min_rowaverage=min_rowaverage, min_non_zeros_perrow=min_non_zeros_perrow,# **parameters
                                                                     )

    # remove rows from the gene_lengths_file that are not in the filtered filtered raw count file
    gene_lengths_file_df = pd.read_csv(gene_lengths_file, sep='\t')
    print(gene_lengths_file_df.shape, ' is gene_lengths_file_df shape')
    # remove rows from the gene_lengths_file by selecting the index of the filtered raw count file (use the gene_id column they are unique)
    gene_lengths_file_df_filtered = gene_lengths_file_df[gene_lengths_file_df['gene_id'].isin(raw_count_file_df_filtered['gene_id'])]
    print('gene_lengths_file_df_filtered shape is (genes filtered from raw_count_file removed)', gene_lengths_file_df_filtered.shape)

    if save_output:
        final_output_dir = os.path.join(output_dir if output_dir is not None else count_dir,)
        raw_count_file_df_filtered.to_csv(os.path.join(final_output_dir, f'salmon.merged.gene_counts.filtered_minrowavg{min_rowaverage}_minnz{min_non_zeros_perrow}.tsv'), index=False, sep='\t')
        length_scaled_count_file_df_filtered.to_csv(os.path.join(final_output_dir, f'salmon.merged.gene_counts_length_scaled.filtered_minrowavg{min_rowaverage}_minnz{min_non_zeros_perrow}.tsv'), index=False, sep='\t')
        gene_lengths_file_df_filtered.to_csv(os.path.join(final_output_dir, f'salmon.merged.gene_lengths.filtered_minrowavg{min_rowaverage}_minnz{min_non_zeros_perrow}.tsv'), index=False, sep='\t')
        print('filtered files saved to ', final_output_dir)
        print('filtered files saved with appended to file name .filtered_bm{min_rowaverage}.tsv')
        print(f'salmon.merged.gene_counts.filtered_minrowavg{min_rowaverage}_minnz{min_non_zeros_perrow}.tsv')
        print(f'salmon.merged.gene_counts_length_scaled.filtered_minrowavg{min_rowaverage}_minnz{min_non_zeros_perrow}.tsv')
        print(f'salmon.merged.gene_lengths.filtered_minrowavg{min_rowaverage}_minnz{min_non_zeros_perrow}.tsv')
        print('##################### end ##########################')
        if also_save_to_src_dir and output_dir is not None:
            raw_count_file_df_filtered.to_csv(os.path.join(count_dir, f'salmon.merged.gene_counts.filtered_minrowavg{min_rowaverage}_minnz{min_non_zeros_perrow}.tsv'), index=False, sep='\t')
            length_scaled_count_file_df_filtered.to_csv(os.path.join(count_dir, f'salmon.merged.gene_counts_length_scaled.filtered_minrowavg{min_rowaverage}_minnz{min_non_zeros_perrow}.tsv'), index=False, sep='\t')
            gene_lengths_file_df_filtered.to_csv(os.path.join(count_dir, f'salmon.merged.gene_lengths.filtered_minrowavg{min_rowaverage}_minnz{min_non_zeros_perrow}.tsv'), index=False, sep='\t')

    else:
        print('filtered files not saved, set save_output=True to save files')
        print('##################### end ##########################')
    if return_dfs:
        return raw_count_file_df_filtered,length_scaled_count_file_df_filtered, gene_lengths_file_df_filtered
    else:
        return None


def make_adata_nfcore_rnaseq(
        nfcore_output_rnaseq_dir='rnaseq/results/',
        metadata_file_path='metadata.csv',
        batch_key='batch',
        batch_label='bulk_RNAseq_batch',
        save_h5ad=False,
        output_h5ad_file_name='nfcore_bulk_rnaseq_adata.h5ad',
        output_dir=None,
        also_save_to_src_dir=True,
        ):
    '''
    This function makes an AnnData object from the nfcore rnaseq output directory
    parameters:
        nfcore_output_rnaseq_dir (str): Path to the nfcore rnaseq output directory.
        metadata_file_path (str): Path to the metadata file.
        batch_key (str): Key in the adata.obs to use for batch information.
        batch_label (str): Label to add to the adata.obs under key default='batch'
        save_h5ad (bool): Whether to save the AnnData object to an h5ad file.
        output_h5ad_file_name (str): Name of the h5ad file.
        output_dir (str): Directory to save the h5ad file. if None, saves to nfcore_output_rnaseq_dir/output_h5ad_file_name
    returns:
        bulk_adata_new (anndata.AnnData): AnnData object with the counts and metadata.
    '''
    import anndata
    ########### generate file paths
    count_dir=os.path.join(nfcore_output_rnaseq_dir,'star_salmon')
    gene_counts_file=os.path.join(count_dir,'salmon.merged.gene_counts.tsv')
    gene_counts_scaled_file=os.path.join(count_dir,'salmon.merged.gene_counts_scaled.tsv')
    gene_counts_length_scaled_file=os.path.join(count_dir,'salmon.merged.gene_counts_length_scaled.tsv')
    gene_tpm_file=os.path.join(count_dir,'salmon.merged.gene_tpm.tsv')
    gene_lengths_file=os.path.join(count_dir,'salmon.merged.gene_lengths.tsv')
    tx2gene_file=os.path.join(count_dir,'tx2gene.tsv')

    ######################################### Make counts array  start
    # read the gene counts file
    gene_counts_df=pd.read_csv(gene_counts_file,sep='\t',index_col=0)
    X = gene_counts_df.iloc[:,1:].values.T.copy()
    ######################################### Make counts array  END

    ######################################### Make gene names and var df  start

    # keep only the index and first colunm of the gene_counts_df
    genes = gene_counts_df[[gene_counts_df.columns[0]]].copy()
    # re set index without droping the existing index
    genes.reset_index(inplace=True)
    # now maker column with unique gene names the second instance of a non unique gene name will have the gene_id appended to the gene_name
    isduplicated=genes.duplicated(subset="gene_name",keep='first')
    unique_gene_names=[genes['gene_name'][i]+'_'+genes['gene_id'][i] if isduplicated[i] else genes['gene_name'][i] for i in range(len(genes))]
    genes['unique_gene_name']=unique_gene_names
    # set the unique_gene_names as the index
    genes.set_index('unique_gene_name',inplace=True)

    ######################################### Make gene names and var df  END

    ######################################### Make obs df start
    # read metadata file and set index to the sample column
    obs=pd.read_csv(metadata_file_path,index_col='sample')
    # drop duplicate rows with the same index
    obs=obs[~obs.index.duplicated(keep='first')]
    obs=obs.reindex(gene_counts_df.columns.tolist()[1:])
    obs.head(2)
    ################################ END ######### Make obs df  END 

    ######################################### Make adata object  
    bulk_adata_new = anndata.AnnData(X = X, var =genes, obs = obs)

    ##########################################  save featureCounts raw counts to adata layer featureCounts_counts
    bulk_adata_new.layers["salmon_counts"] = bulk_adata_new.X.copy()
    ######################################### Make adata object END


    ######################################### Make  gene_lengths array  start
    # read the gene lengths file
    gene_lengths_df=pd.read_csv(gene_lengths_file,sep='\t',index_col=0)
    gene_lengths_df = gene_lengths_df.fillna(0)## repalce nan s with ZEROs
    X0 = gene_lengths_df.iloc[:,1:].values.T.copy()
    ######  set active data aray to Batch_counts_salmon_merge_df array 
    bulk_adata_new.X=X0.copy()
    ##########  save salmon raw counts to adata layer salmon_counts
    bulk_adata_new.layers["salmon_gene_lengths"] = bulk_adata_new.X.copy()
    ######################################### Make gene_lengths array  END



    ######################################### Make  gene_counts_scaled array  start
    # read the gene counts scaled file
    gene_counts_scaled_df=pd.read_csv(gene_counts_scaled_file,sep='\t',index_col=0)
    gene_counts_scaled_df = gene_counts_scaled_df.fillna(0)## repalce nan s with ZEROs
    X1 = gene_counts_scaled_df.iloc[:,1:].values.T.copy()
    ######  set active data aray to Batch_counts_salmon_merge_df array 
    bulk_adata_new.X=X1.copy()
    ##########  save salmon raw counts to adata layer salmon_counts
    bulk_adata_new.layers["salmon_gene_counts_scaled"] = bulk_adata_new.X.copy()
    ######################################### Make gene_counts_scaled array  END



    ######################################### Make  gene_counts_length_scaled_file array  start
    gene_counts_length_scaled_df=pd.read_csv(gene_counts_length_scaled_file,sep='\t',index_col=0)

    gene_counts_length_scaled_df = gene_counts_length_scaled_df.fillna(0)## repalce nan s with ZEROs
    X2 = gene_counts_length_scaled_df.iloc[:,1:].values.T.copy()
    ######  set active data aray to Batch_counts_salmon_merge_df array 
    bulk_adata_new.X=X2.copy()
    ##########  save salmon raw counts to adata layer salmon_counts
    bulk_adata_new.layers["salmon_gene_counts_length_scaled"] = bulk_adata_new.X.copy()
    ######################################### Make gene_counts_length_scaled_file array  END


    ######################################### Make salmon_effective_TPM array  start
    # read the gene tpm file
    gene_tpm_df=pd.read_csv(gene_tpm_file,sep='\t',index_col=0)
    gene_tpm_df = gene_tpm_df.fillna(0)## repalce nan s with ZEROs
    X3 = gene_tpm_df.iloc[:,1:].values.T.copy()
    ######  set active data aray to Batch_counts_salmon_merge_df array 
    bulk_adata_new.X=X3.copy()
    ##########  save salmon raw counts to adata layer salmon_counts
    bulk_adata_new.layers["salmon_effective_TPM"] = bulk_adata_new.X.copy()
    ######################################### Make salmon_effective_TPM array  END
    

    ############ add batch label
    bulk_adata_new.obs[batch_key]=batch_label
    print(bulk_adata_new)

    if save_h5ad:
        final_output_dir = os.path.join(output_dir if output_dir is not None else nfcore_output_rnaseq_dir,)
        output_h5ad_file_path=os.path.join(final_output_dir,output_h5ad_file_name)
        bulk_adata_new.write_h5ad(output_h5ad_file_path)
        print(f'bulk_adata_new saved to {output_h5ad_file_path}')
        if also_save_to_src_dir and output_dir is not None:
            output_h5ad_file_path_src=os.path.join(nfcore_output_rnaseq_dir,output_h5ad_file_name)
            bulk_adata_new.write_h5ad(output_h5ad_file_path_src)
            print(f'bulk_adata_new also saved to source dir {output_h5ad_file_path_src}')

    return bulk_adata_new



########## START ################   functions to edit NF DEseq2 output files


def add_gene_names_all_to_deseq2_tables_nfcore_differentialabundance(
        nfcore_DA_output_dir='differentialabundance/results/',
        nfcore_output_rnaseq_dir='rnaseq/results/',
        input_file_suffix='.deseq2.results.tsv',
        input_file_separator='\t',
        output_file_suffix='.deseq2.results.csv',
        run_make_gene_id2gene_name_file=False,
        add_homologs=True,
        organism='human',
        h2m_map_file='/home/ubuntu/data/ref/gene_lists/h2m_agg.csv',
        m2h_map_file='/home/ubuntu/data/ref/gene_lists/m2h_agg.csv',
        save_output=True,
        output_dir=None,
        also_save_to_src_dir=True,
        ):
    '''This function adds gene names to all the deseq2 tables in the deseq2_tables_dir_differential
    parameters:
        nfcore_DA_output_dir (str): Path to the nfcore differential abundance output directory.
        nfcore_output_rnaseq_dir (str): Path to the nfcore rnaseq output directory.
        input_file_suffix (str): Suffix of the input deseq2 files to add gene names to.
        input_file_separator (str): Separator of the input deseq2 files.
        output_file_suffix (str): Suffix of the output deseq2 files with gene names added.
        run_make_gene_id2gene_name_file (bool): Whether to run the make_gene_id2gene_name_file function if gene_id2gene_name file is not found.
        add_homologs (bool): Whether to add homologs to the gene_id2gene_name file.
        organism (str): 'human' or 'mouse' to indicate the reference organism for homolog annotations.
        h2m_map_file (str): Path to the human to mouse homolog aggregated file.
        m2h_map_file (str): Path to the mouse to human homolog aggregated file.
        save_output (bool): Whether to save the output deseq2 files with gene names added.
        output_dir (str): Directory to save the output deseq2 files. If None, saves to nfcore_DA_output_dir/tables/differential.
    returns:
        None
    '''
    import os
    import pandas as pd

    # read the nfcore rnaseq output files
    count_dir=os.path.join(nfcore_output_rnaseq_dir,'star_salmon')
    #tx2gene_file=os.path.join(count_dir,'tx2gene.tsv')
    deseq_output_dir_tables = os.path.join(nfcore_DA_output_dir, 'tables', )
    deseq2_tables_dir_differential=os.path.join(deseq_output_dir_tables,'differential')
    #deseq2_tables_dir_processed_abundance=os.path.join(deseq_output_dir_tables,'processed_abundance')
    # Prefer output_dir gene map if provided; fall back to star_salmon directory
    if output_dir is not None:
        gene_id2gene_name_file_name = os.path.join(output_dir, 'gene_id2gene_name.csv')
        if not os.path.exists(gene_id2gene_name_file_name):
            fallback = os.path.join(nfcore_output_rnaseq_dir, 'star_salmon', 'gene_id2gene_name.csv')
            if os.path.exists(fallback):
                gene_id2gene_name_file_name = fallback
            else:
                print(f'gene_id2gene_name file not found at {gene_id2gene_name_file_name}')
                print(f'and no fallback at {fallback}; exiting add_gene_names_to_deseq2_tables function')
                return
        gene_id2gene_name_df = pd.read_csv(gene_id2gene_name_file_name)
    else:
        try:
            gene_id2gene_name_file_name= os.path.join(nfcore_output_rnaseq_dir,'star_salmon','gene_id2gene_name.csv')
            gene_id2gene_name_df=pd.read_csv(gene_id2gene_name_file_name)
        except:
            print(f'gene_id2gene_name file not found at {gene_id2gene_name_file_name}')
            if run_make_gene_id2gene_name_file:
                print('running make_gene_id2gene_name_file to generate gene_id2gene_name file')
                make_gene_id2gene_name_file(
                    nfcore_output_rnaseq_dir,
                    add_homologs=add_homologs,
                    organism=organism,
                    h2m_agg_file=h2m_map_file,
                    m2h_agg_file=m2h_map_file,
                )
                ### now read the new gene_id2gene_name.csv file
                gene_id2gene_name_df=pd.read_csv(gene_id2gene_name_file_name)
            else:
                print('gene_id2gene_name file not found and run_make_gene_id2gene_name_file is set to False')
                print('exiting add_gene_names_to_deseq2_tables function')
                return

    # ###############   loop to add gene names to all the files and save them as .csvs
    try:
        # list of all the files in the deseq2_tables_dir_differential
        deseq2_tables_dir_differential_files=os.listdir(deseq2_tables_dir_differential)
        for file in deseq2_tables_dir_differential_files:
            if file.endswith(input_file_suffix):
                # make a new file name that ends csv
                new_file=file.replace(input_file_suffix,output_file_suffix)
                table_df=pd.read_csv(os.path.join(deseq2_tables_dir_differential,file),sep=input_file_separator)
                print( 'the original deseq2 file is named',file, 'and has shape',table_df.shape)
                table_df=pd.merge(table_df,gene_id2gene_name_df,left_on='gene_id',right_on='gene_id',how='left')
                # sort the table by padj
                table_df=table_df.sort_values(by='padj')
                if save_output:
                    final_output_dir = os.path.join(output_dir if output_dir is not None else deseq2_tables_dir_differential,)
                    table_df.to_csv(os.path.join(final_output_dir,new_file),index=False)
                    print( 'the new  deseq2 file is named',new_file, 'and has shape',table_df.shape)
                    if also_save_to_src_dir and output_dir is not None:
                        table_df.to_csv(os.path.join(deseq2_tables_dir_differential,new_file),index=False)
                        print( 'the new  deseq2 file is also saved to source dir named',new_file, 'and has shape',table_df.shape)
    except:
        print('error in add_gene_names_to_deseq2_tables')


def add_Rank_Metric_S(
        deseq2_file_path,
        pvalue_col_label='pvalue',l2fc_col_label='log2FoldChange',
        Rank_Metric_S_col_label='Rank_Metric_S',
        sort_local=False,
        save_output=False,
        output_dir=None,
        return_df=False):
    '''
    This function adds a local statistic add_Rank_Metric_S to the deseq2 file
    Parameters:
        deseq2_file_path (str): Path to the DESeq2 result file.
        pvalue_col_label (str): Column label for p-value in the DESeq2 result file.
        l2fc_col_label (str): Column label for log2 fold change in the DESeq2 result file.
        Rank_Metric_S_col_label (str): Column label for the new Rank Metric S to be added.
        save_output (bool): If True, saves the modified DataFrame back to the original file path.
        output_dir (str): Directory to save the output deseq2 files. If None, saves to nfcore_DA_output_dir/tables/differential.
        return_df (bool): If True, returns the modified DataFrame.
    Returns: 
        pd.DataFrame or None: The modified DataFrame if return_df is True, otherwise None.
    '''
    import numpy as np
    import pandas as pd
    deseq2_df=pd.read_csv(deseq2_file_path)
    deseq2_df[Rank_Metric_S_col_label]=-np.log10(deseq2_df[pvalue_col_label])*(np.sign(deseq2_df[l2fc_col_label]))
    if sort_local:
        deseq2_df=deseq2_df.sort_values(by=Rank_Metric_S_col_label,ascending=False)
    # ensure that the Rank_Metric_S_col_label is inserted following the pvalue_col_label column in the dataframe
    cols = list(deseq2_df.columns)
    if Rank_Metric_S_col_label in cols:
        cols.remove(Rank_Metric_S_col_label)
    pvalue_col_index = cols.index(pvalue_col_label)
    cols.insert(pvalue_col_index + 1, Rank_Metric_S_col_label)
    deseq2_df = deseq2_df[cols]
    #deseq2_df=deseq2_df[list(deseq2_df.columns[:pvalue_col_index+1])+[Rank_Metric_S_col_label]+list(deseq2_df.columns[pvalue_col_index+1:-1])]
    if save_output:
        final_output_dir = output_dir if output_dir is not None else os.path.dirname(deseq2_file_path)
        deseq2_file_path=os.path.join(final_output_dir,os.path.basename(deseq2_file_path))
        deseq2_df.to_csv(deseq2_file_path,index=False)
        print(f'Rank_Metric_S added and saved to {deseq2_file_path}')
    if return_df:
        return deseq2_df
    return None

def add_Rank_Metric_S_to_all_deseq2_tables_nfcore_differentialabundance(
        nfcore_DA_output_dir='differentialabundance/results/',
        input_file_suffix='.deseq2.results.csv',
        output_dir=None,
        sort_local=False,):
    '''
    This function adds a local statistic add_Rank_Metric_S to all the deseq2 files in the 
        <nfcore_DA_output_dir>/differentialabundance/tables/differential directory
    Parameters:
        nfcore_DA_output_dir (str): Directory containing the output of nfcore/differentialabundance DESeq2 results files.
        input_file_suffix (str): Suffix for the input DESeq2 result files to be concatenated.
        sort_local (bool) : if true sorts individual tables by the rank metric column with ascending=False
    Returns:
        None
    '''
    import os
    import pandas as pd
    # #) set paths
    deseq_output_dir_tables = os.path.join(nfcore_DA_output_dir, 'tables', )
    deseq2_tables_dir_differential=os.path.join(deseq_output_dir_tables,'differential')
    #deseq2_tables_dir_processed_abundance=os.path.join(deseq_output_dir_tables,'processed_abundance')
    # ###############   loop to add Rank_Metric_S to all the files and save them as .csvs
    # add try except block to catch errors
    try:
        # list of all the files in the deseq2_tables_dir_differential
        deseq2_tables_dir_differential_files=os.listdir(deseq2_tables_dir_differential)
        for file in deseq2_tables_dir_differential_files:
            if file.endswith('.deseq2.results.csv'):
                file_path=os.path.join(deseq2_tables_dir_differential,file)
                add_Rank_Metric_S(file_path,
                                  sort_local=sort_local,
                                  save_output=True,
                                  )
    except Exception as e:
        print('error in add_Rank_Metric_S_to_all_deseq2_tables:', e) 
    return


####  function to  concatenate all the deseq2 tables in the deseq2 out put directory with addtiona columns to label which comparions

import os
import glob
import pandas as pd

def concat_deseq2_files_nfcore_differentialabundance(
    nfcore_DA_output_dir='differentialabundance/results/',
    non_nfcore_DA_output_dir=None,
    input_file_suffix='.deseq2.results.csv',
    output_file_prefix='differentialabundance',
    save_output=False,
    output_dir=None,
    also_save_to_src_dir=True,
    ):
    '''
    # run last
    Concatenates DESeq2 results files in a specified directory, optionally saving the output.
    Parameters:
        nfcore_DA_output_dir (str): Directory containing the output of nfcore/differentialabundance DESeq2 results files.
        input_file_suffix (str): Suffix for the input DESeq2 result files to be concatenated.
        output_file_prefix (str): Prefix for the output file name. concatenated file will be named '{output_file_prefix}_concat_deseq2.csv'.
        save_output (bool): If True, saves the concatenated DataFrame to a CSV file.
        output_dir (str): Directory to save the output CSV file. If None, saves to nfcore_DA_output_dir/tables/differential.
    Returns:
        pd.DataFrame: A concatenated DataFrame of all DESeq2 results.
    '''
    # #) set paths
    if non_nfcore_DA_output_dir is not None:
        deseq_output_dir = non_nfcore_DA_output_dir
    else:
        deseq_output_dir = os.path.join(nfcore_DA_output_dir, 'tables', 'differential')
    file_name=output_file_prefix+'_concat_deseq2.csv'
    final_output_dir = output_dir if output_dir is not None else deseq_output_dir
    output_filename = os.path.join(final_output_dir, file_name)
    # Construct the full path for DESeq2 result files
    files_pattern = os.path.join(deseq_output_dir, '*'+input_file_suffix)
    files = glob.glob(files_pattern)
    # Read and concatenate all DESeq2 result files
    df_list = []
    for file in files:
        df = pd.read_csv(file)
        # Extract base filename without extension and use as a column to identify the source file
        base_name = os.path.basename(file).replace(input_file_suffix, '')
        df['Deseq2_Comparison'] = base_name
        df_list.append(df)
    # Concatenate all dataframes into one
    df_concat = pd.concat(df_list, ignore_index=True)
    if save_output:
        df_concat.to_csv(output_filename, index=False)
        print(f'Concatenated DESeq2 results saved to {output_filename}')
        if also_save_to_src_dir and output_dir is not None:
            output_filename_src=os.path.join(nfcore_DA_output_dir, 'tables', 'differential', file_name)
            df_concat.to_csv(output_filename_src, index=False)
            print(f'Concatenated DESeq2 results also saved to source dir {output_filename_src}')
            
    return df_concat

#TODO finish this function 
def extract_read_quantification_metrics_nfcore_rnaseq(
        nfcore_output_rnaseq_dir='rnaseq/results/',
        save_output=False,
        output_dir=None,
        ):
    '''
    This function extracts the read quantification metrics from the nfcore rnaseq output directory
    parameters:
        nfcore_output_rnaseq_dir (str): Path to the nfcore rnaseq output directory.
        save_output (bool): Whether to save the read quantification metrics.
        output_dir (str): Directory to save the read quantification metrics. If None, saves to nfcore_output_rnaseq_dir/star_salmon.
        expected files:
            salmon.merged.gene_counts.tsv
            salmon.merged.gene_lengths.tsv
            salmon.merged.gene_counts_scaled.tsv
            salmon.merged.gene_counts_length_scaled.tsv
            salmon.merged.transcript_counts.tsv
            salmon.merged.transcript_lengths.tsv
            salmon.merged.gene_tpm.tsv
            salmon.merged.transcript_tpm.tsv
        note: normalized_expression_* aliases are intentionally not included.
    returns:
        path_dict (dict): Dictionary with paths to the read quantification metrics files.
        If save_output is True, saves the metrics to the specified output directory using the same filenames.
    '''
    import os
    import shutil

    star_salmon_dir = os.path.join(nfcore_output_rnaseq_dir, 'star_salmon')
    if save_output and output_dir is None:
        output_dir = star_salmon_dir

    file_mappings = [
        ('salmon.merged.gene_counts.tsv', 'salmon.merged.gene_counts.tsv'),
        ('salmon.merged.gene_lengths.tsv', 'salmon.merged.gene_lengths.tsv'),
        ('salmon.merged.gene_counts_scaled.tsv', 'salmon.merged.gene_counts_scaled.tsv'),
        ('salmon.merged.gene_counts_length_scaled.tsv', 'salmon.merged.gene_counts_length_scaled.tsv'),
        ('salmon.merged.transcript_counts.tsv', 'salmon.merged.transcript_counts.tsv'),
        ('salmon.merged.transcript_lengths.tsv', 'salmon.merged.transcript_lengths.tsv'),
        ('salmon.merged.gene_tpm.tsv', 'salmon.merged.gene_tpm.tsv'),
        ('salmon.merged.transcript_tpm.tsv', 'salmon.merged.transcript_tpm.tsv'),
    ]

    if save_output and output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)

    path_dict = {}
    for src_name, dest_name in file_mappings:
        src_path = os.path.join(star_salmon_dir, src_name)
        if not os.path.isfile(src_path):
            print(f'[WARN] Read quantification file missing: {src_path}')
            continue

        if save_output and output_dir is not None:
            dest_path = os.path.join(output_dir, dest_name)
            if os.path.abspath(dest_path) != os.path.abspath(src_path):
                shutil.copy2(src_path, dest_path)
            path_dict[dest_name] = dest_path
        else:
            path_dict[dest_name] = src_path

    return path_dict

#### Start ################################ Standard plots the NFcore outputs. ####################################


################################ ## function for Make plots for each NFcore directory 
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings(
    "ignore",
    message=r"The palette list has more values .*",
    category=UserWarning,
)
def make_volcano_MA_plots_from_nf_deseq_output_dir(
        deseq_dir,file_sufix='.deseq2.results.csv',
        plot_dataset_title=None,
        log2FoldChange_threshold=0.5,
        ylimit_volcano=None,
        xlimit_volcano_l2fc=None,
        xlimit_MAonly=None,
        label_gene_name=False,
        figsize=(12,12),
        save_plots=True,
        save_dir='plots'
        ):
    '''
    This function generates volcano and MA plots for each deseq2 results file in the specified nfcore output directory.
    It creates a directory for the plots and saves the plots as PNG files.,'''
    #path to specfic deseq2 run directory 
    #differeential_abundance_dir=os.path.join(nfcore_output_dir,deseq_dir)
    ## paths for deseq2 table results
    deseq2_tables_dir=deseq_dir
    list_of_deseq2_files_csv=[file for file in os.listdir(deseq2_tables_dir)  if file.endswith(file_sufix) ]
    # make new directory for non-nfcore plots in the differeential_abundance_dir directory
    plots_dir=os.path.join(deseq2_tables_dir,save_dir)
    # Create the directory, including any necessary parent directories
    os.makedirs(plots_dir, exist_ok=True)

    # print list of deseq2 files plots
    print(list_of_deseq2_files_csv)
    for deseq2_file in list_of_deseq2_files_csv:
        deseq2_file_path=os.path.join(deseq2_tables_dir,deseq2_file)
        print(deseq2_file_path)
        df=pd.read_csv(deseq2_file_path,)   
        # set comparison to be first part of file base name  
        comparison=os.path.splitext(os.path.splitext(os.path.splitext(os.path.basename(deseq2_file_path))[0])[0])[0]
        # try statement to catch error if volcano plot fails with e as error
        try:
            print(deseq2_file_path,'volcano plot')
            p=volcano_plot_sns_sinlge_comparison(df,
                    title_text=f'{plot_dataset_title}\n',
                    comparison=comparison,
                    sharex=True,sharey=True,
                    ylimit=ylimit_volcano,xlimit=xlimit_volcano_l2fc,
                    log2FoldChange_threshold=log2FoldChange_threshold,
                    figsize=figsize,
                    label_gene_name=label_gene_name
                    )
            if label_gene_name:
                volcano_plot_file_path=os.path.join(plots_dir,f'{comparison}.volcano.labeled.png')
            else:
                volcano_plot_file_path=os.path.join(plots_dir,f'{comparison}.volcano.png')
            if save_plots:
                p.figure.savefig(volcano_plot_file_path, bbox_inches='tight')
                plt.close(p.figure)
            else:
                plt.close(p.figure)
        except Exception as e:
            print(f'Error creating volcano plot for {comparison} \n {e}')

        # try statement to catch error if MA plot fails with e as error
        try:
            print(deseq2_file_path,'MA plot')
            ma=MA_plot_sns_sinlge_comparison(df,
                    title_text=f'{plot_dataset_title}\n',
                    comparison=comparison,
             sharex=True,sharey=True,
            ylimit=xlimit_volcano_l2fc,xlimit=xlimit_MAonly,
            log2FoldChange_threshold=log2FoldChange_threshold,
             figsize=figsize,
             label_gene_name=label_gene_name
             )
            if label_gene_name:
                ma_plot_file_path=os.path.join(plots_dir,f'{comparison}.MAplot.labeled.png')
            else:
                ma_plot_file_path=os.path.join(plots_dir,f'{comparison}.MAplot.png')
            if save_plots:
                ma.figure.savefig(ma_plot_file_path, bbox_inches='tight')
                plt.close(ma.figure)
            else:
                plt.close(ma.figure)
        except Exception as e:
            print(f'Error creating MA plot for {comparison} \n {e}')

        # try statement to catch error if volcano_MA plot fails with e as error
        try:
            print(deseq2_file_path,'volcano_MA plot')
            long_figsize = (figsize[0] * 1, figsize[1] * 2)
            p_ma=volcano_MA_plot_sns_single_comparison(df,
                                title_text=f'{plot_dataset_title}\n',
                                comparison=comparison,
                                sharex=False,sharey=False,
                                ylimit_volcano=ylimit_volcano,xlimit=xlimit_volcano_l2fc,
                                log2FoldChange_threshold=log2FoldChange_threshold,
                                figsize=long_figsize,
                                label_gene_name=label_gene_name
             )
            if label_gene_name:
                p_ma_plot_file_path=os.path.join(plots_dir,f'{comparison}.volcano_MAplot.labeled.png')
            else:
                p_ma_plot_file_path=os.path.join(plots_dir,f'{comparison}.volcano_MAplot.png')
            if save_plots:
                p_ma.figure.savefig(p_ma_plot_file_path, bbox_inches='tight')
                plt.close(p_ma.figure)
            else:
                plt.close(p_ma.figure)
        except Exception as e:
            print(f'Error creating MA plot for {comparison} \n {e}')


################################ volcano plot functions

import seaborn.objects as so



def volcano_plot_sns_from_concat(df, title_text='volcano_plot',comparison=None,
                     facet_col=None,dot_color=None,
                     sharex=True,sharey=True,ylimit=None,xlimit=None,
                     log2FoldChange_threshold=1,
                     figsize=(15, 10),label_gene_name=False):

    """
    Create a volcano plot using the given DataFrame.
    """
    import seaborn as sns
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.pyplot as plt

        # Get the tab10 palette
    tab10_palette = sns.color_palette("tab10")

    # Move gray to the first position
    #custom_palette = [tab10_palette[7]] + tab10_palette[:7] + tab10_palette[8:]
    custom_palette = [tab10_palette[7]] + tab10_palette[:3] #+ tab10_palette[8:]

    print(df.shape)
    # filter df to only include rows where 'Comperison' column equals comparison argument
    # filter df to only include rows where 'Comperison' column equals comparison argument
    if comparison:
        df = df[df['Deseq2_Comparison']==comparison].copy()
        # remove unused categories for each column in df that is a categorical
        # Assuming 'df' is your DataFrame
        # Get the list of categorical columns
        categorical_columns = df.select_dtypes(include='category').columns
        # Remove unused categories for each categorical column
        for column in categorical_columns:
            df[column] = df[column].cat.remove_unused_categories()
        print(df.shape)
        

        

    required_columns = {'log2FoldChange', 'padj', 'Deseq2_Comparison'}
    if not required_columns.issubset(df.columns):
        raise ValueError(f"DataFrame is missing one of the required columns: {required_columns}")

    ###### Prepare df by adding coluomns for '-log10(padj)' and signifigcance level and adjusting out of range values
    # Replace NaN values with a specified value, for example, 1
    df['padj'] = df['padj'].fillna(1)
    # Prepare data by adjusting p-values to avoid log(0) issues
    df['-log10(padj)'] = -np.log10(df['padj'].replace(0, np.nextafter(0, 1)))
    # Assuming df['padj'] contains the adjusted p-values
    # Add column for alpha=0.2 significance level
    df['alpha=0.2'] = ((df['padj'] < 0.2) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))
    # Add column for alpha=0.1 significance level
    df['alpha=0.1'] = ((df['padj'] < 0.1) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))
    # Add column for alpha=0.05 significance level
    df['alpha=0.05'] = ((df['padj'] < 0.05) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))

    # add column for signifigcance hue
    # first combine the alpha columns into one column labeled Signifigcance
    #df['Significance'] = 'alpha>0.2'
    df['Significance'] = 'Not Significant'
    df.loc[df['alpha=0.2'],'Significance'] = 'alpha=0.2'
    df.loc[df['alpha=0.1'],'Significance'] = 'alpha=0.1'
    df.loc[df['alpha=0.05'],'Significance'] = 'alpha=0.05'
    df['Significance'] = df['Significance'].astype('category')



    ######  adjusting out of range values and changing dot type if out of range
    ##### #####  set limits 
    # set ylimit if none to  to 99 percentile of ['-log10(padj)']
    if not ylimit:
        ylimit = df[(df['padj']<0.05)&(df['log2FoldChange'].abs()>log2FoldChange_threshold)]['-log10(padj)'].quantile(0.99)
        if np.isnan(ylimit):
            ylimit=df['-log10(padj)'].quantile(0.99)
    # set xlimit if none to 99 percentile of abs(x)
    if not xlimit:
        xlimit = df[(df['padj']<0.05)&(df['log2FoldChange'].abs()>log2FoldChange_threshold)]['log2FoldChange'].abs().quantile(0.99)
    # if xlimit is nan set to quantile(0.99)
        if np.isnan(xlimit):
            xlimit=df['log2FoldChange'].abs().quantile(0.99)

    # add 'Marker' column for out of range data points with '-log10(padj)' > ylimit or abs('log2FoldChange') > xlimit value of 'In_Range' or 'Out_of_Range'
    df['Marker'] = 'In_Range'
    df.loc[df['-log10(padj)']>=ylimit,'Marker'] = 'Out_of_Range'
    #  abs('log2FoldChange') 
    df.loc[abs(df['log2FoldChange'])>=xlimit,'Marker'] = 'Out_of_Range'
    # order the categories
    # Ensure the required categories are present
    df['Marker'] = df['Marker'].astype('category')
    required_range_cats = {'In_Range', 'Out_of_Range'}
    if not required_range_cats.issubset(df['Marker'].cat.categories):
        df['Marker'] = df['Marker'].cat.set_categories(['In_Range', 'Out_of_Range'])


    # replace values in the -log10(padj) column that above the ylimit with the ylimit
    if ylimit:
        df['-log10(padj)'] = df['-log10(padj)'].apply(lambda x: (ylimit*0.99) if x>=ylimit else x)
    else:
        ylimit = df['-log10(padj)'].max()
    # replace values in the log2FoldChange column that above or below the xlimit with the xlimit
    if xlimit:
        df['log2FoldChange'] = df['log2FoldChange'].apply(lambda x: (xlimit*0.99) if x>=xlimit else x)
        df['log2FoldChange'] = df['log2FoldChange'].apply(lambda x: (-xlimit*0.99)  if x<=-xlimit else x)
    else:
        xlimit = max(abs(df['log2FoldChange'].min()), df['log2FoldChange'].max())

    ### set the marker size relative to number of dots
    rel_size=df.shape[0]/333






    if label_gene_name:
        fig, ax = plt.subplots(figsize=figsize)
        p = sns.scatterplot(data=df, x='log2FoldChange', y='-log10(padj)', hue='Significance', style='Marker', 
                            palette=custom_palette,sizes=(rel_size),  s=rel_size, 
                            ax=ax)
        p.set(xlim=(-xlimit, xlimit), ylim=(0, ylimit))
        p.set_title(f'{title_text}\n{comparison}\n\n')
        p.set_xlabel("log2fc Deseq2 model")
        p.set_ylabel("-log10(padj) Deseq2 model")
        p.legend(title=facet_col)
        # move legend
        plt.legend(bbox_to_anchor=(1.25, 1), 
            loc=1, 
            borderaxespad=0.)
        #label top genes by padj
        for line in range(0,50):
            p.text(df.sort_values(by='padj')['log2FoldChange'].to_list()[line],df.sort_values(by='padj')['-log10(padj)'].to_list()[line],
                   df.sort_values(by='padj').gene_name.to_list()[line],
                      horizontalalignment='left', size='small', color='black')
    else:
        fig, ax = plt.subplots(figsize=figsize)
        p = sns.scatterplot(data=df, x='log2FoldChange', y='-log10(padj)',hue='Significance', style='Marker', 
                             palette=custom_palette, s=rel_size,  
                            ax=ax)
        p.set(xlim=(-xlimit, xlimit), ylim=(0, ylimit))
        p.set_title(f'{title_text}\n{comparison}\n\n')
        p.set_xlabel("log2fc Deseq2 model")
        p.set_ylabel("-log10(padj) Deseq2 model")
        p.legend(title=facet_col)
        # move legend
        plt.legend(bbox_to_anchor=(1.22, 1), 
            loc=1, 
            borderaxespad=0.)


    
    return p


def volcano_plot_sns_sinlge_comparison(df, title_text='volcano_plot',comparison='DeSeq2 Comparison',
                     facet_col=None,dot_color=None,
                     sharex=True,sharey=True,ylimit=None,xlimit=None,
                     log2FoldChange_threshold=1.0,
                     figsize=(15, 10),label_gene_name=False):

    """
    Create a volcano plot using the given DataFrame.
    """
    import seaborn as sns
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.pyplot as plt

        # Get the tab10 palette
    tab10_palette = sns.color_palette("tab10")

    # Move gray to the first position
    #custom_palette = [tab10_palette[7]] + tab10_palette[:7] + tab10_palette[8:]
    custom_palette = [tab10_palette[7]] + tab10_palette[:3] #+ tab10_palette[8:]

    print(df.shape)

    required_columns = {'log2FoldChange', 'padj', }
    if not required_columns.issubset(df.columns):
        raise ValueError(f"DataFrame is missing one of the required columns: {required_columns}")

    ###### Prepare df by adding coluomns for '-log10(padj)' and signifigcance level and adjusting out of range values
    # Replace NaN values with a specified value, for example, 1
    df['padj'] = df['padj'].fillna(1)
    # Prepare data by adjusting p-values to avoid log(0) issues
    df['-log10(padj)'] = -np.log10(df['padj'].replace(0, np.nextafter(0, 1)))
    # Assuming df['padj'] contains the adjusted p-values
    # Add column for alpha=0.2 significance level
    df['alpha=0.2'] = ((df['padj'] < 0.2) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))
    # Add column for alpha=0.1 significance level
    df['alpha=0.1'] = ((df['padj'] < 0.1) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))
    # Add column for alpha=0.05 significance level
    df['alpha=0.05'] = ((df['padj'] < 0.05) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))

    # add column for signifigcance hue
    # first combine the alpha columns into one column labeled Signifigcance
    #df['Significance'] = 'alpha>0.2'
    df['Significance'] = 'Not Significant'
    df.loc[df['alpha=0.2'],'Significance'] = 'alpha=0.2'
    df.loc[df['alpha=0.1'],'Significance'] = 'alpha=0.1'
    df.loc[df['alpha=0.05'],'Significance'] = 'alpha=0.05'
    df['Significance'] = df['Significance'].astype('category')



    ######  adjusting out of range values and changing dot type if out of range


    ##### #####  set limits 
    # set ylimit if none to  to 99 percentile of ['-log10(padj)']
    if not ylimit:
        ylimit = df[(df['padj']<0.05)&(df['log2FoldChange'].abs()>log2FoldChange_threshold)]['-log10(padj)'].quantile(0.99)
        if np.isnan(ylimit):
            ylimit=df['-log10(padj)'].quantile(0.99)
    # set xlimit if none to 99 percentile of abs(x)
    if not xlimit:
        xlimit = df[(df['padj']<0.05)&(df['log2FoldChange'].abs()>log2FoldChange_threshold)]['log2FoldChange'].abs().quantile(0.99)
    # if xlimit is nan set to quantile(0.99)
        if np.isnan(xlimit):
            xlimit=df['log2FoldChange'].abs().quantile(0.99)


    # add 'Marker' column for out of range data points with '-log10(padj)' > ylimit or abs('log2FoldChange') > xlimit value of 'In_Range' or 'Out_of_Range'
    df['Marker'] = 'In_Range'
    df.loc[df['-log10(padj)']>=ylimit,'Marker'] = 'Out_of_Range'
    #  abs('log2FoldChange') 
    df.loc[abs(df['log2FoldChange'])>=xlimit,'Marker'] = 'Out_of_Range'

    # order the categories
    # Ensure the required categories are present
    df['Marker'] = df['Marker'].astype('category')
    required_range_cats = {'In_Range', 'Out_of_Range'}
    if not required_range_cats.issubset(df['Marker'].cat.categories):
        df['Marker'] = df['Marker'].cat.set_categories(['In_Range', 'Out_of_Range'])


    # replace values in the -log10(padj) column that above the ylimit with the ylimit
    if ylimit:
        df['-log10(padj)'] = df['-log10(padj)'].apply(lambda x: (ylimit*0.99) if x>=ylimit else x)
    else:
        ylimit = df['-log10(padj)'].max()
    # replace values in the log2FoldChange column that above or below the xlimit with the xlimit
    if xlimit:
        df['log2FoldChange'] = df['log2FoldChange'].apply(lambda x: (xlimit*0.99) if x>=xlimit else x)
        df['log2FoldChange'] = df['log2FoldChange'].apply(lambda x: (-xlimit*0.99)  if x<=-xlimit else x)
    else:
        xlimit = max(abs(df['log2FoldChange'].min()), df['log2FoldChange'].max())

    ### set the marker size relative to number of dots
    rel_size=df.shape[0]/333


    if label_gene_name:
        fig, ax = plt.subplots(figsize=figsize)
        p = sns.scatterplot(data=df, x='log2FoldChange', y='-log10(padj)', hue='Significance', style='Marker', 
                            palette=custom_palette,sizes=(rel_size),  s=rel_size, 
                            ax=ax)
        p.set(xlim=(-xlimit, xlimit), ylim=(0, ylimit))
        p.set_title(f'{title_text}\n{comparison}\n\n')
        p.axvline(x=log2FoldChange_threshold, color='gray', linestyle='--',label=f'log2fc>|{log2FoldChange_threshold}| ')
        p.axvline(x=-log2FoldChange_threshold, color='gray', linestyle='--')
        p.set_xlabel("log2fc Deseq2 model")
        p.set_ylabel("-log10(padj) Deseq2 model")
        p.legend( #title=facet_col,
            bbox_to_anchor=(1.25, 1), 
            loc=1, 
            borderaxespad=0.)
        #label top genes by padj
        for line in range(0,50):
            p.text(df.sort_values(by='padj')['log2FoldChange'].to_list()[line],df.sort_values(by='padj')['-log10(padj)'].to_list()[line],
                   df.sort_values(by='padj').gene_name.to_list()[line],
                      horizontalalignment='left', size='small', color='black')
        #for line in range(0,df.shape[0]):
         #   p.text(df.log2FoldChange[line], df['-log10(padj)'][line], df.gene_name[line], horizontalalignment='left', size='medium', color='black')
    else:
        fig, ax = plt.subplots(figsize=figsize)
        p = sns.scatterplot(data=df, x='log2FoldChange', y='-log10(padj)',hue='Significance', style='Marker', 
                             palette=custom_palette, s=rel_size,  
                            ax=ax)
        p.set(xlim=(-xlimit, xlimit), ylim=(0, ylimit))
        p.set_title(f'{title_text}\n{comparison}\n\n')
        p.axvline(x=log2FoldChange_threshold, color='gray', linestyle='--',label=f'log2fc>|{log2FoldChange_threshold}| ')
        p.axvline(x=-log2FoldChange_threshold, color='gray', linestyle='--',)
        p.set_xlabel("log2fc Deseq2 model")
        p.set_ylabel("-log10(padj) Deseq2 model")
        p.legend( #title=facet_col
            )
        # move legend
        plt.legend(bbox_to_anchor=(1.25, 1), 
            loc=1, 
            borderaxespad=0.)

    return p


################################ MA plot functions



def MA_plot_sns_from_concat(df, title_text='MA_plot',comparison='DeSeq2 Comparison',
                     facet_col=None,dot_color=None,
                     sharex=True,sharey=True,ylimit=None,xlimit=None,
                     log2FoldChange_threshold=1,
                     figsize=(15, 10),label_gene_name=False):

    """
    Create a volcano plot using the given DataFrame.
    """
    import seaborn as sns
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.pyplot as plt

        # Get the tab10 palette
    tab10_palette = sns.color_palette("tab10")

    # Move gray to the first position
    #custom_palette = [tab10_palette[7]] + tab10_palette[:7] + tab10_palette[8:]
    custom_palette = [tab10_palette[7]] + tab10_palette[:3] #+ tab10_palette[8:]

    print(df.shape)
    # filter df to only include rows where 'Comperison' column equals comparison argument
    # filter df to only include rows where 'Comperison' column equals comparison argument
    if comparison:
        df = df[df['Deseq2_Comparison']==comparison].copy()
        # remove unused categories for each column in df that is a categorical
        # Assuming 'df' is your DataFrame
        # Get the list of categorical columns
        categorical_columns = df.select_dtypes(include='category').columns
        # Remove unused categories for each categorical column
        for column in categorical_columns:
            df[column] = df[column].cat.remove_unused_categories()
        print(df.shape)

    required_columns = {'log2FoldChange', 'padj','baseMean', 'Deseq2_Comparison'}
    if not required_columns.issubset(df.columns):
        raise ValueError(f"DataFrame is missing one of the required columns: {required_columns}")

    ###### Prepare df by adding coluomns for 'log2(baseMean)' and signifigcance level and adjusting out of range values
    df['log2(baseMean)'] = np.log2(df['baseMean'].replace(0, np.nextafter(0, 1)))
    # Calculate log2(baseMean) and ensure non-negative values
    df['log2(baseMean)'] = df['log2(baseMean)'].clip(lower=0)

    # Assuming df['padj'] contains the adjusted p-values
    # Add column for alpha=0.2 significance level
    df['alpha=0.2'] = ((df['padj'] < 0.2) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))
    # Add column for alpha=0.1 significance level
    df['alpha=0.1'] = ((df['padj'] < 0.1) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))
    # Add column for alpha=0.05 significance level
    df['alpha=0.05'] = ((df['padj'] < 0.05) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))

    # add column for signifigcance hue
    # first combine the alpha columns into one column labeled Signifigcance
    #df['Significance'] = 'alpha>0.2'
    df['Significance'] = 'Not Significant'
    df.loc[df['alpha=0.2'],'Significance'] = 'alpha=0.2'
    df.loc[df['alpha=0.1'],'Significance'] = 'alpha=0.1'
    df.loc[df['alpha=0.05'],'Significance'] = 'alpha=0.05'
    df['Significance'] = df['Significance'].astype('category')


    ################## START ######   adjusting out of range values and changing dot type if out of range
    # set ylimit if none to 99 percentile of abs(x)
    if not ylimit:
        ylimit = df[(df['padj']<0.05)&(df['log2FoldChange'].abs()>log2FoldChange_threshold)]['log2FoldChange'].abs().quantile(0.99)
        if np.isnan(ylimit):
            ylimit=df['log2FoldChange'].quantile(0.99)
    # set xlimit if none to  to 99 percentile of 
    if not xlimit:
        #xlimit = df[df['padj']<0.05]['log2(baseMean)'].quantile(0.99)
        xlimit = df[(df['padj']<0.05)&(df['log2FoldChange'].abs()>log2FoldChange_threshold)]['log2(baseMean)'].quantile(0.99)
        # if xlimit is nan set to quantile(0.99)
        if np.isnan(xlimit):
            xlimit=df['log2(baseMean)'].abs().quantile(0.99)

    # add 'Marker' column for out of range data points with log2(baseMean) > ylimit or abs('log2FoldChange') > xlimit value of 'In_Range' or 'Out_of_Range'
    df['Marker'] = 'In_Range'
    df.loc[df['log2(baseMean)']>=xlimit,'Marker'] = 'Out_of_Range'
    #  abs('log2FoldChange') 
    df.loc[abs(df['log2FoldChange'])>=ylimit,'Marker'] = 'Out_of_Range'
    # order the categories
    # Ensure the required categories are present
    df['Marker'] = df['Marker'].astype('category')
    required_range_cats = {'In_Range', 'Out_of_Range'}
    if not required_range_cats.issubset(df['Marker'].cat.categories):
        df['Marker'] = df['Marker'].cat.set_categories(['In_Range', 'Out_of_Range'])


    # replace values in the log2FoldChange column that above the ylimit with the ylimit
    if xlimit:
         df['log2(baseMean)'] = df['log2(baseMean)'].apply(lambda x: (xlimit*0.99) if x>=xlimit else x)
    # replace values in the log2FoldChange column that above or below the ylimit with the ylimit
    if ylimit:
        df['log2FoldChange'] = df['log2FoldChange'].apply(lambda x: (ylimit*0.99) if x>=ylimit else x)
        df['log2FoldChange'] = df['log2FoldChange'].apply(lambda x: (-ylimit*0.99)  if x<=-ylimit else x)

    ################## END ######   adjusting out of range values and changing dot type if out of range


    ### set the marker size relative to number of dots
    rel_size=df.shape[0]/333

    if label_gene_name:
        fig, ax = plt.subplots(figsize=figsize)
        p = sns.scatterplot(data=df, x='log2(baseMean)', y='log2FoldChange', hue='Significance', style='Marker', 
                            palette=custom_palette,sizes=(rel_size),  s=rel_size, 
                            ax=ax)
        p.set(xlim=(0, xlimit), ylim=(-ylimit, ylimit))
        p.set_title(f'{title_text}\n{comparison}\n\n')
        p.set_ylabel("log2fc Deseq2 model")
        p.set_xlabel("log2(baseMean) Deseq2 model")
        p.legend(#title=facet_col,
            bbox_to_anchor=(1.25, 1), 
            loc=1, 
            borderaxespad=0.)
        #label top genes by padj
        for line in range(0,50):
            p.text(df.sort_values(by='padj')['log2(baseMean)'].to_list()[line],df.sort_values(by='padj')['log2FoldChange'].to_list()[line],
                   df.sort_values(by='padj').gene_name.to_list()[line],
                      horizontalalignment='left', size='small', color='black')
        #for line in range(0,df.shape[0]):
        #    p.text(df.log2FoldChange[line], df['-log10(padj)'][line], df.gene_name[line], horizontalalignment='left', size='medium', color='black')
    else:
        fig, ax = plt.subplots(figsize=figsize)
        p = sns.scatterplot(data=df, x='log2(baseMean)', y='log2FoldChange',hue='Significance', style='Marker', 
                             palette=custom_palette, s=rel_size,  
                            ax=ax)
        p.set(xlim=(0, xlimit), ylim=(-ylimit, ylimit))
        p.set_title(f'{title_text}\n{comparison}\n\n')
        p.set_ylabel("log2fc Deseq2 model")
        p.set_xlabel("log2(baseMean) Deseq2 model")
        p.legend(#title=facet_col,
            bbox_to_anchor=(1.25, 1), 
            loc=1, 
            borderaxespad=0.)


    return p


def MA_plot_sns_sinlge_comparison(df, title_text='MA_plot',comparison='DeSeq2 Comparison',
                     facet_col=None,dot_color=None,
                     sharex=True,sharey=True,ylimit=None,xlimit=None,
                     log2FoldChange_threshold=1.0,
                     figsize=(15, 10),label_gene_name=False):

    """
    Create a volcano plot using the given DataFrame.
    """
    import seaborn as sns
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.pyplot as plt

        # Get the tab10 palette
    tab10_palette = sns.color_palette("tab10")

    # Move gray to the first position
    #custom_palette = [tab10_palette[7]] + tab10_palette[:7] + tab10_palette[8:]
    custom_palette = [tab10_palette[7]] + tab10_palette[:3] #+ tab10_palette[8:]

    print(df.shape)

    required_columns = {'log2FoldChange', 'padj','baseMean', }
    if not required_columns.issubset(df.columns):
        raise ValueError(f"DataFrame is missing one of the required columns: {required_columns}")

    ###### Prepare df by adding coluomns for 'log2(baseMean)' and signifigcance level and adjusting out of range values

    df['log2(baseMean)'] = np.log2(df['baseMean'].replace(0, np.nextafter(0, 1)))
    # Calculate log2(baseMean) and ensure non-negative values
    df['log2(baseMean)'] = df['log2(baseMean)'].clip(lower=0)

    # Assuming df['padj'] contains the adjusted p-values
    # Add column for alpha=0.2 significance level
    df['alpha=0.2'] = ((df['padj'] < 0.2) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))
    # Add column for alpha=0.1 significance level
    df['alpha=0.1'] = ((df['padj'] < 0.1) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))
    # Add column for alpha=0.05 significance level
    df['alpha=0.05'] = ((df['padj'] < 0.05) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))

    # add column for signifigcance hue
    # first combine the alpha columns into one column labeled Signifigcance
    #df['Significance'] = 'alpha>0.2'
    df['Significance'] = 'Not Significant'
    df.loc[df['alpha=0.2'],'Significance'] = 'alpha=0.2'
    df.loc[df['alpha=0.1'],'Significance'] = 'alpha=0.1'
    df.loc[df['alpha=0.05'],'Significance'] = 'alpha=0.05'
    df['Significance'] = df['Significance'].astype('category')



    ################## START ######   adjusting out of range values and changing dot type if out of range
    # set ylimit if none to 99 percentile of abs(x)
    if not ylimit:
        ylimit = df[(df['padj']<0.05)&(df['log2FoldChange'].abs()>log2FoldChange_threshold)]['log2FoldChange'].abs().quantile(0.99)
        if np.isnan(ylimit):
            ylimit=df['log2FoldChange'].quantile(0.99)
    # set xlimit if none to  to 99 percentile of 
    if not xlimit:
        #xlimit = df[df['padj']<0.05]['log2(baseMean)'].quantile(0.99)
        xlimit = df[(df['padj']<0.05)&(df['log2FoldChange'].abs()>log2FoldChange_threshold)]['log2(baseMean)'].quantile(0.99)
        # if xlimit is nan set to quantile(0.99)
        if np.isnan(xlimit):
            xlimit=df['log2(baseMean)'].abs().quantile(0.99)

    # add 'Marker' column for out of range data points with log2(baseMean) > ylimit or abs('log2FoldChange') > xlimit value of 'In_Range' or 'Out_of_Range'
    df['Marker'] = 'In_Range'
    df.loc[df['log2(baseMean)']>=xlimit,'Marker'] = 'Out_of_Range'
    #  abs('log2FoldChange') 
    df.loc[abs(df['log2FoldChange'])>=ylimit,'Marker'] = 'Out_of_Range'

    # order the categories
    # Ensure the required categories are present
    df['Marker'] = df['Marker'].astype('category')
    required_range_cats = {'In_Range', 'Out_of_Range'}
    if not required_range_cats.issubset(df['Marker'].cat.categories):
        df['Marker'] = df['Marker'].cat.set_categories(['In_Range', 'Out_of_Range'])


    # replace values in the log2FoldChange column that above the ylimit with the ylimit
    if xlimit:
         df['log2(baseMean)'] = df['log2(baseMean)'].apply(lambda x: (xlimit*0.99) if x>=xlimit else x)
    # replace values in the log2FoldChange column that above or below the ylimit with the ylimit
    if ylimit:
        df['log2FoldChange'] = df['log2FoldChange'].apply(lambda x: (ylimit*0.99) if x>=ylimit else x)
        df['log2FoldChange'] = df['log2FoldChange'].apply(lambda x: (-ylimit*0.99)  if x<=-ylimit else x)

    ################## END ######   adjusting out of range values and changing dot type if out of range


    ### set the marker size relative to number of dots
    rel_size=df.shape[0]/333

    if label_gene_name:
        fig, ax = plt.subplots(figsize=figsize)
        p = sns.scatterplot(data=df, x='log2(baseMean)', y='log2FoldChange', hue='Significance', style='Marker', 
                            palette=custom_palette,sizes=(rel_size),  s=rel_size, 
                            ax=ax)
        # add verticle line at y=xpercentile
        percentile_line_BM=df['log2(baseMean)'].quantile([0.7,]).values[0]
        percentile_line_BM_counts=df['baseMean'].quantile([0.7,]).values[0]
        ax.axvline(x=percentile_line_BM, color='black', linestyle=':',
                   label=f'70th percentile \n  counts= {int(percentile_line_BM_counts)} ')
        p.axhline(y=log2FoldChange_threshold, color='gray', linestyle='--',label=f'log2fc>|{log2FoldChange_threshold}| ')
        p.axhline(y=-log2FoldChange_threshold, color='gray', linestyle='--',)
        p.set(xlim=(0, xlimit), ylim=(-ylimit, ylimit))
        p.set_title(f'{title_text}\n{comparison}\n\n')
        p.set_ylabel("log2fc Deseq2 model")
        p.set_xlabel("log2(baseMean) Deseq2 model")
        p.legend(bbox_to_anchor=(1.25, 1), 
            loc=1, 
            borderaxespad=0.)
        #label top genes by padj
        for line in range(0,50):
            p.text(df.sort_values(by='padj')['log2(baseMean)'].to_list()[line],df.sort_values(by='padj')['log2FoldChange'].to_list()[line],
                   df.sort_values(by='padj').gene_name.to_list()[line],
                      horizontalalignment='left', size='small', color='black')
    else:
        fig, ax = plt.subplots(figsize=figsize)
        p = sns.scatterplot(data=df, x='log2(baseMean)', y='log2FoldChange',hue='Significance', style='Marker', 
                             palette=custom_palette, s=rel_size,  
                            ax=ax)
        # add verticle line at y=xpercentile
        percentile_line_BM=df['log2(baseMean)'].quantile([0.7,]).values[0]
        percentile_line_BM_counts=df['baseMean'].quantile([0.7,]).values[0]
        ax.axvline(x=percentile_line_BM, color='black', linestyle=':',
                   label=f'70th percentile \n  counts= {int(percentile_line_BM_counts)} ')
        p.axhline(y=log2FoldChange_threshold, color='gray', linestyle='--',label=f'log2fc>|{log2FoldChange_threshold}| ')
        p.axhline(y=-log2FoldChange_threshold, color='gray', linestyle='--',)
        p.set(xlim=(0, xlimit), ylim=(-ylimit, ylimit))
        p.set_title(f'{title_text}\n{comparison}\n\n')
        p.set_ylabel("log2fc Deseq2 model")
        p.set_xlabel("log2(baseMean) Deseq2 model")
        p.legend(bbox_to_anchor=(1.25, 1), 
            loc=1, 
            borderaxespad=0.)
    
    return p


################################ Volcano plot with MA plot function

def volcano_MA_plot_sns_single_comparison(df, title_text='volcano_plot',comparison='DeSeq2 Comparison',
                     facet_col=None,dot_color=None,
                     sharex=True,sharey=True,
                     ylimit_volcano=None,ylimit_MA=None,
                     xlimit=None,
                     log2FoldChange_threshold=1.0,
                     figsize=(12, 12),label_gene_name=False):

    """
    Create a volcano plot using the given DataFrame.
    """
    import seaborn as sns
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.pyplot as plt

        # Get the tab10 palette
    tab10_palette = sns.color_palette("tab10")

    # Move gray to the first position
    #custom_palette = [tab10_palette[7]] + tab10_palette[:7] + tab10_palette[8:]
    custom_palette = [tab10_palette[7]] + tab10_palette[:3] #+ tab10_palette[8:]

    print(df.shape)

    required_columns = {'log2FoldChange', 'padj','baseMean', }
    if not required_columns.issubset(df.columns):
        raise ValueError(f"DataFrame is missing one of the required columns: {required_columns}")

    ###### Prepare df by adding coluomns for 'log2(baseMean)' 
    df['log2(baseMean)'] = np.log2(df['baseMean'].replace(0, np.nextafter(0, 1)))
    # Calculate log2(baseMean) and ensure non-negative values
    df['log2(baseMean)'] = df['log2(baseMean)'].clip(lower=0)

    ###### Prepare df by adding coluomns for '-log10(padj)' and signifigcance level and adjusting out of range values
    # Replace NaN values with a specified value, for example, 1
    df['padj'] = df['padj'].fillna(1)
    # Prepare data by adjusting p-values to avoid log(0) issues
    df['-log10(padj)'] = -np.log10(df['padj'].replace(0, np.nextafter(0, 1)))
    # Assuming df['padj'] contains the adjusted p-values
    # Add column for alpha=0.2 significance level
    df['alpha=0.2'] = ((df['padj'] < 0.2) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))
    # Add column for alpha=0.1 significance level
    df['alpha=0.1'] = ((df['padj'] < 0.1) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))
    # Add column for alpha=0.05 significance level
    df['alpha=0.05'] = ((df['padj'] < 0.05) & (abs(df['log2FoldChange'])>=log2FoldChange_threshold))

    # add column for signifigcance hue
    # first combine the alpha columns into one column labeled Signifigcance
    #df['Significance'] = 'alpha>0.2'
    df['Significance'] = 'Not Significant'
    df.loc[df['alpha=0.2'],'Significance'] = 'alpha=0.2'
    df.loc[df['alpha=0.1'],'Significance'] = 'alpha=0.1'
    df.loc[df['alpha=0.05'],'Significance'] = 'alpha=0.05'
    df['Significance'] = df['Significance'].astype('category')



    ######  adjusting out of range values and changing dot type if out of range

    ##### #####  set limits 
    # set ylimit_MA if none to  to 99 percentile of ['log2(baseMean)']
    if not ylimit_MA:
        #ylimit_MA = df[df['padj']<0.05]['log2(baseMean)'].quantile(0.99)
        ylimit_MA = df[(df['padj']<0.05)&(df['log2FoldChange'].abs()>log2FoldChange_threshold)]['log2(baseMean)'].quantile(0.99)
        if np.isnan(ylimit_MA):
            ylimit_MA=df['log2(baseMean)'].quantile(0.99)
    # set ylimit_volcano if none to  to 99 percentile of ['-log10(padj)']
    if not ylimit_volcano:
        #ylimit_volcano = df[df['padj']<0.05]['-log10(padj)'].quantile(0.99)
        ylimit_volcano = df[(df['padj']<0.05)&(df['log2FoldChange'].abs()>log2FoldChange_threshold)]['-log10(padj)'].quantile(0.99)
        if np.isnan(ylimit_volcano):
            ylimit_volcano=df['-log10(padj)'].quantile(0.99)
    # set xlimit if none to 99 percentile of abs(x)
    if not xlimit:
        #xlimit = df[df['log2FoldChange'].abs()>0.5]['log2FoldChange'].abs().quantile(0.99)
        #xlimit = df[(df['log2FoldChange'].abs()>0.5)]['log2FoldChange'].abs().quantile(0.95)
        xlimit = df[(df['padj']<0.05)&(df['log2FoldChange'].abs()>log2FoldChange_threshold)]['log2FoldChange'].abs().quantile(0.99)
        # if xlimit is nan set to quantile(0.99)
        if np.isnan(xlimit):
            xlimit=df['log2FoldChange'].abs().quantile(0.99)

    # add 'Marker_volcano' column for out of range data points with '-log10(padj)' > ylimit or abs('log2FoldChange') > xlimit value of 'In_Range' or 'Out_of_Range'
    df['Marker_pvalue'] = 'In_Range'
    df.loc[df['-log10(padj)']>=ylimit_volcano,'Marker_pvalue'] = 'Out_of_Range'
    #  abs('log2FoldChange') 
    df.loc[abs(df['log2FoldChange'])>=xlimit,'Marker_pvalue'] = 'Out_of_Range'

    # Ensure the required categories are present
    df['Marker_pvalue'] = df['Marker_pvalue'].astype('category')
    required_range_cats = {'In_Range', 'Out_of_Range'}
    if not required_range_cats.issubset(df['Marker_pvalue'].cat.categories):
        df['Marker_pvalue'] = df['Marker_pvalue'].cat.set_categories(['In_Range', 'Out_of_Range'])

    # add 'Marker_MA' column for out of range data points with 'log2(baseMean)' > ylimit or abs('log2FoldChange') > xlimit value of 'In_Range' or 'Out_of_Range'
    df['Marker_Counts'] = 'In_Range'
    df.loc[df['log2(baseMean)']>=ylimit_MA,'Marker_Counts'] = 'Out_of_Range'
    #  abs('log2FoldChange') 
    df.loc[abs(df['log2FoldChange'])>=xlimit,'Marker_Counts'] = 'Out_of_Range'
    # order the categories
    # Ensure the required categories are present
    df['Marker_Counts'] = df['Marker_Counts'].astype('category')
    required_range_cats = {'In_Range', 'Out_of_Range'}
    if not required_range_cats.issubset(df['Marker_Counts'].cat.categories):
        df['Marker_Counts'] = df['Marker_Counts'].cat.set_categories(['In_Range', 'Out_of_Range'])



    # replace values in the 'log2(baseMean)' column that above the ylimit_MA with the ylimit_volcano
    if ylimit_MA:
        df['log2(baseMean)'] = df['log2(baseMean)'].apply(lambda x: (ylimit_MA*0.99) if x>=ylimit_MA else x)
    else:
        ylimit_MA = df['log2(baseMean)'].max()

    # replace values in the -log10(padj) column that above the ylimit_volcano with the ylimit_volcano
    if ylimit_volcano:
        df['-log10(padj)'] = df['-log10(padj)'].apply(lambda x: (ylimit_volcano*0.99) if x>=ylimit_volcano else x)
    else:
        ylimit_volcano = df['-log10(padj)'].max()
    # replace values in the log2FoldChange column that above or below the xlimit with the xlimit
    if xlimit:
        df['log2FoldChange'] = df['log2FoldChange'].apply(lambda x: (xlimit*0.99) if x>=xlimit else x)
        df['log2FoldChange'] = df['log2FoldChange'].apply(lambda x: (-xlimit*0.99)  if x<=-xlimit else x)
    else:
        xlimit = max(abs(df['log2FoldChange'].min()), df['log2FoldChange'].max())

    print(f'ylimit_MA={ylimit_MA}, ylimit_volcano={ylimit_volcano}, xlimit={xlimit}')
     
    ### set the marker size relative to number of dots
    rel_size=df.shape[0]/333

    fig, axes = plt.subplots(2, 1,sharex=sharex,sharey=sharey,figsize=figsize)
    # sup title
    fig.suptitle(f'{title_text}\n{comparison}\n', fontsize=16)

    if label_gene_name:
        p = sns.scatterplot(data=df, x='log2FoldChange', y='-log10(padj)', hue='Significance', style='Marker_pvalue', 
                            palette=custom_palette,sizes=(rel_size),  s=rel_size, 
                            ax=axes[0])
        # add vertical  lines at 1 and -1 for log2fc
        axes[0].axvline(x=log2FoldChange_threshold, color='gray', linestyle='--',label=f'log2fc>|{log2FoldChange_threshold}| ')
        axes[0].axvline(x=-log2FoldChange_threshold, color='gray', linestyle='--',)
        p.set(xlim=(-xlimit, xlimit), ylim=(0, ylimit_volcano))
        p.set_title(f'\n-log10(padj) vs log2fc Deseq2\n')
        p.set_xlabel("log2fc Deseq2 model")
        p.set_ylabel("-log10(padj) Significance Deseq2 model")
        p.legend(#title=facet_col, 
            bbox_to_anchor=(1.25, 1), 
            loc=1, 
            borderaxespad=0.)
        #label top genes by padj
        for line in range(0,50):
            axes[0].text(df.sort_values(by='padj')['log2FoldChange'].to_list()[line],
                   df.sort_values(by='padj')['-log10(padj)'].to_list()[line],
                   df.sort_values(by='padj').gene_name.to_list()[line],
                      horizontalalignment='left', size='small', color='black')
    else:
        p = sns.scatterplot(data=df, x='log2FoldChange', y='-log10(padj)',hue='Significance', style='Marker_pvalue', 
                             palette=custom_palette, s=rel_size,  
                            ax=axes[0])
        # add vertical  lines at 1 and -1 for log2fc
        axes[0].axvline(x=log2FoldChange_threshold, color='gray', linestyle='--',label=f'log2fc>|{log2FoldChange_threshold}| ')
        axes[0].axvline(x=-log2FoldChange_threshold, color='gray', linestyle='--',)
        p.set(xlim=(-xlimit, xlimit), ylim=(0, ylimit_volcano))
        p.set_title(f'\n-log10(padj) vs log2fc Deseq2\n')
        p.set_xlabel("log2fc Deseq2 model")
        p.set_ylabel("-log10(padj) Significance Deseq2 model")
        p.legend(#title=facet_col,
            bbox_to_anchor=(1.25, 1), 
            loc=1, 
            borderaxespad=0.)
        
    if label_gene_name:
        ma = sns.scatterplot(data=df, x='log2FoldChange', y='log2(baseMean)', hue='Significance', style='Marker_Counts', 
                            palette=custom_palette,sizes=(rel_size),  s=rel_size, 
                            ax=axes[1])
        # add verticle line at y=xpercentile
        percentile_line_BM=df['log2(baseMean)'].quantile([0.7,]).values[0]
        percentile_line_BM_counts=df['baseMean'].quantile([0.7,]).values[0]
        axes[1].axhline(y=percentile_line_BM, color='black', linestyle=':',
                   label=f'70th percentile \n  counts= {int(percentile_line_BM_counts)} ')
        # add vertical  lines at 1 and -1 for log2fc
        axes[1].axvline(x=log2FoldChange_threshold, color='gray', linestyle='--',label=f'log2fc>|{log2FoldChange_threshold}| ')
        axes[1].axvline(x=-log2FoldChange_threshold, color='gray', linestyle='--',)
        ma.set(xlim=(-xlimit, xlimit), ylim=(0, ylimit_MA))
        ma.set_title(f'Counts vs log2fc Deseq2 \n')
        ma.set_ylabel("log2(baseMean) Counts Deseq2")
        ma.set_xlabel("log2fc Deseq2 model")
        ma.legend(#title=facet_col, 
            bbox_to_anchor=(1.25, 1), 
            loc=1, 
            borderaxespad=0.5)
        #label top genes by padj
        for line in range(0,50):
            axes[1].text(df.sort_values(by='padj')['log2FoldChange'].to_list()[line],
                    df.sort_values(by='padj')['log2(baseMean)'].to_list()[line],
                   df.sort_values(by='padj').gene_name.to_list()[line],
                      horizontalalignment='left', size='small', color='black')
    else:
        ma = sns.scatterplot(data=df, x='log2FoldChange', y='log2(baseMean)',hue='Significance', style='Marker_Counts', 
                             palette=custom_palette, s=rel_size,  
                            ax=axes[1])
        # add verticle line at y=xpercentile
        percentile_line_BM=df['log2(baseMean)'].quantile([0.7,]).values[0]
        percentile_line_BM_counts=df['baseMean'].quantile([0.7,]).values[0]
        axes[1].axhline(y=percentile_line_BM, color='black', linestyle=':',
                   label=f'70th percentile \n  counts= {int(percentile_line_BM_counts)} ')
        # add vertical  lines at 1 and -1 for log2fc
        axes[1].axvline(x=log2FoldChange_threshold, color='gray', linestyle='--',label=f'log2fc>|{log2FoldChange_threshold}| ')
        axes[1].axvline(x=-log2FoldChange_threshold, color='gray', linestyle='--',)
        ma.set(xlim=(-xlimit, xlimit), ylim=(0, ylimit_MA))
        ma.set_title(f'Counts vs log2fc Deseq2 \n')
        ma.set_ylabel("log2(baseMean) Counts Deseq2")
        ma.set_xlabel("log2fc Deseq2 model")
        ma.legend(#title=facet_col,
            bbox_to_anchor=(1.25, 1), 
            loc=1, 
            borderaxespad=0.)

    return fig
