def sort_data(data_path, coord_path, genes_path, meta_filt, min_cell_per_gene=0, norm_factor=1e6):
    """
    This function loads plaque data, gene expression data and cluster labels.

    Inputs: 
        data_path (str): path to plaque data
        coord_path (str): path to spot coordinates
        genes_path (str): path to gene expression data
        meta_filt (dataframe): dataframe or ordered seurat processed spots
        min_cell_per_gene (int): minimum number of cells per gene
        norm_factor (int): normalization factor

    Returns:
        spot_df (dataframe): dataframe of spot metadata
        gene_df (dataframe): dataframe of gene expression data


    """
    import scanpy
    import json
    import numpy as np
    import pandas as pd

    # load seurat data
    UMI = meta_filt['cell'].values 


    # Load plq data
    with open(data_path, 'r') as j:
        plaq_df = json.loads(j.read())
    #extract row and column info
    row_l = np.asarray([plaq_df['oligo'][i]['row'] for i in range(len(plaq_df['oligo']))])
    col_l = np.asarray([plaq_df['oligo'][i]['col'] for i in range(len(plaq_df['oligo']))])

    #Load coord data
    coord_df = pd.read_csv(coord_path)

    #Load gene data
    gene_ = scanpy.read_10x_h5(genes_path)
    gene_.var_names_make_unique()

    #loop over each seurat ordered UMI
    plq_bool = []
    coord_l = []
    keepcell_ind = []

    for x,u in enumerate(UMI):
        #Find current position of UMI in coord file
        curr_coord = coord_df[coord_df['barcode']==u]
        ind = np.intersect1d(np.where(row_l == curr_coord['array_row'].values[0]),np.where(col_l==curr_coord['array_col'].values[0]))[0]
        assert (row_l[ind], col_l[ind]) == (curr_coord['array_row'].values[0], curr_coord['array_col'].values[0])

        #Filter coords by UMI
        if x == 0: coord_l = curr_coord
        else: coord_l = np.vstack((coord_l, curr_coord))
        #Label plaques
        if 'tissue' in plaq_df['oligo'][ind].keys(): plq_bool = np.append(plq_bool, 1)
        else: plq_bool = np.append(plq_bool, 0)

        #Filter gene file by UMIs
        keepcell_ind = (np.append(keepcell_ind, int(np.where(gene_.obs.index == u)[0]))).astype(int)
        
    #Put spot metadata into dict
    d={}
    for n,k in enumerate(coord_df.keys()): d.update({k: coord_l[:,n]})
    d.update({'plaque': plq_bool.astype(int), 'cluster': meta_filt['cluster'], 'anot': meta_filt['anot']})
    spot_df = pd.DataFrame(data=d)

    #Put gene data into dict
    #preprocess genes
    #scanpy.pp.filter_genes(gene_, min_cells=min_cell_per_gene)
    scanpy.pp.normalize_total(gene_, target_sum=norm_factor)
    cell_n = gene_.obs.index[keepcell_ind]
    gene_m = gene_.X.toarray()[keepcell_ind,:]
    gene_n = gene_.var_names
    e={}
    for n,k in enumerate(gene_n): e.update({k: gene_m[:,n]})
    gene_df = pd.DataFrame(data=e, index=cell_n)

    #Remove transgenes 
    gene_df = gene_df.drop(['Thy1', 'humanAPP'],axis=1)

    #sanity check
    assert sum(spot_df['in_tissue'] == 0) == 0
    assert sum([spot_df['barcode'].values[i] != UMI[i] for i in range(len(UMI))]) == 0
    assert sum([cell_n[i] !=UMI[i]for i in range(len(UMI))])==0
    assert sum([gene_df.index[i] !=UMI[i]for i in range(len(UMI))])==0

    return(spot_df, gene_df)



#========================================
def report_class_acc(true, pred):
#========================================

    """
    This function reports the accuracy of a classifier for each class.
    Input: true = true labels, pred = predicted labels
    """

    non_plq = sum(true[np.where(true==0)] == pred[np.where(true==0)]) / sum(true==0)
    plq = sum(true[np.where(true==1)] == pred[np.where(true==1)]) / sum(true==1)
    print('Non-plaque accuracy = ' + str(np.round(non_plq,3)) + ' , ' + str(int(non_plq*sum(true==0))) + ' of ' + str(sum(true==0)))
    print('Plaque accuracy = ' + str(np.round(plq,2)) + ' , ' + str(int(plq*sum(true==1))) + ' of ' + str(sum(true==1)))