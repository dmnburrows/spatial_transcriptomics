import numpy as np
from collections import defaultdict

#==================================================================================================
def sort_data(data_path, coord_path, genes_path, scale_path, meta_filt, min_cell_per_gene=0):
#==================================================================================================

    """
    This function loads plaque data, gene expression data and cluster labels.

    Inputs: 
        data_path (str): path to plaque data
        coord_path (str): path to spot coordinates
        genes_path (str): path to gene expression data
        scale_path (str): path to scale data
        meta_filt (dataframe): dataframe or ordered seurat processed spots
        min_cell_per_gene (int): minimum number of cells per gene

    Returns:
        spot_df (dataframe): dataframe of spot metadata
        gene_df (dataframe): dataframe of gene expression data


    """
    import scanpy
    import json
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
    d.update({'plaque': plq_bool.astype(int), 'cluster': meta_filt['cluster']})
    spot_df = pd.DataFrame(data=d)

    #Put gene data into dict
    #preprocess genes
    scanpy.pp.filter_genes(gene_, min_cells=min_cell_per_gene)
    scanpy.pp.normalize_total(gene_, target_sum=1e6)
    cell_n = gene_.obs.index[keepcell_ind]
    gene_m = gene_.X.toarray()[keepcell_ind,:]
    gene_n = gene_.var_names
    e={}
    for n,k in enumerate(gene_n): e.update({k: gene_m[:,n]})
    gene_df = pd.DataFrame(data=e, index=cell_n)

    #Remove transgenes 
    gene_df = gene_df.drop(['Thy1', 'humanAPP'],axis=1)

    #Convert to micron space
    with open(scale_path, 'r') as j:
        scale_dict = json.loads(j.read())
    micron_per_pix = 65/scale_dict['spot_diameter_fullres']
    micron_x = micron_per_pix*spot_df['pxl_col_in_fullres']
    micron_y = micron_per_pix*spot_df['pxl_row_in_fullres']
    spot_df['micron_x'] = micron_x
    spot_df['micron_y'] = micron_y
    spot_df = spot_df[['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres', 'micron_x', 'micron_y', 'plaque', 'cluster']]

    #sanity check
    assert sum(spot_df['in_tissue'] == 0) == 0
    assert sum([spot_df['barcode'].values[i] != UMI[i] for i in range(len(UMI))]) == 0
    assert sum([cell_n[i] !=UMI[i]for i in range(len(UMI))])==0
    assert sum([gene_df.index[i] !=UMI[i]for i in range(len(UMI))])==0

    return(spot_df, gene_df)

#==================================================================================================
def dist_nearest_plaq(df):
#==================================================================================================
    """
    This function calculates the distance of each spot to the nearest plaque.
    Inputs:
        df: dataframe with spot coordinates and plaque coordinates
    Outputs:
        df: dataframe with spot coordinates, plaque coordinates, and distance to nearest plaque
    
    """

    #Find cortex plaque coords
    plq_pos = df[df['plaque']==1]['micron_x'].values,df[df['plaque']==1]['micron_y'].values
    plq_pos = np.asarray(plq_pos).T

    #Loop through each spot and find the nearest plaque
    #Convert into um distances!!!!
    dist_v = []
    for c in range(len(df)):
        spot = df.iloc[c]['micron_x'], df.iloc[c]['micron_y']
        dist = np.sqrt(((spot[0]-plq_pos[:,0])**2 +(spot[1]-plq_pos[:,1])**2).astype(float))
        dist_v = np.append(dist_v,np.min(dist))
        dist_v = np.asarray(dist_v)
    df['dist_nearest_plaq'] = dist_v
    return(df)

#==================================================================================================
def mean_bin_plqdist(nbins, curr_gene, curr_lab):
#==================================================================================================
    """
    This function takes as input a set of genes, and 2 dataframes of gene expression and distance to nearest plaque for a spot, and outputs
    the mean gene expression as a function of distance for over a set of bins. 

    Inputs:
        nbins: number of bins to divide the data into
        curr_gene: dataframe of gene expression for a set of genes
        curr_lab: dataframe of metadata for a set of genes
    
    Outputs:
        bin_means: mean gene expression for each bin
        bin_std: standard deviation of gene expression for each bin
        bin_edges: edges of each bin

    """
    from scipy import stats

    dist_flat = curr_lab['dist_nearest_plaq'].values
    umi_flat = curr_lab.index.values

    #sort by umis by distance
    sort_dist, sort_umi = adm.sort_2list(dist_flat, umi_flat)
    #Define number of bins
    bins = (np.linspace(np.min(sort_dist), np.max(sort_dist), nbins+1)).astype(int)
    bin_means, bin_edges, binnumber = stats.binned_statistic(sort_dist, curr_gene[sort_umi], statistic='mean', bins=bins)
    bin_std, bin_edges, binnumber = stats.binned_statistic(sort_dist, curr_gene[sort_umi], statistic='std', bins=bins)
    return(bin_means, bin_std, bin_edges)

    
    
#==================================================================================================
class Graph():
#==================================================================================================
    """
    This class builds a symmetric graph object with nodes, edges and edge weights.
    
    """
    
    def __init__(self):
        """
        self.edges is a dict of all possible next nodes
        e.g. {'X': ['A', 'B', 'C', 'E'], ...}
        self.weights has all the weights between two nodes,
        with the two nodes as a tuple as the key
        e.g. {('X', 'A'): 7, ('X', 'B'): 2, ...}
        """
        self.edges = defaultdict(list)
        self.weights = {}

    
    def add_edge(self, from_node, to_node, weight):
        # Note: assumes edges are bi-directional
        self.edges[from_node].append(to_node)
        self.edges[to_node].append(from_node)
        self.weights[(from_node, to_node)] = weight
        self.weights[(to_node, from_node)] = weight
        
        
        
#==================================================================================================
def dijsktra(graph, initial, end):
#==================================================================================================
    """
    This function uses the dijsktra method to find the shortest path lengths between two points.
    
    Input:
        graph (Graph object): Graph object of Graph class 
        initial (str): name of starting node
        end (str) : name of end node
        
    Output:
        path (list): list of nodes 
    """

    # shortest paths is a dict of nodes
    # whose value is a tuple of (previous node, weight)
    shortest_paths = {initial: (None, 0)}
    current_node = initial
    visited = set()
    
    while current_node != end:
        visited.add(current_node)
        destinations = graph.edges[current_node]
        weight_to_current_node = shortest_paths[current_node][1]

        for next_node in destinations:
            weight = graph.weights[(current_node, next_node)] + weight_to_current_node
            if next_node not in shortest_paths:
                shortest_paths[next_node] = (current_node, weight)
            else:
                current_shortest_weight = shortest_paths[next_node][1]
                if current_shortest_weight > weight:
                    shortest_paths[next_node] = (current_node, weight)
        
        next_destinations = {node: shortest_paths[node] for node in shortest_paths if node not in visited}
        if not next_destinations:
            return "Route Not Possible"
        # next node is the destination with the lowest weight
        current_node = min(next_destinations, key=lambda k: next_destinations[k][1])
    
    # Work back through destinations in shortest path
    path = []
    while current_node is not None:
        path.append(current_node)
        next_node = shortest_paths[current_node][0]
        current_node = next_node
    # Reverse path
    path = path[::-1]
    return path


#=======================
def extract_p(de_genes, data_type):
#=======================
    """
    This function extracts p values from mast de_genes matrix for continuous or discrete components

    Inputs:
        de_genes (dataframe): mast dictionary of de_genes
        data_type (str): C for continuous or D for discrete
    Returns:
        p_mat (dataframe): dataframe with p values for each gene 
    """

    p_mat = de_genes[de_genes['component'] == data_type]
    p_mat = p_mat[p_mat['contrast'] == 'adj_plq']
    p_mat['Pr(>Chisq)'] = p_mat['Pr(>Chisq)'].fillna(1)
    p_mat['p_val'] = p_mat['Pr(>Chisq)']
    return(p_mat)

#=====================================
def filter_by_genes(n, mode, de_genes, logcpm):
#=====================================
    """
    This function filters a gene matrix by defined gene list.
    
    Inputs:
        n (int): number of de genes
        mode (str): rank by MAST coefficients for continuous (C) or discrete (D)
        de_genes (dataframe): MAST gene output
        logcpm (dataframe): gene expression matrix
    
    Outputs:

    """
    
    #Find top n DE genes
    de_dict = {'C': extract_p(de_genes, 'C'), 'D': extract_p(de_genes,'D')}
    top_de = de_dict[mode].sort_values(by='p_val', ascending=True).head(n)
    #Filter CPMmat with DE genes
    sub_gm = logcpm.loc[top_de['primerid']].T
    return(sub_gm)


#========================================
def report_metrics(true, pred, pred_prob):
#========================================

    """
    This function reports the accuracy of a classifier for each class.
    Input: true = true labels, pred = predicted labels
    """
    import sklearn.metrics as metrics

    acc = sum(true == pred)/len(true)
    print('TOTAL ACCURACY (#correct predictions/#total predictions) = ' + str(np.round(acc,3)) + ' , ' + str(int(acc*len(true))) + ' of ' + str(len(true)))
    non_plq = sum(true[np.where(true==0)] == pred[np.where(true==0)]) / sum(true==0)
    plq = sum(true[np.where(true==1)] == pred[np.where(true==1)]) / sum(true==1)
    print('ACCURACY (#correct predictions/#total true non-plaque) non-plaque = ' + str(np.round(non_plq,3)) + ' , ' + str(int(non_plq*sum(true==0))) + ' of ' + str(sum(true==0)))
    print('ACCURACY (#correct predictions/#total true plaque) plaque = ' + str(np.round(plq,2)) + ' , ' + str(int(plq*sum(true==1))) + ' of ' + str(sum(true==1)))

    prec_no, prec_plq = metrics.precision_score(true, pred, average=None)
    print('PRECISION (TP/TP+FP) non-plaque = ' + str(np.round(prec_no,3)))
    print('PRECISION (TP/TP+FP) plaque= ' + str(np.round(prec_plq,3)))

    rec_no, rec_plq = metrics.recall_score(true, pred, average=None)
    print('RECALL (TP/TP+FN) Non-plaque = ' + str(np.round(rec_no,3)))
    print('RECALL (TP/TP+FN) plaque = ' + str(np.round(rec_plq,3)))

    f1_no, f1_plq = metrics.f1_score(true, pred, average=None)
    print('F1 SCORE (2*(Pre/Pre+Rec)) non-plaque = ' + str(np.round(f1_no,3)))
    print('F1 SCORE (2*(Pre/Pre+Rec)) plaque = ' + str(np.round(f1_plq,3)))

    #true positive rate
    print('TRUE POSITIVE RATE (TP/TP+FN) non-plaque = ' + str(np.round(rec_no,3)))
    print('TRUE POSITIVE RATE (TP/TP+FN) plaque = ' + str(np.round(rec_plq,3)))

    #false positive rate
    nplq_fpr = sum(true[np.where(pred == 0)]) / (sum(true[np.where(pred == 0)]) + sum(true[np.where(pred == 1)]))
    print('FALSE POSITIVE RATE (FP/FP+TN) non-plaque  = ' + str(np.round(nplq_fpr,3)))
    plq_fpr = sum(true[np.where(pred == 1)] !=1) / (sum(true[np.where(pred == 1)] !=1) + sum(true[np.where(pred==0)] == 0))
    print('FALSE POSITIVE RATE (FP/FP+TN) plaque = ' + str(np.round(plq_fpr,3)))

    #roc auc noplq
    fpr, tpr, thresholds = metrics.roc_curve(true, pred_prob[:,0], pos_label=0)
    roc_auc = metrics.auc(fpr, tpr)
    print('ROC AUC non-plaque = ' + str(np.round(roc_auc,3)))

    #roc auc plq
    fpr, tpr, thresholds = metrics.roc_curve(true, pred_prob[:,1], pos_label=1)
    roc_auc = metrics.auc(fpr, tpr)
    print('ROC AUC plaque = ' + str(np.round(roc_auc,3)))

    
#======================================
def model_cvs(model, X, y, cv):
#======================================
    """
    This function computes model cross validation scores across range of metrics.
    
    Inputs:
        model (object): sklearn model object to fit
        X (np array): training data
        y (np array): training data labels
        cv (int): number of folds
    
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import train_test_split
    import sklearn.metrics
    from sklearn.model_selection import cross_val_score
    
    #NB default pos_label is 1
    scores = cross_val_score(model, X, y, cv=cv, scoring='accuracy')
    print("Accuracy (#correct predictions/#total predictions): %0.3f (+/- %0.3f)" % (scores.mean(), scores.std()))
    print(scores)
    
    scores = cross_val_score(model, X, y, cv=cv, scoring='precision')
    print("Plaque precision (TP/TP+FP): %0.3f (+/- %0.3f)" % (scores.mean(), scores.std()))
    print(scores)
    
    scores = cross_val_score(model, X, y, cv=cv, scoring='recall')
    print("Plaque recall (TP/TP+FN): %0.3f (+/- %0.3f)" % (scores.mean(), scores.std()))
    print(scores)
    
    scores = cross_val_score(model, X, y, cv=cv, scoring='f1')
    print("Plaque f1 (2*(Pre/Pre+Rec)): %0.3f (+/- %0.3f)" % (scores.mean(), scores.std()))
    print(scores)


    scores = cross_val_score(model, X, y, cv=cv, scoring='roc_auc')
    print("ROC AUC: %0.3f (+/- %0.3f)" % (scores.mean(), scores.std()))
    print(scores)

    
    
          


