
#====================================
def sample_rates(n_clusts, n_cells, n_genes, rate_range):
#====================================

    """
        This function samples random rates from a uniform distribution for every cell class, 
        and then samples from a Poisson distribution with rates for all cells in a given class, 
        returning a n x d array of simulated transcript counts, where n = cells, d = genes. 
        
   
    Inputs:
        n_clusts (int): number of cell types
        n_cells (int): number of cells
        n_genes (int): number of genes
        rate_range (tuple): min and max of uniform distribution for sampling the rates
    
    Outputs:
        cell_rate_mat (np array): n x d array of simulated transcript rates which inputs to Poisson, where n = cells, d = genes. 
        cell_counts (np array): n x d array of simulated transcript counts, where n = cells, d = genes. 
        clust_vec (np array): 1 x n array of cell types for each cell, where n = spots.
        """
    import numpy as np

    #Sample random rates from a uniform distribution for each class and each gene
    rate_mat = np.zeros((n_clusts, n_genes))
    rate_mat = np.random.uniform(low=rate_range[0], high=rate_range[1], size=(rate_mat.shape)).astype(int)
    cell_rate_mat = np.repeat(rate_mat, int(n_cells/n_clusts), axis=0) #get the rates for each cell - n cells in each group with same rates
    cell_rate_mat = np.vstack((cell_rate_mat, (np.repeat(np.reshape(cell_rate_mat[-1],(1,cell_rate_mat[-1].shape[0])), n_cells % n_clusts, axis=0)))) #Add any residual cells into the last cluster
    clust_vec = np.repeat(np.arange(0,n_clusts), int(n_cells/n_clusts), axis=0) #vector of cells labelled by their cluster
    clust_vec = np.concatenate((clust_vec, np.full(n_cells % n_clusts, n_clusts-1))) #add residual cells

    #Sample from Poisson for each cell given its cell type
    cell_counts = np.random.poisson(cell_rate_mat, size=cell_rate_mat.shape)

    return(cell_rate_mat, cell_counts, clust_vec)


#=============================
def proportions(inp, n_clusts):
#=============================
    """
    This function converts a vector of cell type numbers into proportions.
    
    Inputs:
        n_clusts (int): number of cell types
        inp (np array): a gx1 array of cell type numbers, where g is the number of cells in this group 
        
    returns:
        props (np array): a 1xd array of proportions of each cell type, from 0 to d
    
    """
    import numpy as np

    props = np.zeros(n_clusts) #create an array of 0s as long as the number of cell types
    unq = np.unique(inp, return_counts=True) #find unique cell types in this group and their counts
    partial_props = unq[1]/sum(unq[1]) #find the proportions of each cell type
    count=0
    for o in range(len(props)): #loop over each of the cell types and add in proportions in order
        if o in unq[0]: 
            props[o] = partial_props[count]
            count+=1
    return(props)


#====================================
def generate_spots(n_cells, n_clusts, cell_counts, clust_vec):
#====================================

    """
        This function groups cells together into random groups (spots) and calculates the total 
        expression and proportion of each cell type in that spot. 
        
   
    Inputs:
        n_cells (int): number of cells
        n_clusts (int): number of celltypes
        cell_counts (np array): k x d array of transcript counts, where k = cells, d = genes. 
        clust_vec (np array): 1 x n array of cell types for each cell, where n = spots.

    Outputs:
        n_spots (int): number of simulated spots
        spots (np array): y x d array of simulated spot counts, where y = spots, d = genes. 
        prop_vec (np array): y x k array of calculated proportions at each spot, where y = spots, k = celltypes.
        """

    import numpy as np
    import random

    #Generate groupings for each spot
    #=============================================
    orig_cells = np.arange(0,cell_counts.shape[0]) 
    random.shuffle(orig_cells) #randomly shuffle indeces

    #split into n spots
    chunk_size = []

    #loop over 100 group sizes
    for f in range(1,100): 
        chunk_size = np.append(chunk_size, np.full(int(np.random.uniform((n_cells/20), (n_cells/10))), f)) #randomly sample the number of groups of size f - range is defined as n_cells/10 so that ~10 cells per group is max
        if sum(chunk_size) >= n_cells:
            break

    chunk_size = chunk_size[np.cumsum(chunk_size) <= n_cells] #reduce groups to match number of cells
    groups = np.split(orig_cells, np.cumsum(chunk_size).astype(int)) #group by chunk size
    if sum(groups[-1]) == 0: groups=groups[:-1] #Remove empty groups at end due to chunking residuals

    if int(np.sum(chunk_size)) != n_cells:
        if int(np.sum(chunk_size) + len(groups[-1])) != n_cells:
            print('Grouping error - number of grouped cells /= number of total cells')

    n_spots = len(groups)
    #Mix cells together for spots and calculate true proportions
    #==================================================================
    spots = np.zeros((n_spots, cell_counts.shape[1])) #spots x genes
    prop_vec = np.zeros((len(groups),n_clusts)) #array of proportions for spots x celltypes
    #loop over each group
    for g in range(len(groups)):
        spots[g] = np.sum(cell_counts[groups[g]], axis=0) #generate sum of gene expression
        cell_vec = clust_vec[groups[g]] #vector of celltypes for each spot
        prop_vec[g] = proportions(cell_vec, n_clusts) #calculate the true proportions

    return(n_spots, spots, prop_vec)

#=============================
def mean_exp(n_clusts, n_genes, cell_counts, clust_vec):
#=============================
    """
    This function calculates mean expression of celltypes, given observed transcript counts and cell type identities.
    
    Inputs:
        n_clusts (int): number of celltypes
        n_genes (int): number of genes
        cell_counts (np array): k x d array of transcript counts, where k = cells, d = genes. 
        clust_vec (np array): 1 x n array of cell types for each cell, where n = spots.
        
    returns:
        mean_ref (np array): a k x d array of mean expression, where k = cells, d = genes. 
    
    """
    import numpy as np

    #Calculate mean gene expression per group
    mean_ref = np.zeros((n_clusts, n_genes))
    for i in range(n_clusts):
        mean_ref[i] = np.mean(cell_counts[clust_vec == i], axis=0)

    return(mean_ref)


#================================    
class simulate_cell_mix: 
#================================    
    """
    Class to simulate spatial transciptomic spot data, where each spot contains a mixture of cells of different gene expression patterns.
    
    """
    
    #========================
    def __init__(self, n_clusts, n_cells, n_genes):
    #========================
        self.n_clusts = n_clusts # number of cell types
        self.n_cells = n_cells # number of cells
        self.n_genes = n_genes # number of genes
        print('Loaded parameters: ' + str(n_clusts) + ' cell types , ' + str(n_cells) + ' cells, & ' + str(n_genes) + ' genes.' )


    #====================================
    def simulate_gene_exp(self, rate_range):
    #====================================
        
        """
        This functions creates spots from simulated gene expression data. 
   
        Inputs:
            rate_range (tuple): min and max of uniform distribution for sampling the rates

    
        """
        #Simulate transcript counts
        _, self.cell_counts , self.clust_vec = sample_rates(self.n_clusts, self.n_cells, self.n_genes, rate_range)

        #Randomly mix into spots
        self.n_spots, self.spots, self.prop_vec = generate_spots(self.n_cells, self.n_clusts, self.cell_counts, self.clust_vec)

        #Calculate mean gene expression from reference
        self.mean_exps = mean_exp(self.n_clusts, self.n_genes, self.cell_counts, self.clust_vec)

        print('Created spot mixtures from simulated data: ' + str(self.n_spots) + ' spots.')

        return(self)
    
    
    #====================================
    def real_gene_exp(self, cell_counts, clust_vec):
    #====================================
        
        """
        This functions creates spots from simulated gene expression data. 
   
        Inputs:
            rate_range (tuple): min and max of uniform distribution for sampling the rates
    
        """

        #Randomly mix into spots
        self.n_spots, self.spots, self.prop_vec = generate_spots(self.n_cells, self.n_clusts, cell_counts, clust_vec)

        #Calculate mean gene expression from reference
        self.mean_exps = mean_exp(self.n_clusts, self.n_genes, self.cell_counts, self.clust_vec)
        print('Created spot mixtures from real data: ' + str(self.n_spots) + ' spots.')
        return(self)
    

#==============================================================
def add_noise(spots, per, a_std, g_std, e_std):
#==============================================================

    """
    This function adds different types of noise to simulated spots data. 
    
    Inputs:
        spots (np array): y x d array of simulated spot counts, where y = spots, d = genes.
        per (int): percentage of dropout genes per spot
        a_std (int): pixel shared gaussian noise - standard deviations 
        g_std (int): gene shared gaussian noise - standard deviations
        e_std (int): independent spot and gene noise - standard deviations
        
    returns:
        spots (np array): y x d array of simulated noisy spot counts, where y = spots, d = genes.
    
    """
    import numpy as np

    #Dropout a certain percentage of genes
    if per!=None:
        for i in range(spots.shape[0]):
            rand_ind = np.random.choice(np.arange(spots.shape[1]), size = int((per/100) * spots.shape[1]), replace=False) #random index for selecting
            spots[i,rand_ind]=0
    
    if e_std!=None:
        #Add random noise and make int and remove negatives
        spots = spots+np.random.normal(0, e_std, (spots.shape)) ##EXPONENTIAL OR GAUSSIAN?

    if g_std!=None:
        #gamma - over each gene
        gamma = np.random.normal(0, g_std, (spots.shape[1])) ##EXPONENTIAL OR GAUSSIAN?
        gamma_mat = np.asarray([gamma for i in range(spots.shape[0])]) #Repeat across columns for elementwise addition
        spots = spots+gamma_mat

    if a_std!=None:
        #alpha - over each spot
        alpha = np.random.normal(0, a_std, (spots.shape[0])) ##EXPONENTIAL OR GAUSSIAN?
        alpha_mat = np.asarray([alpha for i in range(spots.shape[1])]).T #Repeat across columns for elementwise addition
        spots = spots+alpha_mat

    spots = spots.astype(int)  #REMOVE?
    spots[spots < 0] = 0
    spots +=1 #remove any zeros

    return(spots)
