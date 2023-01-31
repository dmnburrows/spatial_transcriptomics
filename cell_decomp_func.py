
#=============================
def proportions(inp, n_clusts):
#=============================
    """
    This function creates a series of empty lists of the same dimension.
    
    
    Inputs:
        n_clusts (int): number of cell types
        inp (np array): a gx1 array of cell type numbers, where g is the number of cells in this group 
        
    returns:
        out_l (list of list): list of list
    
    """
    import numpy as np

    props = np.zeros(n_clusts) #create an array of 0s as long as the number of cell types
    unq = np.unique(inp, return_counts=True) #find unique cell types in this group and theor counts
    partial_props = unq[1]/sum(unq[1]) #find the proportions of each cell type
    count=0
    for o in range(len(props)): #loop over each of the cell types and add in proportions in order
        if o in unq[0]: 
            props[o] = partial_props[count]
            count+=1
    return(props)