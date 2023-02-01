
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