import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import scvelo as scv
from sklearn.neighbors import NearestNeighbors
from scipy.optimize import fsolve
from scipy.optimize import lsq_linear
from scipy.spatial.distance import cdist
from scipy.linalg import lstsq
from scipy.sparse import csr_matrix
from scipy.signal import correlate
from scipy.signal import correlation_lags
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx

def graph_fatecompass(adata, mode, basis, components, n_neighbors):
    """ Function used to compute the nearest neighbor graph in the reduced space. 
    
    Arguments
    ---------
    adata: class:`~anndata.AnnData`
        Annotated data matrix.
    mode: `str`
        which drift will be used for the transition probabilities, it can be 'potential' or 'velocity'.
    basis: `str`
        dimensionality reduction method compute the nearest neighbor graph.
    components: `int`
        number of dimensions to compute the nearest neighbor graph. 
    n_neighbors: `int`
        number of neighbors to compute the graph.
        
    Returns
    -------
    indices_fatecompass: `.obsm`
        indices of the nearest neighbors.    
    """
    n_neighbors = n_neighbors + 1
    if mode == 'potential':
        print("\nComputing dimensions in reduced space\n")
        if basis == 'pca':
            dimensions = sc.tl.pca(adata, copy= True)
            dimensions = dimensions.obsm['X_pca'][:,0:components]
        elif basis == 'umap': 
            dimensions = sc.tl.umap(adata, n_components=components, copy= True)
            dimensions = dimensions.obsm['X_umap']
        else:
            print("Please specify a 'pca' or 'umap' as basis for dimensionality reduction.")
            sys.exit()
    elif mode == 'velocity': 
        print("\nComputing dimensions in reduced space\n")
        if basis == 'pca':
            dimensions = sc.tl.pca(adata, copy= True)
            scv.tl.velocity_embedding(dimensions, basis='pca')
            velo_dimensions = dimensions.obsm['velocity_pca'][:,0:components]
            dimensions = dimensions.obsm['X_pca'][:,0:components]
            adata.obsm['dimensions_fatecompass'] = dimensions
            adata.obsm['velo_dimensions_fatecompass'] = velo_dimensions
            print("-->added\n'dimensions_fatecompass' and 'velo_dimensions_fatecompass', dimension and embedded velocities (adata.obsm)")
        elif basis == 'umap': 
            dimensions = sc.tl.umap(adata, n_components=components, copy= True)
            scv.tl.velocity_embedding(dimensions, basis='umap')
            velo_dimensions = dimensions.obsm['velocity_umap']
            dimensions = dimensions.obsm['X_umap']
            adata.obsm['dimensions_fatecompass'] = dimensions
            adata.obsm['velo_dimensions_fatecompass'] = velo_dimensions
            print("-->added\n'dimensions_fatecompass' and 'velo_dimensions_fatecompass', dimension and embedded velocities (adata.obsm)")
        else:
            print("Please specify a 'pca' or 'umap' as basis for dimensionality reduction.")
            sys.exit()

    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='ball_tree',leaf_size=50,radius=1).fit(dimensions)
    distances, indices = nbrs.kneighbors(dimensions)
    
    adata.obsm['indices_fatecompass'] = indices[:,1:]
    
    print("\nFinished -->added\n 'indices_fatecompass', knn graph (adata.obsm)\n")
    
    return  

def average_number_connections(Dv,D,cutoff):
    """ Function used in the heuristic to fit the diffusion coefficient for the transition probabilities using RNA velocity. 
    
    Arguments
    ---------
    Dv: class:`~numpy.array`
        pair-wise distance between the current state and the future state.
    D: function handle 
    cutoff: `int`
    
    Returns
    -------
    n: average number of connections   
    """
    Pv = np.exp(-(np.power(Dv,2)) / D)
    Pv = Pv / np.sum(Pv, axis=0)
    n  = np.mean(np.sum(Pv>cutoff,axis=0))
    return n

def randp(p):
    """ Function to generate random number used in the stochastic simulations. 
    
    Arguments
    ---------
    p: class:`~numpy.array`
        vector with the normalized weight of all the possible transitions for a given cell.
    
    Returns
    -------
    i: `int`
        random number. 
    """
    c = p.cumsum(axis=0)
    r = c[-1] * np.random.rand()
    i = 0
    while r > c[i]:
        i += 1        
    return i

def rna_velocity_driven_transition_probabilities(adata,cutoff,D0):
    """ Function to get transition probabilities using RNA velocity as drift. 
    
    Arguments
    ---------
    adata: class:`~anndata.AnnData`
        Annotated data matrix.
        
    Returns
    -------
    W: class:`list`
        list with number of elements equal to number of cells. Each entry is an array with the normalized weight of all the possible transitions for a given cell.
    N: class:`list`
        list with number of elements equal to number of cells. Each entry is an array with the indices of all the possible transitions for a given cell.    
    """
    
    # 1. get nearest neighbor graph distances, dimensions and embeed velocities 
    
    n_neighbors     = np.array(adata.obsm['indices_fatecompass'].shape)[1] 
    dimensions      = adata.obsm['dimensions_fatecompass'] 
    distances       = cdist(dimensions,dimensions)
    distances.sort(axis=0)
    distances       = distances[0:n_neighbors,:]
    velo_dimensions = adata.obsm['velo_dimensions_fatecompass']
    
   
    # 3. Fit delta t and diffusion coeffient for the kernel.
    
    #cutoff = 0.0095 #0.001
    #D0 = 1e-3 #0.05

    dt = np.mean(distances) / np.mean(np.linalg.norm(velo_dimensions, axis=1))
    Dv = cdist(dimensions, dimensions + velo_dimensions*dt)
    
    avercon = lambda D: average_number_connections(Dv,D,cutoff) - 2*n_neighbors
    D = fsolve(avercon, D0)
    
    print("\nFitted value of D: ", D, "\nFitted value of dt: ", dt)
    
    # 4. Getting transition probabilities.
    
    Pv = np.exp(-(np.power(Dv,2)) / D)
    Pv = Pv / np.sum(Pv, axis=0)
    P  = Pv
    
    W = list()
    N = list()

    for i in range(len(P)):
        N.append(P[:,i].nonzero())
        W.append(P[N[i],i])
        W[i] = W[i]/np.sum(W[i])
        
    print("\n*****Transition probabilities successfully computed.*****\n")
        
    return W,N

def rna_velocity_driven_stochastic_simulations(adata, root, cell_types_key, numiter=1e3, numsimcells=1e3, cutoff=0.0095,D0=1e-3):# callable STEP 1a of FateCompass
    """ Function to perform stochastic simulations using the velocity driven transition probabilities. 
    
    Arguments
    ---------
    adata: class:`~anndata.AnnData`
        Annotated data matrix.
    root: `int`
        root cell to initiate the stochastic simulations. 
    cell_types_key: `str`
        key of the observations clustering to consider.
    numiter: `int`
        number of iterations for the Monte Carlo sampling algortihm. Default: 1e3
    numsimcells: `int`
        number of trajectories to simulate. Default: 1e3
    cutoff: `int`
        parameter for transition probabilities. Default: 9.5e-3
    D0: `int`
        parameter for transition probabilities. Default: 1e-3
         
    Returns
    -------
    states: class:`~numpy.array`
        simulated differentiation trajectories.  
    num_trajectories: class:`~dictionary` 
        number of trajectories ending in specific fates. 
    """
    
    # 1. Estimate transition probabilities using RNA velocity kernel
    
    W , N = rna_velocity_driven_transition_probabilities(adata,cutoff,D0)
    
    # 2. Perform stochastic simulations using Monte Carlo algorithm.
    
    print("\n*****Performing stochastic simulations.*****\n")
    
    keys = adata.obs[cell_types_key].cat.categories
    num_trajectories = {keys[i]: 0 for i in range(len(keys))}    
    states = np.zeros([numiter,numsimcells], dtype=int)
    for f in range(numsimcells):
        if f%100 == 0:
            print("sim cell\t",f)
            
        states[0,f] = root
        for i in range(numiter-1):
            ind = states[i,f]
            n = randp(W[ind][0])
            states[i+1, f] = N[ind][0][n]
            
        end_fate = adata.obs[cell_types_key].values[states[numiter-1,f]]
        num_trajectories[end_fate] += 1
        
    adata.uns['states'] = states
    adata.uns['num_trajectories']  = num_trajectories
    
    print("\nFinished -->added\n 'states' and 'num_trajectories', stochastic trajectories (adata.uns)\n")
        
    return  

def diff_potential_driven_transition_probabilities(adata, cell_types_key, mode, initial_fate, end_fates):
    """ Function to get transition probabilities using differentiation potential as drift. 
    
    Arguments
    ---------
    adata: class:`~anndata.AnnData`
        Annotated data matrix.
    cell_types_key: `str`
        key of the observations clustering to consider.
    mode: `str`
        way for providing prior biological knowledge. It can be 'cell_types', 'marker_genes', or 'prior_knowledge_indices'
    initial_fate: `str`
        info about the initial fate, it'll depend on the selected `mode`. It can be `['name_initial_cell_type']`, `['marker_gene_initial_cell_type']`, or `[]`.  
    end_fates: `str`
        infor about the final fate, it'll depend on the selected `mode`. It can be `['name_final_cell_type_1',...]`, `['marker_gene_final_cell_type_1',...]`, or `[]`. 
        
    Returns
    -------
    W: class:`~numpy.array`
        array where each entry is the potential energy of a given cell.   
    root: `int`
        index of the initial fate
    """
    
    indices = adata.obsm['indices_fatecompass']
    
    if adata.n_obs < 5e3: 
        den = int(1/0.1) 
    elif adata.n_obs > 5e3 and adata.n_obs < 1e4:
        den = int(1/0.05)
    elif adata.n_obs > 1e4:
        den = int(1/0.02)
    
    adata.X = csr_matrix(adata.X).copy()
    if mode == 'cell_types':    
        root       = np.flatnonzero(adata.obs[cell_types_key] == initial_fate)[9:10]
        end_nodes  = [np.flatnonzero(adata.obs[cell_types_key] == i)[0:den] for i in end_fates]
        end_nodes  = np.concatenate((end_nodes))
        nodes = np.concatenate((root,end_nodes))
    elif mode == 'marker_genes':
        root       = np.flatnonzero(pd.Categorical(list(map(str,list(adata[:,[initial_fate]].X > np.max(adata[:,[initial_fate]].X.todense()-0.5))))))[0]
        end_nodes  = [np.flatnonzero(pd.Categorical(list(map(str,list(adata[:,[i]].X > np.max(adata[:,[i]].X.todense()-0.5))))))[0:den] for i in end_fates]
        end_nodes  = np.concatenate((end_nodes))
        nodes = np.zeros([len(end_nodes)+1,])
        nodes[0] = root
        nodes[1:] = end_nodes
    elif mode == 'prior_knowledge_indices':
        nodes = adata.uns['nodes_fatecompass']#IMPROVE!!!
        root  = nodes[0]
    
    node_type = 1/den * np.ones([1,len(nodes)])
    node_type[0][0] = -1
    
    M = np.zeros(shape=(len(indices),len(indices)))
    for i in range(len(indices)):
        M[i,indices[i,:]] = 1
        
    G = nx.from_numpy_matrix(M)
    
    W = np.zeros([1,len(indices)])
    for i in range(len(nodes)):
        D = nx.shortest_path_length(G, source=nodes[i])
        d = np.array([D[i] for i in range(len(indices))])
        #D = np.array(list(D.values()))
        W = W + node_type[0][i]*100 / np.power((d + 1),0.5)
        
    W = W - np.max(W)
    
    print("\n*****Differentiation potential successfully computed.*****\n")
    
    return W, root

def diff_potential_driven_stochastic_simulations(adata, cell_types_key, mode, initial_fate, end_fates, numiter=1e3, numsimcells=1e3):# callable STEP 1b of FateCompass
    """ Function to perform stochastic simulations using the differentiation potential driven transition probabilities. 
    
    Arguments
    ---------
    adata: class:`~anndata.AnnData`
        Annotated data matrix.
    cell_types_key: `str`
        key of the observations clustering to consider.
    mode: `str`
        way for providing prior biological knowledge. It can be 'cell_types', 'marker_genes', or 'prior_knowledge_indices'
    initial_fate: `str`
        info about the initial fate, it'll depend on the selected `mode`. It can be `['name_initial_cell_type']`, `['marker_gene_initial_cell_type']`, or `[]`.  
    end_fates: `str`
        infor about the final fate, it'll depend on the selected `mode`. It can be `['name_final_cell_type_1',...]`, `['marker_gene_final_cell_type_1',...]`, or `[]`. 
    numiter: `int`
        number of iterations for the Monte Carlo sampling algortihm. Default: 1e3
    numsimcells: `int`
        number of trajectories to simulate. Default: 1e3
         
    Returns
    -------
    states: class:`~numpy.array`
        simulated differentiation trajectories.    
    num_trajectories: class:`~dictionary` 
        number of trajectories ending in specific fates. 
    """
   
    W, root = diff_potential_driven_transition_probabilities(adata, cell_types_key, mode, initial_fate, end_fates)
    
    indices = adata.obsm['indices_fatecompass']
    
    n_neighbors = np.array(indices.shape)[1]

    pot = np.exp(W)
    pot = pot / np.sum(pot)
    
    adata.obs['potential'] = np.log(pot.T)
    
    print("\n'potential', differentiation potential gradient (adata.obs)\n")
    
    print("\n*****Performing stochastic simulations.*****\n")
    
    keys = adata.obs[cell_types_key].cat.categories
    num_trajectories = {keys[i]: 0 for i in range(len(keys))}   
    states = np.zeros([numiter,numsimcells], dtype=int)
    for f in range(numsimcells):
        if f%100 == 0:
            print("sim cell\t",f)
            
        states[0,f] = root
        state = root
        energy = W[0][state]
        for i in range(numiter-1):
            new_state  = indices[state,np.random.randint(n_neighbors)]
            new_energy = W[0][new_state]
            
            if (np.exp(new_energy - energy) > np.random.rand()):
                state  = new_state
                energy = new_energy
                
            states[i+1,f] = state
            
        end_fate = adata.obs[cell_types_key].values[states[numiter-1,f]]
        num_trajectories[end_fate] += 1
                   
    adata.uns['states'] = states
    adata.uns['num_trajectories']  = num_trajectories
    
    print("\nFinished -->added\n 'states' and 'num_trajectories', stochastic trajectories (adata.uns)\n")
            
    return            

def color_scatter_sim_cell(vector,d):
    if d == 2:
        c = np.linspace(1,10,len(vector))
    else:
        c = np.zeros(len(vector),3)
        colormap = matplotlib.cm.get_cmap('viridis')
        for i in range(c):
            n = np.floor(len(colormap)*(i-1)/len(c)) + 1
            c[i,:] = colormap[n,:]
            
    return c

def plot_sim_cell(adata, basis, mycell, color):# callable
    """ Function to plot one example of the stochastic trajectories. 
    
    Arguments
    ---------
    adata: class:`~anndata.AnnData`
        Annotated data matrix.
    basis: `str`
        dimensionality reduction method compute the nearest neighbor graph.
    mycell: class:`~numpy.array`
        vector with the states of a simulated cell. It can be obtain by running: `mycell = adata.uns['states'][:,0]`
    color: `str`
        key of the observations clustering for coloring the cells.
         
    Returns
    -------
    plot    
    """
    
    plt.rcdefaults()
    color_cell = color_scatter_sim_cell(mycell,2)
    clusters = adata.obs[color].values
    keys = adata.obs[color].cat.categories
    values = adata.uns[color+'_colors']
    lut1 = {keys[i]: values[i] for i in range(len(keys))}
    col_cluster = clusters.map(lut1)
    if basis == 'umap':       
        umap2 = adata.obsm['X_umap']
        plt.scatter(umap2[:,0], umap2[:,1], c=col_cluster, alpha = 0.1, s=5)
        plt.scatter(umap2[mycell,0], umap2[mycell,1], c=color_cell)
        plt.plot(umap2[mycell,0], umap2[mycell,1],c='k')
        plt.show()
    elif basis == 'pca':
        pca2 = adata.obsm['X_pca'][:,0:2]
        plt.scatter(pca2[:,0], pca2[:,1], c=col_cluster, alpha = 0.1, s=5)
        plt.scatter(pca2[mycell,0], pca2[mycell,1], c=color_cell)
        plt.plot(pca2[mycell,0], pca2[mycell,1],c='k')
        plt.show()
    
    return

def fate_probabilities(adata, cell_types_key):
    
    [numiter, numsimcells]  = np.array(adata.uns['states'].shape)
    numcells = adata.shape[0]
    keys = adata.obs[cell_types_key].cat.categories
    fate_pbb   = {keys[i]: np.zeros(numcells,) for i in range(len(keys))}
    count_cell = np.zeros(numcells,)
    
    states = adata.uns['states']
    for c in range(numcells):
        for i in range(numsimcells):
            end_fate = adata.obs[cell_types_key].values[states[numiter-1,i]]
            if c in states[:,i]:
                fate_pbb[end_fate][c] = fate_pbb[end_fate][c] + 1
                count_cell[c] = count_cell[c] + 1
                
    for f in keys:
        fate_pbb[f] = fate_pbb[f]/(count_cell+1)
        
    adata.uns['fate_pbbs'] = fate_pbb
    
    return

def plot_fate_pbb(adata, basis, color):# callable
    """ Function to plot the distribution of fate probabilities. 
    
    Arguments
    ---------
    adata: class:`~anndata.AnnData`
        Annotated data matrix.
    basis: `str`
        dimensionality reduction method compute the nearest neighbor graph.
    color: `str`
        fate for which you wish to plot the probabilities, the key must belong to the clustering.
         
    Returns
    -------
    plot    
    """
    
    plt.rcdefaults()
    col = adata.uns['fate_pbbs'][color]
    if basis == 'umap':       
        umap2 = adata.obsm['X_umap']
        plt.scatter(umap2[:,0], umap2[:,1], c=col, alpha = 0.5, s=5)
        plt.show()
    elif basis == 'pca':
        pca2 = adata.obsm['X_pca'][:,0:2]
        plt.scatter(pca2[:,0], pca2[:,1], c=col, alpha = 0.5, s=5)
        plt.show()
    
    return

def bootstrapping(gene_indexes, E, N, diff_matrix):
    """ Fuction to get the distribution of the TF activities estimate 
    
    Arguments
    ---------
    gene_indexes: list of `int`
        List with the indexes of genes.
    E: class:`~numpy.array`
        cell- and gene- normalized gene expression matrix.
    N: class:`~numpy.array`
        site counts normalized.
    diff_matrix: sparse matrix class:`~numpy.array`
        transition matrix
        
    Returns
    -------
    activities_distribution: class:`~dictionary`
        distribution of the estimate for TF activities calculates using bootstrapping.
    """
    
    print("\nInitializing bootstrapping to build the distribution of the estimate for the TF activities.\n")
    numsamples = 100
    estimate_distribution = {i: np.zeros((np.shape(E)[1],np.shape(N)[1])) for i in range(numsamples)}
    for i in range(numsamples):
        training_idx = np.sort(np.random.choice(gene_indexes,int(len(gene_indexes)*0.8),replace=False))
        E_training   = E[:,training_idx]
        N_training   = N[training_idx,:]
        
        Astar = np.linalg.lstsq(N_training, E_training.T) # Activities = [motifs x cells]
        
        estimate_distribution[i] = diff_matrix @ Astar[0].T # [cells x factors] 
        if i%10 == 0:
            print("sample\t",i)
        
    num_factors = np.shape(N)[1]
    activities_distribution = {i: np.zeros((np.shape(E)[0],numsamples)) for i in range(num_factors)}
        
    for f in range(num_factors):
        for i in range(numsamples):
            activities_distribution[f][:,i] = np.array(estimate_distribution[i][:,f]).reshape([np.shape(E)[0],])
        
    return activities_distribution

def data_diffusion_regularization(adata, gene_indexes, E, N, tolerance):
    """ Gets the maximum likelihood of the TF activities and regularize them using data diffusion. 
    
    Arguments
    ---------
    adata: class:`~anndata.AnnData`
        Annotated data matrix.
    gene_indexes: list of `int`
        List with the indexes of genes.
    E: class:`~numpy.array`
        cell- and gene- normalized gene expression matrix.
    N: class:`~numpy.array`
        site counts normalized.
    tolerance: `int`
        minimum difference between the test and training datasets used in the cross-validation scheme.
        
    Returns
    -------
    A_regularized: class:`~numpy.array`
        Regularized TF activities 
    activities_distribution: class:`~dictionary`
        distribution of the estimate for TF activities calculates using bootstrapping.
    """
    
    # 1. Setting training and testing datasets
    
    np.random.seed(0)
    
    training_idx = np.sort(np.random.choice(gene_indexes,int(len(gene_indexes)*0.8),replace=False))
    test_idx     = np.array(sorted(set(range(min(gene_indexes), max(gene_indexes)+1)).difference(training_idx)))
    #training_idx, test_idx = gene_indexes[:int(len(gene_indexes)*0.8)], gene_indexes[int(len(gene_indexes)*0.8):]
    E_training, E_test = E[:,training_idx], E[:,test_idx]
    N_training, N_test = N[training_idx,:], N[test_idx,:]
       
    # 2. Maximum- likelihood estimates of the activities
    
    #Astar = np.linalg.lstsq(N_training, E_training.T) # Activities = [motifs x cells]
    Astar = lstsq(N_training, E_training.T) # Activities = [motifs x cells]
    
    # 3. Nearest neighbor graph for regularization using data diffusion 

    indices = adata.obsm['indices_fatecompass']

    M = np.zeros(shape=(len(indices),len(indices)))
    for i in range(len(indices)):
        M[i,indices[i,:]] = 1
    
    # We want a cell's own observed values to have the highest impact on the imputation of its own values; 
    # therefore, our transition matrix M allows for self-loops, and these are the most probable steps in the random walk.
    M = M + 10 * np.eye(len(indices)) 
    M = M / M.sum(axis=1)
    M = M.T
    
    # 4. Cross-validation scheme to fit the value of t that minimizes the mean squared error (MSE)
    
    print("\nInitializing cross-validation scheme to fit the value of t for data diffusion regularization\n")
    
    A_test = Astar[0].T # [cells x motifs]
    
    MSE_training, MSE_test = [], []
    MSE_training_aux, MSE_test_aux = np.power((E_training.T - N_training @ A_test.T),2), np.power((E_test.T - N_test @ A_test.T),2)  
    MSE_training.append(np.mean(MSE_training_aux))
    MSE_test.append(np.mean(MSE_test_aux))    
    residual = MSE_test[0] - MSE_training[0]
    print("Diffusion time: ", 0, "\t Residual: ", residual)
    
    t = 0    
    while residual >= tolerance:
        t += 1
        print("Diffusion time: ", t, "\t Residual: ", residual)
        A_test = M @ A_test
        MSE_training_aux, MSE_test_aux = np.power((E_training.T - N_training @ A_test.T),2), np.power((E_test.T - N_test @ A_test.T),2)  
        MSE_training.append(np.mean(MSE_training_aux))
        MSE_test.append(np.mean(MSE_test_aux))
        residual = MSE_test[t] - MSE_training[t]
        
    # 5. Computing regularized activites with the fitted diffusion time
    
    diff_matrix = csr_matrix(np.linalg.matrix_power(M,t))

    A_regularized = diff_matrix @ Astar[0].T # [cells x factors]      
    
    activities_distribution = bootstrapping(gene_indexes, E, N, diff_matrix) 
            
    return A_regularized, activities_distribution

def tf_activities(adata, bs, tolerance=1e-2):# callable STEP 2 FateCompass
    """ Infers TF activities using a linear model of gene regulation and data diffusion. 
    
    Arguments
    ---------
    adata: class:`~anndata.AnnData`
        Annotated data matrix.
    bs: class:`~pandas.DataFrame`
        Binding sites matrix. Data frame with rows equal genes and columns equal TFs. Index: gene names, Columns: factor names.
    tolerance: `int`
        minimum difference between the test and training datasets used in the cross-validation scheme.
        
    Returns
    -------
    tf_activities: class:`~pandas.DataFrame`
        Inferred TF activities. Stored in `adata.uns`
    tf_activities_distribution: class:`~dictionary`
        distribution of the estimate for TF activities computed using bootstrapping. Stored in `adata.uns`
    """
    
    # 1. Reduce binding sites and expression matrices to filtered genes
    
    promoter_list = list(bs.index)
    adata_genes = list(adata.var.index)
    
    genes = []
    for g in adata_genes:
        if g in promoter_list:
            genes.append(g)
 
    gene_indexes = []
    for g in genes:
        gene_indexes.append( genes.index(g) )
        
    bs = bs.loc[genes,:]
    N = csr_matrix(bs) # Binding sites matrix = [genes x motifs]
    
    # Reduce expression matrix to genes present in the biding sites matrix
    adata.X = csr_matrix(adata.X).copy()
    E = adata.X[:,np.isin(adata.var_names,bs.index)]#.todense() # Expression matrix = [cells x genes]
       
    # 2. Getting and printing dimensions of the problem
    
    num_cells, num_genes = E.shape
    num_factors = np.shape(N)[1]
    
    print("Number of cells:",num_cells)
    print("Number of genes:",num_genes)
    print("Number of motifs:",num_factors)
    
    # 3. Normalization of binding sites and expression matrices
    
    # a. Site counts are normalized to sum to zero across genes     
    #N = N.to_numpy()
    Ntilde = N - N.mean(axis=0)
    print("Normalized binding site matrix [genes x motifs]: " + str(Ntilde.shape))
    
    # b. Each cell is normalized by substracting the mean expression 
    Eprime = (E.T - E.mean(axis=1).T).T
    # c. Each gene is normalized by substracting the mean expression 
    Etilde = Eprime - Eprime.mean(axis=0)
    print("Cell- and Gene- Normalized expression matrix [cells x genes]: " + str(Etilde.shape))
    
    # 4. fitting activities using data diffusion regularization 
    A_regularized, activities_distribution = data_diffusion_regularization(adata, gene_indexes, Etilde, Ntilde, tolerance) # [cells x motifs]
    
    for i in range(len(bs.columns)):
        activities_distribution[bs.columns[i]] = activities_distribution.pop(i)
    
    adata.uns['tf_activities'] = pd.DataFrame(A_regularized.T,index=bs.columns,columns=adata.obs.index).transpose()
    
    adata.uns['tf_activities_distribution'] = activities_distribution
    
    print("\nFinished -->added\n 'tf_activities' and 'tf_activities_distribution', TF activities regularized (adata.uns)\n")
    
    return 

def avg_profiles_over_trajectories(adata, cell_types_key):
    """ Function to calculate the average profiles over stochastic trajectories. 
    
    Arguments
    ---------
    adata: class:`~anndata.AnnData`
        Annotated data matrix.
    cell_types_key: `str`
        key of the observations clustering to consider.
        
    Returns
    -------
    mean_E: class:`~dictionary`
        average gene expression profile over stochastic trajectories. Stored in `adata.uns`  
    sem_E: class:`~dictionary`
        average standard error of the mean for gene expression profile over stochastic trajectories. Stored in `adata.uns` 
    mean_A: class:`~dictionary`
        average TF activity profile over stochastic trajectories. Stored in `adata.uns`  
    sem_A: class:`~dictionary`
        average standard error of the mean for TF activity profile over stochastic trajectories. Stored in `adata.uns` 
    """
    
    [numiter, numsimcells] = np.array(adata.uns['states'].shape)
    
    keys = adata.obs[cell_types_key].cat.categories
    adata.X = csr_matrix(adata.X).copy()
    E = adata.X
    A = csr_matrix(adata.uns['tf_activities'])
    
    mean_E  = {keys[i]: csr_matrix((numiter,np.size(E,axis=1))) for i in range(len(keys))}
    mean_E2 = {keys[i]: csr_matrix((numiter,np.size(E,axis=1))) for i in range(len(keys))}
    var_E   = {keys[i]: csr_matrix((numiter,np.size(E,axis=1))) for i in range(len(keys))}
    std_E   = {keys[i]: csr_matrix((numiter,np.size(E,axis=1))) for i in range(len(keys))}
    sem_E   = {keys[i]: csr_matrix((numiter,np.size(E,axis=1))) for i in range(len(keys))}

    mean_A  = {keys[i]: csr_matrix((numiter,np.size(A,axis=1))) for i in range(len(keys))}
    mean_A2 = {keys[i]: csr_matrix((numiter,np.size(A,axis=1))) for i in range(len(keys))}
    var_A   = {keys[i]: csr_matrix((numiter,np.size(A,axis=1))) for i in range(len(keys))}
    std_A   = {keys[i]: csr_matrix((numiter,np.size(A,axis=1))) for i in range(len(keys))}
    sem_A   = {keys[i]: csr_matrix((numiter,np.size(A,axis=1))) for i in range(len(keys))}

    for c in range(numsimcells):
        if c%100 == 0:
            print("sim cell\t",c)

        end_fate = adata.obs[cell_types_key].values[adata.uns['states'][numiter-1,c]]

        # Average gene expression profiles over stochastic trajectories    
        mean_E[end_fate]  = mean_E[end_fate] + E[adata.uns['states'][:,c],:]
        mean_E2[end_fate] = mean_E2[end_fate] + E[adata.uns['states'][:,c],:].power(2)
        # Average TF acttivity profiles over stochastic trajectories    
        mean_A[end_fate]  = mean_A[end_fate] + A[adata.uns['states'][:,c],:]
        mean_A2[end_fate] = mean_A2[end_fate] + A[adata.uns['states'][:,c],:].power(2)
        
    num_trajectories = adata.uns['num_trajectories']
    for i in keys:

        # Average gene expression profiles over stochastic trajectories
        mean_E[i]  = mean_E[i]/(num_trajectories[i]+1)
        mean_E2[i] = mean_E2[i]/(num_trajectories[i]+1)
        var_E[i]   = mean_E2[i] - mean_E[i].power(2)
        std_E[i]   = var_E[i].sqrt()
        sem_E[i]   = std_E[i]/np.sqrt((num_trajectories[i]+1))
        # Average TF activity profiles over stochastic trajectories
        mean_A[i]  = mean_A[i]/(num_trajectories[i]+1)
        mean_A2[i] = mean_A2[i]/(num_trajectories[i]+1)
        var_A[i]   = mean_A2[i] - mean_A[i].power(2)
        std_A[i]   = var_A[i].sqrt()
        sem_A[i]   = std_A[i]/np.sqrt((num_trajectories[i]+1))
        
    adata.uns['mean_E'] = mean_E
    adata.uns['sem_E']  = sem_E
    adata.uns['mean_A'] = mean_A
    adata.uns['sem_A']  = sem_A
    
    
    print("\nFinished --> added\n 'mean_E' and 'sem_E', average gene expression profiles over trajectories (adata.uns)\n 'mean_A' and 'sem_A', average TF activity profiles over trajectories (adata.uns)\n")
    
    return 

def plot_trajectory(adata, mode, variable, cell_types_key, trajectory):
    
    plt.rcdefaults()
    scv.settings.set_figure_params("scvelo",fontsize=16)
    
    map_tf_id = list(adata.uns['tf_activities_distribution'].keys())
    map_gene_id = list(adata.var.index)
    keys   = adata.obs[cell_types_key].cat.categories
    values = adata.uns[cell_types_key+'_colors']
    col_dict = {keys[i]: values[i] for i in range(len(keys))}
    
    if mode == 'mRNA':
    
        for j in trajectory:
            x = np.linspace(0,1,(adata.uns['mean_E'][j][10:,:]).shape[0])
            col = col_dict[j]
            for i in variable:         
                y = (adata.uns['mean_E'][j][:,map_gene_id.index(i)][10:]).toarray().reshape([len(x),])
                error = (adata.uns['sem_E'][j][:,map_gene_id.index(i)][10:]).toarray().reshape([len(x),])

                plt.plot(x,y,c=col)
                plt.fill_between(x, y+error, y-error, alpha=0.3, facecolor=col) #edgecolor = col
        plt.show()
        
    elif mode == 'activity':
        
        for j in trajectory:
            x = np.linspace(0,1,(adata.uns['mean_A'][j][10:,:]).shape[0])
            col = col_dict[j]
            for i in variable:         
                y = (adata.uns['mean_A'][j][:,map_tf_id.index(i)][10:,]).toarray().reshape([len(x),])
                error = (adata.uns['sem_A'][j][:,map_tf_id.index(i)][10:,]).toarray().reshape([len(x),])

                plt.plot(x,y,c=col)
                plt.fill_between(x, y+error, y-error, alpha=0.3, facecolor=col) #edgecolor = col
        #plt.legend()#New line
        plt.show()
           
    return

def differential_tf_activity(adata,cell_types_key):
    
    cell_types = adata.obs[cell_types_key].cat.categories
    factors    = adata.uns['tf_activities_distribution'].keys()
    
    # 1. z-score  
    
    z_score = {i: 0 for i in factors}
    for f in factors:
        mu = np.zeros([len(adata),])
        sigma = np.zeros([len(adata),])

        for c in range(len(adata)):
            mu[c] = np.mean(adata.uns['tf_activities_distribution'][f][c,:])
            sigma[c] = np.std(adata.uns['tf_activities_distribution'][f][c,:])

        z_score[f] = np.sqrt(1/len(adata) * np.sum((mu/sigma)**2))
      
    print("\n z-score --> computed\n")
        
    # 2. Variability over time 
    
    std_tf_time = {i: {f: 0 for f in factors} for i in cell_types}
    for i in cell_types:
        count = 0
        for f in factors:
            std_tf_time[i][f] = np.std(adata.uns['mean_A'][i][10:,count].toarray())##CHECK
            count += 1
            
    print("variability over time --> computed\n")
            
    # 3. Dynamic correlation 
    
    cross_corr = {i: {f: 0 for f in factors} for i in cell_types}
    time_lags  = {i: {f: 0 for f in factors} for i in cell_types}
    
    for i in cell_types:    
        for f in factors:        
            tf = f.split('_')        
            if len(tf) != 1:             
                cross_corr[i][f] = {j: 0 for j in tf}
                time_lags[i][f]  = {j: 0 for j in tf}
                
    gene_list = adata.var_names.to_list()
    tf_list   = list(factors)
    
    for i in cell_types: 
        for f in factors:        
            tf = f.split('_')
            if len(tf) == 1: 
                if tf[0] in gene_list:
                    tf = tf[0]
                    x1 = adata.uns['mean_A'][i][10:,tf_list.index(tf)].toarray()
                    x2 = adata.uns['mean_E'][i][10:,gene_list.index(tf)].toarray()
                    x1 = x1 - np.mean(x1)
                    x2 = x2 - np.mean(x2)
                    cc, lag = correlate(x1, x2)/(len(x1)*x1.std()*x2.std()), correlation_lags(len(x1),len(x2))
                    cc_max, lag_max = cc[cc.argmax()], lag[lag.argmax()]
                    cross_corr[i][f] = cc
                    time_lags[i][f]  = lag
                else:
                    cross_corr[i][f] = 'NaN'
                    time_lags[i][f]  = 'NaN'
            else:
                k = 1
                while k <= len(tf):
                    if tf[k-1] in gene_list:
                        x1 = adata.uns['mean_A'][i][10:,tf_list.index(f)].toarray()
                        x2 = adata.uns['mean_E'][i][10:,gene_list.index(tf[k-1])].toarray()
                        x1 = x1 - np.mean(x1)
                        x2 = x2 - np.mean(x2)
                        cc, lag = correlate(x1, x2)/(len(x1)*x1.std()*x2.std()), correlation_lags(len(x1),len(x2))
                        cc_max, lag_max = cc[cc.argmax()], lag[lag.argmax()]
                        cross_corr[i][f][tf[k-1]] = cc
                        time_lags[i][f][tf[k-1]]  = lag
                    else:
                        cross_corr[i][f][tf[k-1]] = 'NaN'
                        time_lags[i][f][tf[k-1]]  = 'NaN'
                    k += 1
                    
    print("dynamical correlation --> computed\n")
                    
    
    adata.uns['z_score'] = z_score
    
    adata.uns['std_tf_time'] = std_tf_time
    
    adata.uns['cross_corr'] = cross_corr
    
    adata.uns['time_lags'] = time_lags
    
    return 

def ksdensity_fatecompass(adata, criterion, cell_types_key=[], trajectory=[]):
    
    keys   = adata.obs[cell_types_key].cat.categories
    values = adata.uns[cell_types_key+'_colors']
    col_dict = {keys[i]: values[i] for i in range(len(keys))}
    factors = adata.uns['z_score'].keys()
    
    fig, axs = plt.subplots(1, len(criterion), figsize=(15, 5))
    
    count = 0
    for i in criterion:
        if i == 'z_score':
            data = np.array([adata.uns['z_score'][f] for f in factors])
            kde = stats.gaussian_kde(data)
            x = np.linspace(data.min(), data.max(), 100)
            p = kde(x)
            axs[count].plot(x,p,color='k')
            axs[count].axvline(x[p.argmax()],color='r')
            axs[count].text(x[p.argmax()],0,str(x[p.argmax()]),rotation=90)
            axs[count].set_title(i)

        elif i == 'std_tf_time':
            count_t = 0
            x_t     = np.zeros((len(trajectory),))
            for t in trajectory:
                col = col_dict[t]
                data = np.array([adata.uns['std_tf_time'][t][f] for f in factors])
                kde = stats.gaussian_kde(data)
                x = np.linspace(data.min(), data.max(), 100)
                p = kde(x)
                x_t[count_t] = x[p.argmax()]
                axs[count].plot(x,p,c=col,)
                count_t += 1 
            axs[count].legend(trajectory)
            axs[count].axvline(x_t.mean(),color='r')
            axs[count].text(x_t.mean(),0,str(x_t.mean()),rotation=90)
            axs[count].set_title(i)
        count += 1
    return

def get_df_differential_tf_activity(adata, fates, thresholds):
    
    # Format data frame
    factors = adata.uns['z_score'].keys()
    col_names = ['TFs']
    for i in fates:
        col_names.append(str('std_tf_time_'+i))
    for i in fates:
        col_names.append(str('max_cross_corr_'+i))
        col_names.append(str('time_max_cross_corr_'+i))
    col_names.append('z_score')
    col_names.append('FateCompass_prediction')
    df = pd.DataFrame(columns=col_names,index=factors)
    
    # Pass values of variability over time
    for i in fates:
        df[str('std_tf_time_'+i)] = adata.uns['std_tf_time'][i].values()
        
    # Pass values of z-score
    df['z_score'] = adata.uns['z_score'].values()
    
    # Pass values of TFs - Separate each motif family by individual TFs
    count = 0
    for i in factors:
        position = np.ones(len(df), dtype=int)
        tf = i.split('_')
        if len(tf) == 1: 
            df['TFs'].loc[i] = tf[0]
        else: 
            position[count] = len(tf)
            df = df.iloc[np.arange(len(df)).repeat(position)]
            df['TFs'].loc[i] = tf   
        count = count + len(tf) 
        
    # Pass values of cross-correlation 
    for f in fates:
        for i in range(len(df['TFs'])):
            row = df.iloc[i]
            motif = row.name
            tf = row['TFs']
            if len(motif.split('_')) == 1:
                cc   = adata.uns['cross_corr'][f][tf]
                tlag = adata.uns['time_lags'][f][tf]
                if cc != 'NaN':
                    df[str('max_cross_corr_'+f)].iloc[i] = cc[cc.argmax()][0]
                    df[str('time_max_cross_corr_'+f)].iloc[i] = tlag[cc.argmax()]
                else: 
                    df[str('max_cross_corr_'+f)].iloc[i] = 'NaN'
                    df[str('time_max_cross_corr_'+f)].iloc[i] = 'NaN'
            else: 
                cc   = adata.uns['cross_corr'][f][motif][tf]
                tlag = adata.uns['time_lags'][f][motif][tf]
                if cc != 'NaN':
                    df[str('max_cross_corr_'+f)].iloc[i] = cc[cc.argmax()][0]
                    df[str('time_max_cross_corr_'+f)].iloc[i] = tlag[cc.argmax()]
                else: 
                    df[str('max_cross_corr_'+f)].iloc[i] = 'NaN'   
                    df[str('time_max_cross_corr_'+f)].iloc[i] = 'NaN'
                    
    for f in fates: 
        var_th = thresholds['variability']
        zsc_th = thresholds['z_score']
        cor_th = thresholds['correlation']
        for i in range(len(df['TFs'])):
            var_tf = float(df[str('std_tf_time_'+f)].iloc[i])
            zsc_tf = float(df['z_score'].iloc[i])
            cor_tf = float(df[str('max_cross_corr_'+f)].iloc[i])
            if var_tf >= var_th and zsc_tf >= zsc_th and cor_tf >= cor_th:
                if pd.isna(df['FateCompass_prediction'].iloc[i]) == True:
                    df['FateCompass_prediction'].iloc[i] = str(f)
                else:
                    prediction = df['FateCompass_prediction'].iloc[i]
                    df['FateCompass_prediction'].iloc[i] = str(prediction)+str(f)
    
    return df