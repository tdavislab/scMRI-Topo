"""Helper funtions for inference with Structural Correlation Graphs."""
from operator import itemgetter
from os.path import join
import numpy as np
import scipy.sparse.csgraph as csg


def loadSCGData(ICN):
    """Load gray matter intensities for specified SCGs.

        The dataset contains 49 ASD and 49 TDC subjects.
        Number of ROIs varies by ICN:
            'Global': number of ROIs is 7266
            'DMN'   : number of ROIs is 39
            'SN'    : number of ROIs is 32
            'ECN'   : number of ROIs is 39
        The function returns two 49 x m numpy arrays containing gray matter
        densities of ROIs of the specified ICN for ASD and TDC subjects.

    Parameters
    ----------
    ICN : String
        Four possible choices are: 'Global', 'DMN', 'SN', 'ECN'

    Returns
    -------
    icn_asd: 2D numpy array
        49 x m array of floating point values.
        Contains gray matter densities for m ROIs corresponding to the
        specified ICN, for 49 ASD subjects in the dataset.
    icn_tdc: 2D numpy array
        49 x m array of floating point values.
        Contains gray matter densities for m ROIs corresponding to the
        specified ICN, for 49 TDC subjects in the dataset.

    """
    dataDir = './data'
    if(ICN == 'Global'):
        fname = join(dataDir, '7266_Raw_Densities.csv')
        data = np.loadtxt(fname, delimiter=",", skiprows=0)
        (m, n) = np.shape(data)
        icn_asd = data[:49, :]
        icn_tdc = data[49:, :]

    elif(ICN in ['DMN', 'SN', 'ECN']):
        fname = join(dataDir, 'ICN_SUB_Densities.csv')
        data = np.loadtxt(fname, delimiter=",", skiprows=1)
        (m, n) = np.shape(data)
        if(ICN == 'DMN'):
            colrange = range(2, 41)
        elif(ICN == 'SN'):
            colrange = range(41, 73)
        else:
            colrange = range(73, 112)

        icn_asd = data[np.ix_(range(0, 49), colrange)]
        icn_tdc = data[np.ix_(range(49, 98), colrange)]

    else:
        print("Error: No such ICN\n")
        print("possible choices: 'Global', 'DMN', 'SN', 'ECN'")

    return(icn_asd, icn_tdc)


def getMST(icn_asd, icn_tdc):
    """Construct Minimum Spanning Tree of SCG.

    We use (1. - abs(Corr(x,y))) as the weight for an edge between nodes x, y.

    Parameters
    ----------
    icn_asd : 2D numpy array of floats
        49 x m array of gray matter densities of m ROIs of specified ICN for
        the 49 ASD subjects in the dataset.
    icn_tdc : 2D numpy array of floats
        49 x m array of gray matter densities of m ROIs of specified ICN for
        the 49 TDC subjects in the dataset.

    Returns
    -------
    mst_asd : list of 3-tuples (int, int, float)
        sorted list of (source, target, weight) tuples corresponding to
        the edges in the minimum spanning tree of SCG for ASD subjects.
    mst_tdc : list of 3-tuples (int, int, float)
        sorted list of (source, target, weight) tuples corresponding to
        the edges in the minimum spanning tree of SCG for TDC subjects.
    """
    cor_asd = abs(np.corrcoef(icn_asd, rowvar=False))
    cor_tdc = abs(np.corrcoef(icn_tdc, rowvar=False))
    if not(np.shape(cor_tdc) == np.shape(cor_asd)):
        print('Error : correlation matrices have different shape')
        exit(1)

    mat_asd = 1.0 - cor_asd
    np.fill_diagonal(mat_asd, np.inf)
    G_asd = csg.csgraph_from_dense(mat_asd, null_value=np.inf)
    T_asd = csg.minimum_spanning_tree(G_asd)
    mst_asd = sorted(zip(T_asd.nonzero()[0], T_asd.nonzero()[1],
                         (1.0 - T_asd.data)), key=itemgetter(2))

    mat_tdc = 1.0 - cor_tdc
    np.fill_diagonal(mat_tdc, np.inf)
    G_tdc = csg.csgraph_from_dense(mat_tdc, null_value=np.inf)
    T_tdc = csg.minimum_spanning_tree(G_tdc)
    mst_tdc = sorted(zip(T_tdc.nonzero()[0], T_tdc.nonzero()[1],
                         (1.0 - T_tdc.data)), key=itemgetter(2))


    return(mst_asd, mst_tdc)


def getBettiCurve(mst_asd, mst_tdc):
    """Compute Betti-0 curves from minimum spanning trees.

    Parameters
    ----------
    mst_asd : list of 3-tuples (int, int, float)
        sorted list of (source, target, weight) tuples corresponding to
        the edges in the minimum spanning tree of SCG for ASD subjects.
    mst_tdc : list of 3-tuples (int, int, float)
        sorted list of (source, target, weight) tuples corresponding to
        the edges in the minimum spanning tree of SCG for TDC subjects.

    Returns
    -------
    beta_asd: list of Betti-0 numbers (number of connected components) for ASD
    beta_tdc: list of Betti-0 numbers (number of connected components) for TDC
    thresholds: list of threshold at which Betti-0 numbers are computed
    """
    w_asd = [e[2] for e in mst_asd]
    w_tdc = [e[2] for e in mst_tdc]
    # thresholds = sorted(np.unique(w_tdc + w_asd))
    thresholds = np.linspace(1.0/10000, 1.0, 10000)
    beta_asd = list()
    beta_tdc = list()
    for t in thresholds:
        beta_asd.append(1+sum([1 for w in w_asd if w <= t]))
        beta_tdc.append(1+sum([1 for w in w_tdc if w <= t]))

    return(beta_asd, beta_tdc, thresholds)
