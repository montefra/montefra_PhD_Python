"""Collection of utilities for the analysis of the power spectrum"""

import numpy as np

#================================================================
# measured power spectrum                                        
#================================================================

class PS_header(object):
    """
    The power spectrum files have a header containing 
    #        sum(w)          sum(w^2n(z))    sum(w^2)
    #data   75288.1803      3.4065  19024.8771
    #random 3398210.1471    132.4419        749371.2610
    This class reads the header, stores the values and provide few methods to
    get alpha, normalisation and shot noise
    Parameters
    ----------
    f: string or file object
        name of the power spectrum file 
    """
    def __init__(self, f):
        def _read(fo):
            """
            read the interesting parts from the file object fo
            """
            self.data = fo.readline()  #first line ignored
            self.data = [float(s) for s in fo.readline().split()[1:]]
            self.random = [float(s) for s in fo.readline().split()[1:]]

        try:
            with open(f, 'r') as fo: 
                _read(fo)
        except TypeError: # if it is a file object
            initial_position = f.tell()  # get the current file position
            f.seek(0)  # go back to the beginning of the file
            _read(f)
            f.seek(initial_position) # put back the file in the previous position


    def get_alpha(self):
        """
        returns the value of alpha=sum_{data}(w)/sum_{ran}(x)
        """
        return self.data[0]/self.random[0]

    def get_N2randoms(self):
        """
        Returns the square of the normalisation computed from the randoms
        N^2 = alpha*sum_{ran}(w^2n(z))
        """
        return self.get_alpha()*self.random[1]
    def get_N2data(self):
        """
        Returns the square of the normalisation computed from the data 
        N^2 = sum_{data}(w^2n(z))
        """
        return self.data[1]

    def get_SNrandoms(self):
        """
        Get the shot noise from the randoms alone:
        P_{sn} = alpha*(alpha+1)*sum_{ran}(w^2)
        """
        alpha = self.get_alpha()
        return alpha*(1+alpha)*self.random[2]
    def get_SNdata(self):
        """
        Get the shot noise from the data and the randoms:
        P_{sn} = sum_{data}(w^2) + alpha**2*sum_{ran}(w^2)
        """
        alpha = self.get_alpha()
        return self.data[2] + alpha**2*self.random[2]
# end class PS_header(object):

def average_bins(k, pk, n_modes, merge_nbins):
    """
    Average the power spectrum 'merge_bins' at a time using 'n_modes' as weights
    Parameters
    ----------
    k, pk, n_modes: 1D arrays
        wavenumber, power spectrum and number of modes per k-bin
    merge_bins: int
        number of k-bins to merge to obtain the output
    output
    outk, outpk: 1D arrays:
        averaged wavenumber and power spectrum in larger k-bins
    """
    pk_size = pk.size
    if k.size != pk_size and n_modes != pk_size:
        ValueError("k, pk and n_modes must have the same size")

    out_bins = pk_size//merge_nbins  # number of block of size 'merge_nbins'
    remaining_bins = merge_nbins*out_bins - pk_size

    # convert k, pk and n_modes into a merge_nbins*out_bins 2D array
    if remaining_bins == 0:
        k2D = k.reshape([out_bins, merge_nbins])
        pk2D = pk.reshape([out_bins, merge_nbins])
        n_modes2D = n_modes.reshape([out_bins, merge_nbins])
    else:
        k2D = k[:remaining_bins].reshape([out_bins, merge_nbins])
        pk2D = pk[:remaining_bins].reshape([out_bins, merge_nbins])
        n_modes2D = n_modes[:remaining_bins].reshape([out_bins, merge_nbins])
    
    #average of k and weighted average of pk along the 1st axis
    outk = np.average(k2D, axis=1)
    outpk = np.average(pk2D, axis=1, weights=n_modes2D)
    if remaining_bins != 0:
        outk = np.r_[outk, np.average(k[remaining_bins:])] #add the last elements
        outpk = np.r_[outpk, np.average(pk[remaining_bins:], weights=n_modes[remaining_bins:])] #add the last elements
    
    return outk, outpk                  



#================================================================
# model power spectrum                                            
#================================================================
class WinmatException(Exception):
    pass
class PkException(Exception):
    pass

def pk_model(plin, p1loop, parameters, k=None, winmat=None):
    """
    Return the model power spectrum 
    P(k)= b^2*(exp(-k^2/kstar^2)*P_{lin}(k) + A_MC*P_{1loop}(k))

    Parameters
    ----------
    plin: 2D numpy array
        kl, P_{lin}(kl). If *k* nor *winmat* are given, the output 
        is evaluated at *kl*
    p1loop: 2D numpy array
        k1, P_{1loop}(k1)
    parameters: list
        model parameters: b, kstar, A_MC [others if needed]
    k: 1D numpy array
        if given the model is evaluated in 'k'
    winmat: list of numpy arrays
        if given the model is convolved with the window function.
        if given it must contain two or 4 arrays
        winmat = kj, Wij [, W0j, G02i].
        If the last two are given, the integral constraint is computed

    output
    ------
    modelpk: numpy array
        1xN array with N=kl.shape, if *k* and *winmat* are *None*,
        N=k.shape if *k* is given and *winmat* is *None*,
        N=winmat[1].shape[0] if *winmat* is not *None*

    All interpolations are performed with *numpy.interp*. The interpolated 
    values outside the desired range are set to 0
    """

    do_plin_interp = True # do the linear power spectrum
    # get the wavenumber where to evaluate the model
    if winmat is None:
        if k is None:
            do_plin_interp = False
            k = plin[:,0] 
    else:
        nwinmat = len(winmat) # get the number of elements in winmat
        if nwinmat != 2 and nwinmat != 4:
            raise WinmatException("The list *winmat* must contain 2 or 4 numpy arrays instead of {}".format(nwinmat))
        k = winmat[0]

    # interpolate the linear and 1loop power spectra
    if do_plin_interp:
        plin_interp = np.interp(k, plin[:,0], plin[:,1], left=0, right=0)
    else:
        plin_interp = plin[:,1]
    p1loop_interp = np.interp(k, p1loop[:,0], p1loop[:,1], left=0, right=0)

    b, kstar, AMC = parameters[:3] #get the parameter
    modelpk = b**2 * (np.exp(-(k/kstar)**2)*plin_interp + AMC*p1loop_interp)   #model power spectrum

    # 
    if winmat is not None:
        if nwinmat == 2:
            wG = None
        else:
            wG = winmat[2:]
        modelpk = convolve(modelpk, winmat[1], wG=wG)
    return modelpk
# end def pk_model(plin, p1loop, parameters, k=None, winmat=None):

def convolve(pk, wij, k=None, wG=None):
    """
    convolve ps with wij matrix.
    Parameters
    ----------
    pk: 1D array
        power spectrum to convolve. size L or N
    wij: 2D array
        NxM matrix with the actual window matrix
    k: list of two 1D arrays
        if ps.size == L, then k[0].size==L, k[1].size==M.
        `pk` is interpolated in k[1] 
    wG: list of two 1D arrays
        W0j, G02i: if given the integral constrain computed
        
    output
    ------
    convolved_ps: 1D array of dimension N
    """
    # check k
    if k is not None:
        if len(k) != 2:
            raise PkException("The list *k* must have two elements")
        if k[0].size != pk.size:
            raise PkException("*pk* and *k[0]* must have the same size")
        if k[1].size != wij.shape[1]:
            raise WinmatException("""The size of *k[1]* and of the second dimension\
                    of *wij* must be the same""")

    if wG is not None:
        if len(wG) != 2:
            raise WinmatException("The list *wG* must contain 2 numpy arrays")
        if wG[0].size != wij.shape[1]:
            raise WinmatException("""The size of *wG[0]* and of the second dimension\
                    of *wij* must be the same""")
        if wG[1].size-1 != wij.shape[0]:
            raise WinmatException("""The size of *wG[1]* must be larger by 1 than the\
                    first dimension of *wij*""")

    if k is not None:
        pk_interp = np.interp(k[1], k[0], pk, left=0, right=0)
    else:
        pk_interp = pk

    # convolve the power spectrum with the window matrix
    pk_convolved = np.dot(wij, pk_interp)
    # apply the integral constraint
    if wG is not None:
        w0j, g20i = wG
        pk_convolved -= np.sum(w0j*pk_interp)*g20i[1:]/g20i[0]

    return pk_convolved
#end def convolve(pk, wij, k=None, wG=None):


