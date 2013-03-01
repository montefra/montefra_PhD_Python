"""Collection of utilities for the analysis of the power spectrum"""


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
    fname: string
        name of the power spectrum file 
    """
    def __init__(self, fname):
        with open(fname, 'r') as f: 
            self.data = f.readline()  #first line ignored
            self.data = [float(s) for s in f.readline().split()[1:]]
            self.random = [float(s) for s in f.readline().split()[1:]]

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

