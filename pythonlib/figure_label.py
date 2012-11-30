# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import docstring
import numpy as np

# <codecell>

class MyFigure( Figure ):
    """test set_xlabel and set_ylabel"""
    def __init__(self, *args, **kwargs):
        Figure.__init__(self, *args, **kwargs)

    @docstring.dedent_interpd
    def set_xlabel(self, xlabel, fontdict=None, labelpad=0.08, **kwargs):
        """
        Call signature::

          set_xlabel(xlabel, fontdict=None, labelpad=None, **kwargs)

        Set the label for the xaxis.

        *labelpad* is the spacing in points between the label and the x-axis

        Valid kwargs are :class:`~matplotlib.text.Text` properties:
        %(Text)s

        ACCEPTS: str

        .. seealso::

            :meth:`text`
                for information on how override and the optional args work

        WARNING: as now it works only if used after all the axes/subplots 
        have been created and require fine tuning (there is no integration
        with tight_layout)
        """
        if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
            kwargs['horizontalalignment'] = 'center'
        if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
            kwargs['verticalalignment'] = 'bottom'

        #get the lower and leftmost border of the axes/subplots
        left, bottom, right, up = np.NAN, np.NAN, np.NAN, np.NAN
        for ax in self.get_axes():
            bbox = ax.get_position()
            left = min( bbox.xmin, left )
            right = max( bbox.xmax, right )
            bottom = min( bbox.ymin, bottom )
            up = max( bbox.ymax, up )
        
        #write the xlabel as text below the lowest axes. labelpad should leave enough space
        #for the default xticklabels
        t = self.text( (right+left)/2, bottom-labelpad, xlabel, fontdict=fontdict, **kwargs )
        return t

    @docstring.dedent_interpd
    def set_ylabel(self, ylabel, fontdict=None, labelpad=0.10, **kwargs):
        """
        Call signature::

          set_ylabel(ylabel, fontdict=None, labelpad=None, **kwargs)

        Set the label for the xaxis.

        *labelpad* is the spacing in points between the label and the y-axis

        Valid kwargs are :class:`~matplotlib.text.Text` properties:
        %(Text)s

        ACCEPTS: str

        .. seealso::

            :meth:`text`
                for information on how override and the optional args work

        WARNING: as now it works only if used after all the axes/subplots 
        have been created and require fine tuning (there is no integration
        with tight_layout)
        """

        if( 'rotation' not in kwargs ):   #default rotation
            kwargs['rotation'] = 'vertical'

        if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
            kwargs['horizontalalignment'] = 'left'
        if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
            kwargs['verticalalignment'] = 'center'
        
        #get the lower and leftmost border of the axes/subplots
        left, bottom, right, up = np.NAN, np.NAN, np.NAN, np.NAN
        for ax in self.get_axes():
            bbox = ax.get_position()
            left = min( bbox.xmin, left )
            right = max( bbox.xmax, right )
            bottom = min( bbox.ymin, bottom )
            up = max( bbox.ymax, up )
        
        #write the ylabel as text on the left the leftmost axes. labelpad should leave enough space
        #for the default yticklabels with one or two digits
        t = self.text( left-labelpad, (up+bottom)/2, ylabel, fontdict=fontdict, **kwargs )
        return t

if __name__ == "__main__":

    fig, ax = plt.subplots(nrows=3, ncols=2, sharex=False, sharey=False, figsize=(5,5),
	    FigureClass=MyFigure )

#    for a in ax[:,0]:
#	a.set_ylabel( "test y label", )
#    for a in ax[-1,:]:
#	a.set_xlabel( "test x label", )
    fig.set_xlabel( "test x label", labelpad=0.08 )
    fig.set_ylabel( "test y label", labelpad=0.10 )

    plt.show()

