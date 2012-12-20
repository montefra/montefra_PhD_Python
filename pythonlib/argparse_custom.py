"""This library contains functions and command line arguments that are commong among codes"""

import argparse as ap
import itertools as it
import my_functions as mf

#The following functions add arguments to argparse
def version_verbose(p, version):
    """
    Add the version and verbose options
    Parameters
    ----------
    p: argparse instance
        object containing the command line parsing arguments
    output
    ------
    p: argparse instance
        object containing the updated command line parsing arguments 
    """
    p.add_argument('--version', action='version', version=version)
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")

    return p

def insert_or_replace1(p, print_def=False):
    """
    Add exclusive group for creating the output file name substituting or inserting
    Parameters
    ----------
    p: argparse instance
        object containing the command line parsing arguments
    print_def: bool (optional)
        include default values in the help of the arguments
        (useful if 'formatter_class' is not set to 'argparse.ArgumentDefaultsHelpFormatter')
    output
    ------
    p: argparse instance
        object containing the updated command line parsing arguments 
    group: 
        group containing mutual exclusive options
    """
    group = p.add_mutually_exclusive_group()
    if( print_def ):
        group.add_argument("-i", "--insert", action="store", nargs=2,
                default=[".z", ".dat"], help="""Output file name created
                inserting '%(dest)s[0]' before '%(dest)s[1]' in the input file
                name. (default: %(default)s)""")
    else:
        group.add_argument("-i", "--insert", action="store", nargs=2,
                default=[".z", ".dat"], help="""Output file name created
                inserting '%(dest)s[0]' before '%(dest)s[1]' in the input file
                name.""")
    group.add_argument("-r", "--replace", action="store", nargs=2, 
            help="""Output file name created replacing '%(dest)s[0]' with
            '%(dest)s[1]' in the input file name. If none of these two options
            given, 'insert' is assumed""") 
    return p, group

def overwrite_or_skip(p):
    """
    Add exclusive group for overwriting or skipping the existing file names
    Parameters
    ----------
    p: argparse instance
        object containing the command line parsing arguments
    output
    ------
    p: argparse instance
        object containing the updated command line parsing arguments 
    group: 
        group containing mutual exclusive options
    """
    group = p.add_mutually_exclusive_group()
    group.add_argument("--overwrite", action="store_true", help="""If given does
            not check if the output file exists or not before saving it.""")
    group.add_argument("--skip", action="store_true", help="""Skip already
            existing output file names. No operation done on the corresponding
            input file.""")
    return p, group


def cosmology_group(p, description=None, print_def=False, h0_def=None):
    """
    Add group with cosmological parameters
    Parameters
    ----------
    p: argparse instance
        object containing the command line parsing arguments
    description: string (optional)
        description of the group
    print_def: bool (optional)
        include default values in the help of the arguments
        (useful if 'formatter_class' is not set to 'argparse.ArgumentDefaultsHelpFormatter')
    h0_def: float (optional)
        default value for h0
    output
    ------
    p: argparse instance
        object containing the updated command line parsing arguments 
    cosmo: 
        group containing the cosmological options
    """
    cosmo = p.add_argument_group(title="Cosmology", 
            description=description)
    if( print_def ):
        cosmo.add_argument("--om", action="store", type=float, default='0.274',
                help='Omega_matter. (default: %(default)s)')
        cosmo.add_argument("--ok", action="store", type=float, default='0.',
                help='Omega_curvature. (default: %(default)s)')
        cosmo.add_argument("--wde", action="store", type=float, default='-1.',
                help='Dark energy equation of state. (default: %(default)s)')
        cosmo.add_argument("--h0", action="store", type=float, default=h0_def,
                help='''Reduced Hubble parameter. (default: %(default)s)''')
    else:
        cosmo.add_argument("--om", action="store", type=float, default='0.274',
                help='Omega_matter.')
        cosmo.add_argument("--ok", action="store", type=float, default='0.',
                help='Omega_curvature.')
        cosmo.add_argument("--wde", action="store", type=float, default='-1.', 
                help='Dark energy equation of state.')
        cosmo.add_argument("--h0", action="store", type=float, default=h0_def,
                help='''Reduced Hubble parameter.''')
    return p, cosmo


def parallel_group(p, description=None):
    """
    Add group with parameters for parallel computation
    Parameters
    ----------
    p: argparse instance
        object containing the command line parsing arguments
    description: string (optional)
        description of the group
    output
    ------
    p: argparse instance
        object containing the updated command line parsing arguments 
    parallel: 
        group containing the parallel options
    """
    parallel = p.add_argument_group(title="Parallel computation", 
            description=description)

    parallel.add_argument("-p", "--parallel", action="store_true",
            help="""Enable parallel computing as in Ipython > 0.12. The
            ipcluster must be already up. Can be started simply typing:
            ipcluster start --n=#engines""")

    parallel.add_argument("-u", "--update", action="store", default=30,
            type=int, help="""Update rate for the queue status. Set to *-1* to
            disable the printout""")

    return p, parallel


# ==========================
# Type classes and functions
# ==========================

def pow2(string):
    """
    Custom argparse type
    This function check if the input is integer and if it a power of 2
    Parameters
    ----------
    string: string
        string to parse

    output: int
        return the integer
    """
    try:
        string = int(string)
    except:
        msg = "%r is not an integer" % string
        raise ap.ArgumentTypeError(msg)
    if( ((string & (string - 1)) == 0) and string > 0 ):
        return string   #return the integer
    else:
        msg = "%r is a power of 2" % string
        raise ap.ArgumentTypeError(msg)

class FileExists(object):
    """
    Check if the file exists
    """
    def __call__(self, parser, args, values, option_string=None):
        if hasattr(values, '__iter__'):
            fnames = [f.name for f in values] 
            [f.close() for f in values]
        else:
            fnames = values.name 
            values.close()

        setattr(namespace, self.dest, fnames)

class FileType(object):
    """Factory for creating file object types

    Instances of FileType are typically passed as type= arguments to the
    ArgumentParser add_argument() method.

    Keyword Arguments:
        - mode -- A string indicating how the file is to be opened. Accepts the
            same values as the builtin open() function.
        - bufsize -- The file's desired buffer size. Accepts the same values as
            the builtin open() function.
    """

    def __init__(self, mode='r', bufsize=-1):
        self._mode = mode
        self._bufsize = bufsize

    def __call__(self, string):
        # the special argument "-" means sys.std{in,out}
        if string == '-':
            if 'r' in self._mode:
                return _sys.stdin
            elif 'w' in self._mode:
                return _sys.stdout
            else:
                msg = _('argument "-" with mode %r') % self._mode
                raise ValueError(msg)

        # all other arguments are used as file names
        try:
            return open(string, self._mode, self._bufsize)
        except IOError as e:
            message = _("can't open '%s': %s")
            raise ArgumentTypeError(message % (string, e))

    def __repr__(self):
        args = self._mode, self._bufsize
        args_str = ', '.join(repr(arg) for arg in args if arg != -1)
        return '%s(%s)' % (type(self).__name__, args_str)

# ==========================
# Action classes and functions
# ==========================

class store_fmt(ap.Action):
    """
    Check the lenghts of the format list. If has one element, is transformed to
    a string 
    """
    def __call__(self, parser, namespace, values, option_string=None):
        if( len(values) ==1 ):
            setattr(namespace, self.dest, values[0])
        else:
            setattr(namespace, self.dest, values)

class store_cycle(ap.Action):
    "Substitute the value or the array of values with a itertool cycle"
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, it.cycle(values))

def required_length(nmin,nmax):
    """
    Check if an argument has a number of values between nmin and nmax

    Parameters
    ----------
    nmin, nmax: integer
        minimum and maximum number of arguments

    output
    ------
    RequredLength class
    """
    class RequiredLength(ap.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin<=len(values)<=nmax:
                msg='argument "{f}" requires between {nmin} and {nmax} arguments'.format(f=self.dest,nmin=nmin,nmax=nmax)
                raise ap.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength


