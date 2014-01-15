"""This library contains functions and command line arguments that are commong among codes"""

import argparse as ap
import io_custom as ioc
import itertools as it
import my_functions as mf

#formatter_class that allows raw description and argument default
class RawDescrArgDefHelpFormatter(ap.RawDescriptionHelpFormatter,
        ap.ArgumentDefaultsHelpFormatter):
    """
    Enable features of 'RawDescriptionHelpFormatter' and 'ArgumentDefaultsHelpFormatter' united
    1) both description and epilogue are already correctly formatted and should not be line-wrapped
    2) automatically adds information about default values to each of the argument help messages
    """
#formatter_class that allows raw text and argument default
class RawTextArgDefHelpFormatter(ap.RawTextHelpFormatter,
        ap.ArgumentDefaultsHelpFormatter):
    """
    Enable features of 'RawTextHelpFormatter' and 'ArgumentDefaultsHelpFormatter' united
    1) both description and epilogue are already correctly formatted and should not be line-wrapped
    2) maintains whitespace for all sorts of help text, including argument descriptions
    3) automatically adds information about default values to each of the argument help messages
    """

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

def insert_or_replace(p, print_def=False):
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

def pandas_group(p, description=None):
    """
    Add group with pandas keys
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
    pandas: 
        group containing pandas options
    """
    pandas = p.add_argument_group(title="Pandas", description=description)
    pandas.add_argument("--pandas", action="store_true", 
            help="Use `pandas.read_table` instead of `numpy.loadtxt` to read the files")
    pandas.add_argument("--chunks", action="store", type=int,
            help="If pandas used, read the input file in chunks of '%(dest)s' lines")
    return p, pandas

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

def int_or_list(string):
    """
    Check if the input is a int or list of ints
    Parameters
    ----------
    string: string
        string to parse

    output: list of ints
        return a list of integers
    """
    splitted_string = string.split()
    return [to_int(ss) for ss in splitted_string]

def to_int(string):
    """
    Convert the string into integer. Raise an error if fails
    Parameters
    ----------
    string:
        input string to evaluate

    output
    ------
    integer: int
    """
    try:
        int_string = int(string)
        return int_string
    except ValueError:
        msg = "{0} is not an integer".format(string)
        raise ap.ArgumentTypeError(msg)

def int_or_str(string):
    """
    Check if the input can be converted to int. 
    If yes returs an integer, otherwise a string
    Parameters
    ----------
    string: string
        string to parse

    output: int or string
    """
    try: 
        return int(string)
    except ValueError:
        return string

def outfile(string):
    """
    Avoid overwriting output files
    Parameters
    ----------
    string: string
        string to parse

    output: string
    """
    if ioc.file_exists(string):
        raise ap.ArgumentTypeError("file '{}' alread exists".format(string))
    else:
        return string

# ==========================
# Action classes and functions
# ==========================

class StoreFmt(ap.Action):
    """
    Check the lenghts of the format list. If has one element, is transformed to
    a string 
    """
    def __call__(self, parser, namespace, values, option_string=None):
        if(len(values) == 1):
            setattr(namespace, self.dest, values[0])
        else:
            setattr(namespace, self.dest, values)

class StoreCycle(ap.Action):
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
                msg='''argument "{f}" requires between {nmin} and {nmax}
                arguments'''.format(f=self.dest,nmin=nmin,nmax=nmax)
                raise ap.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength

def required_range(nmin,nmax):
    """
    Check if the values of the arguments are withing the given limits
    Make sure that the 'type' is numeric (not string)

    Parameters
    ----------
    nmin, nmax: numbers
        minimum and maximum values for the arguments

    output
    ------
    RequredRange class
    """
    def check_value(v, nmin, nmax):
        """Check that nmin>=v>=nmax
        Throw a ArgumentTypeError execption if this is not true
        """
        if i<nmin and i>nmax:
            msg='''argument "{f}" requires numbers between {nmin}
            and {nmax} '''.format(f=self.dest,nmin=nmin,nmax=nmax)
            raise ap.ArgumentTypeError(msg)

    class RequiredRange(ap.Action):
        def __call__(self, parser, args, values, option_string=None):
            try:   # assume that values is a list and loop trough it
                for v in values:
                    check_value(v, nmin, nmax)
            except TypeError:  #if it's not a list, pass the value directly
                check_value(values, nmin, nmax)

            setattr(args, self.dest, values)
    return RequiredRange

def file_exists(warning=False, remove=False):
    """
    Check that files exist. Build similarly to argparse.FileType

    Keyword Arguments, used only if more than one files provided:
        - warning -- bool. If *True* throw a warning instead of an error if file
          does not exists
        - remove -- bool. If *True* remove the non existing files

    output
    ------
    FileExists class
    """
    class FileExists(ap.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            # Try to open the file name
            to_be_removed = [] #list of non existing file names
            # message to send to the error handling 
            message = "can't open '%s'"
            # if only one file name passed
            if isinstance(values, str):
                if not ioc.file_exists(values):
                    raise ap.ArgumentTypeError(message % (values))
            else:
                for fn in values:
                    if not ioc.file_exists(fn):
                        if warning is False:
                            raise ap.ArgumentTypeError(message % (fn))
                        else:
                            print("File '{0}' does not exists".format(fn))
                            to_be_removed.append(fn)
                if remove is True:  # remove non existing files
                    for tbr in to_be_removed:
                        values.remove(tbr)
            setattr(namespace, self.dest, values)
    return FileExists

def multiple_of(multiple, reshape=False):
    """
    Check that an argument is given 'N' times multiple

    Parameters
    ----------
    multiple: integer
        the number of arguments must be multiple of this
    reshape: bool
        if given and len(values) == N*multiple, the list is reshape to
        len(value)=N and len(value[i]) =multiple

    output
    ------
    MultipleOf class 
    """
    class MultipleOf(ap.Action):
        def __call__(self, parser, args, values, option_string=None):
            if len(values)%multiple != 0:
                msg='''The number of argument of "{f}" must be multiple of {m}
                '''.format(f=self.dest,m=multiple)
                raise ap.ArgumentTypeError(msg)
            if reshape:
                output = [values[multiple*i:multiple*(i+1)] for i in
                        range(len(values)//multiple)]
            else:
                output = values
            setattr(args, self.dest, output)
    return MultipleOf

class Cm2Inch(ap.Action):
    "convert the input from cm to inces. Useful for matplotlib figure sizes"
    def __call__(self, parser, namespace, values, option_string=None):
        from myplotmodule import inc2cm
        values = [v/inc2cm for v in values]
        setattr(namespace, self.dest, values)

