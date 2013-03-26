"""collection of functions that enable parallel computation using the parallel
environment in IPython"""

class Load_balanced_view(object):
    """class that implements the initialisation of a ipython parallel
    load_ballance_view performing some check. It also execute allows to execute
    some python command on all engines, submit pieces of code, and print a
    progress log. A cleanup function provided
    """

    def __init__(self, client=None):
        """
        Start a load_balanced_view from IPython.parallel.  If a client is not
        given, checks are run to see if ipcluster exists and engines are
        running If none of this happen the computation is switched to single
        serial.  Otherwise the client and the load balanced view are
        initialised.

        Parameters
        ----------
        *client*: an IPython parallel client
            if *None* a new object created
        """

        self.do_parallel = True   #everything ok 
        try:  #try to import Client
            from IPython.parallel import Client, error
            if client != None:
                self.c = client
            else:
                self.c = Client() 
            self.engines_id = self.c.ids  #get the id of the engines
            self.dview = self.c[:]
            self.lbview = self.c.load_balanced_view() #load view
        except ImportError:  #if the import fails
            print( """Ipython.parallel.Client cannot be imported.\
 Make sure to have Ipython version > 0.11 installed.""" )
            self.do_parallel = self._continue_serial()
        except IOError:  #if Ipython is not present
            print( """The Ipython cluster has not been started start it\
 before executing the code. e.g. 'ipcluster start --n=4'.""")
            self.do_parallel = self._continue_serial()

    def _continue_serial(self):
        """asks if the user wants to continue in serial mode or quit"""
        import stdin_stdout as sio
        message = "Do you want to continue in serial mode"
        if( sio.yes_or_not( message, 'y' ) ):
            return False  #disable the paraller computation
        else:
            exit()

    def is_parallel_enabled(self):
        """Returns *True* if the initialization went fine, othewise *False*
        output
        ------
        *parallel*: bool
            *True* if the paraller environment has been set up without
            problems, *False* otherwise
        """
        return self.do_parallel

    def exec_on_engine( self, code, block=True ):
        """
        Execute the given code on all engines 
        
        Parameters
        ----------
        to_execute: string or list of strings 
            command(s) to execute on all the nodes. Thought for short tasks,
            like importing modules
        block: bool
            whether or not to wait until done to return. default: True
        """
        if isinstance(code, basestring): # if it's a string
            code = [code,]  # convert to list
        # execute the required commands 
        # (better to do in block mode, avoids errors if command is slow)
        for te in code:
            try:
                self.dview.execute(te, block=block)
            except error.CompositeError, e:  # if an error occurs, print a single one, not one per engine
                e.raise_exception()

    def push(self, variables):
        """
        wrapper around dview.push(dict)
        push a list of variables to the ipython engines
        Parameters
        ----------
        variables: dictionary
            dictionary of variables
        """
        self.dview.push(variables)

    def apply( self, f, *args, **kwargs):
        """
        wrapper around 'lview.apply(self, f, *args, **kwargs)'

        Docstring:
            calls f(*args, **kwargs) on remote engines, returning the result.

        This method sets all apply flags via this View's attributes.

        if self.block is False:
            returns AsyncResult
        else:
            returns actual result of f(*args, **kwargs)
        """
        return self.lbview.apply( f, *args, **kwargs )

    def get_queue_status(self):
        """
        get the status of the queue
        """
        return self.lbview.queue_status()  

    def advancement_jobs(self, jobs, update=30, init_status=None):
        """Print the advancement of the jobs in the queue.  
        This functions returns when all jobs are finished

        Parameters
        ----------
        jobs: list of AsyncResult objects
            list of jobs submitted to the task scheduler
        update: float or int
            update the status every 'update' seconds. If negative, only the initial and
            final status are written
        init_status: dict
            dictionary returned from load_balanced_view.queue_status(). If given the
            number of jobs per processors is returned
        """

        import numpy as np
        tot_jobs = len(jobs)
        print( "Starting {0} jobs using {1} engines".format( tot_jobs,
            len(self.engines_id) ) ) #start message
        if(update > 0):  #if: advancement status
            import io_custom as sio
            while not self.wait(jobs=jobs, timeout=update):
                status = self.get_queue_status()
                #get the number of running jobs
                totrunning = np.sum( [status[i]['tasks'] for i in self.engines_id ] )
                tot_torun = status['unassigned'] 
                already_run = tot_jobs - (totrunning + tot_torun) 
                percentage_run = already_run / float(tot_jobs)
                #print the status message
                message = """{0:.1%} done. {1} finished {2} running, {3} pending.""".format(percentage_run, already_run, totrunning, tot_torun) 
                sio.printer( message )
            #end while not lbview.wait( ... )
            sio.printer( "Finished" )
        else: #else if: advancement status
            self.wait( jobs=jobs )  #wait until it finishes
            print( "Finished" )
        #end if: advancement status

        #if details about the jobs per processor are wanted
        print("")
        if( init_status is not None ):
            final_status = self.get_queue_status() #get the final status
            print( "{0:<5}: # processes".format("id") )
            for i in self.engines_id:
                print( "{0:<5}: {1}".format(i,
                    final_status[i]['completed']-init_status[i]['completed'] )
                    )
        # end def advancement_jobs( ... )

    def wait(self, jobs=None, timeout=-1):
        """wrapper around lview.wait(self, jobs=None, timeout=-1)
        waits on one or more `jobs`, for up to `timeout` seconds.

        Parameters
        ----------

        jobs : int, str, or list of ints and/or strs, or one or more AsyncResult objects
            ints are indices to self.history
            strs are msg_ids
            default: wait on all outstanding messages
        timeout : float
            a time in seconds, after which to give up.
            default is -1, which means no timeout

        Returns
        -------

        True : when all msg_ids are done
        False : timeout reached, some msg_ids still outstanding
        """
        return self.lbview.wait(jobs=jobs, timeout=timeout)

    def clear_cache(self):
        """
        clear the cache of the parallel computation to avoid memory overload.
        from: http://mail.scipy.org/pipermail/ipython-user/2012-December/011874.html
        check if something like this will be implemented eventually
        """
        self.c.purge_results('all') #clears controller
        self.c.results.clear()
        self.c.metadata.clear()
        self.dview.results.clear()
        self.lbview.results.clear()
        assert not self.c.outstanding, "don't clear history when tasks are outstanding"
        self.c.history = []
        self.dview.history = []
        self.lbview.history = []
