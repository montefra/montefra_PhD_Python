"""collection of functions that enable parallel computation using the parallel
environment in IPython"""

def start_load_balanced_view( ):
  """
  Start a load_balanced_view from IPython.parallel. 
  If the module is not found or ipcluster is not initialised this function does
  nothing and the computation is switched to single serial
  Otherwise the client and the load balanced view are initialised. numpy and
  my_functions are imported in all the engines as np and mf, respectively

  output
  ------
  parallel: bool
    The paraller environment has been set up without problems
  lbview: load_balanced_view object
    use lview.apply( ... ) to apply a job to the balanced view. 
  """
  try:  #try to import Client
    from IPython.parallel import Client
  except ImportError:  #if the import fails
    print( """Ipython.parallel.Client cannot be imported.
	Make sure to have Ipython version > 0.xx installed.
	Serial computation enabled.""" )
    return False, None  #disable the paraller computation

  #if the import works, initialise the client.
  try:
    c = Client()   #start the client
  except IOError:
    print( """The Ipython cluster has not been started
	start it before executing the code. 
	'ipcluster start --n=4'.""")
    import stdin_stdout as sio
    sio.message = "Do you want to continue in seria mode"
    if( yes_or_not( message, 'y' ) ):
      return False, None  #disable the paraller computation
    else:
      print( "Start ipcluster and restart this script" )
      exit

  lbview = c.load_balanced_view() #and the load view
  #import numpy on all engines
  c[:].execute( 'import numpy as np; import my_functions as mf' )  

  return True, lbview  #return true and the balanced view
#end def start_load_balanced_view( )


def advancement_jobs( lbview, tot_jobs, enginesid, update=30,
    init_status=None):
  """Print the advancement of the jobs in the queue.  
  This functions returns when all jobs are finished
  Parameters
  ----------
  lbview: load_balanced_view object
    scheduler that runs the jobs
  tot_jobs: int
    number of jobs submitted to the task scheduler
  enginesid: list of int
    ids of the engines used in the computation
  update: float or int
    update the status every 'update' seconds. If negative, only the initial and
    final status are written
  init_status: dict
    dictionary returned from load_balanced_view.queue_status(). If given the
    number of jobs per processors is returned
  """

  print( "Starting {0} jobs using {1} engines".format( 
    tot_jobs, len(enginesid) ) ) #start message
  if( update > 0 ):  #if: advancement status
    import stdin_stdout as sio
    while not lbview.wait( timeout=update ):
      status = lbview.queue_status()
      #get the number of running jobs
      totrunning = np.sum( [status[i]['tasks'] for i in enginesid ] )
      tot_torun = status['unassigned'] 
      already_run = tot_jobs - (totrunning + tot_torun) 
      percentage_run = already_run / float(tot_jobs)
      #print the status message
      sio.printer( """{0:.1%} ({1} out of {2}) of submitted job finished. {3}
	  running, {4} pending. """.format( percentage_run, already_run,
	    tot_jobs, totrunning , tot_torun ) )
    sio.printer( "Finished" )
  else:
  #end if: advancement status
    lbview.wait()  #wait until it finishes
    print( "Finished" )
  #end if: advancement status
  #if details about the jobs per processor are wanted
  print("")
  if( init_status is not None ):
    final_status = lbview.queue_status() #get the final status
    print( "{0:<5}: # processes".format("id") )
    for i in enginesid:
      print( "{0:<5}: {1}".format(i,
	final_status[i]['completed']-init_status[i]['completed'] ) )
  return
# end def advancement_jobs( ... )
