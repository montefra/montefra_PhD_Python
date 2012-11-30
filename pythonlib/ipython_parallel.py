"""collection of functions that enable parallel computation using the parallel
environment in IPython"""

def start_load_balanced_view( to_execute=None ):
  """
  Start a load_balanced_view from IPython.parallel. 
  If the module is not found or ipcluster is not initialised this function does
  nothing and the computation is switched to single serial
  Otherwise the client and the load balanced view are initialised. numpy and
  my_functions are imported in all the engines as np and mf, respectively

  Parameters
  ----------
  to_execute: string or list of strings
    command(s) to execute on all the nodes. Thought for short tasks, like importing modules
  output
  ------
  parallel: bool
    The paraller environment has been set up without problems
  lbview: load_balanced_view object
    use lview.apply( ... ) to apply a job to the balanced view. 
  """
  try:  #try to import Client
    from IPython.parallel import Client, error
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
    message = "Do you want to continue in serial mode"
    if( sio.yes_or_not( message, 'y' ) ):
      return False, None  #disable the paraller computation
    else:
      print( "Start ipcluster and restart this script" )
      exit()

  dview = c[:]   #get the direct view
  #execute the required code on all engines
  if to_execute is not None :
    if( type( to_execute ) == str ): #if it's just a string
      to_execute = [to_execute, ]  #convert to list
    #execute the required commands in block mode, avoids errors if command is slow
    for te in to_execute:
      try:
        dview.execute( te, block=True )  #try to execute the command
      except error.CompositeError, e:  #if an error occurs, print a single one, not one per engine
        e.raise_exception()

  lbview = c.load_balanced_view() #and the load view

  return True, lbview  #return true and the balanced view
#end def start_load_balanced_view( )


def advancement_jobs( lbview, jobs, enginesid, update=30,
    init_status=None):
  """Print the advancement of the jobs in the queue.  
  This functions returns when all jobs are finished
  Parameters
  ----------
  lbview: load_balanced_view object
    scheduler that runs the jobs
  jobs: list of AsyncResult objects
    list of jobs submitted to the task scheduler
  enginesid: list of int
    ids of the engines used in the computation
  update: float or int
    update the status every 'update' seconds. If negative, only the initial and
    final status are written
  init_status: dict
    dictionary returned from load_balanced_view.queue_status(). If given the
    number of jobs per processors is returned
  """

  import numpy as np
  tot_jobs = len(jobs)
  print( "Starting {0} jobs using {1} engines".format( 
    tot_jobs, len(enginesid) ) ) #start message
  if( update > 0 ):  #if: advancement status
    import stdin_stdout as sio
    while not lbview.wait( jobs=jobs, timeout=update ):
      status = lbview.queue_status()
      #get the number of running jobs
      totrunning = np.sum( [status[i]['tasks'] for i in enginesid ] )
      tot_torun = status['unassigned'] 
      already_run = tot_jobs - (totrunning + tot_torun) 
      percentage_run = already_run / float(tot_jobs)
      #print the status message
      message = """{0:.1%} done. {1} finished {2} running, {3} pending.""".format( 
        percentage_run, already_run, totrunning , tot_torun ) 
      sio.printer( message )
    sio.printer( "Finished" )
  else:
  #end if: advancement status
    lbview.wait( jobs=jobs )  #wait until it finishes
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
