"""Container for functions and classes regarding interaction with stdin and stdout"""

import sys

####################
#       stdin      # 
####################

def yes_or_not( message, default ):
  """asks the question in *message* and check if the input is *y* or *n*
  Parameter
  ---------
  message: string
    question for the user
  default: 'y'|'n'
    default value in case the input is left empty
  output
  ------
  True|False if input from stdin is 'y'|'n'
  """
  #check the default value
  if( default == 'y' ):
    choise = ([y]/n)
  elif( default == 'n' ):
    choise = (y/[n])
  else:
    print( "The possible choises are 'y' or 'n' only" )
    exit()
  #if message has a '?' at the end remove it
  message = message.strip()
  if( message[-1] == '?' ):
    message = message[:-1]
  #ask the question and parse the answer
  while True:
    from_stdin = raw_input( "{0} {1}? " ).strip()
    if( len(from_stdin)==0 ):  #if the input is empty, set default
      from_stdin = default
    #check and return if std matches
    if( from_stdin == 'y' ):
      return True  
    elif( from_stdin == 'n' ):
      return False


####################
#      stdout      # 
####################


def printer(data):
  """
  Print things to stdout on one line dynamically
  Parameters
  ----------
  data: 
    what has to go on screen
  """
  sys.stdout.write("\r\x1b "+data.__str__())
  sys.stdout.flush()
