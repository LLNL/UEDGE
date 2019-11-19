#
# Python file with some parallel operations
#
parallel_version = "$Id: parallel.py,v 7.1 2019/11/01 22:35:38 meyer8 Exp $"

from Numeric import *
from types import *
# --- Try import mpi - if not found, then run in serial mode
try:
  import mpi
  me = mpi.rank
  npes = mpi.procs
except ImportError:
  me = 0
  npes = 0

if npes > 0: lparallel = 1
else:        lparallel = 0

# --- The interface has changed some in the newest version of pyMPI.
# --- Check the interface to the mpi.recv command. The newer versions
# --- return a tuple instead of just the data itself.
# --- Is there a better way of doing this?
if lparallel:
  mpi.send(me,me)
  _i = mpi.recv(me)
  if type(_i) == TupleType: _newpympi = 1
  else:                     _newpympi = 0
else:
  _newpympi = 1

if _newpympi:
  def mpirecv(pe=0,ms=0):
    result,stat = mpi.recv(pe,ms)
    return result
else:
  def mpirecv(pe=0,ms=0):
    return mpi.recv(pe,ms)

# ---------------------------------------------------------------------------
# --- By default, set so that only PE0 sends output to the terminal.
if lparallel:
  mpi.synchronizeQueuedOutput('/dev/null')

def number_of_PE():
  if not lparallel: return 0
  return mpi.procs
def get_rank():
  if not lparallel: return 0
  return mpi.rank

# ---------------------------------------------------------------------------
# Enable output from all processors
def EnableAll():
  if not lparallel: return
  mpi.synchronizeQueuedOutput(None)

# ---------------------------------------------------------------------------
# Disable output from all but processor 0
def DisableAll():
  if not lparallel: return
  mpi.synchronizeQueuedOutput('/dev/null')

# ---------------------------------------------------------------------------
# Print object on all processors
def pprint(obj):
  if not lparallel:
    print(str(obj))
    return
  # Ignore all exceptions to make sure that there is not a lock up.
  try:
    ss = str(obj)
  except:
    ss = ''
  mpi.synchronizedWrite(ss+'\n')
  if mpi.rank == 0: print()

# ---------------------------------------------------------------------------
# Print array (or list) from all processors
def aprint(obj):
  if not lparallel:
    print(str(obj))
    return
  mpi.synchronizedWrite(str(obj))

# ---------------------------------------------------------------------------
# Get address of processor
def self_address():
  if not lparallel: return 0
  return mpi.rank

# ---------------------------------------------------------------------------
# Copy an array from processor i to processor 0
def getarray(src,v,dest=0):
  if not lparallel: return v
  if mpi.rank == src:
    mpi.send(v,dest)
  elif mpi.rank == dest:
    return mpirecv(src)
  return v

# ---------------------------------------------------------------------------
# Gather an object from all processors into a list
def gather(obj,dest=0):
  if not lparallel: return [obj]
  if mpi.rank == dest:
    result = []
    for i in range(mpi.procs):
      if i == dest:
        result.append(obj)
      else:
        result.append(mpirecv(i))
    return result
  else:
    mpi.send(obj,dest)
    return [obj]

# ---------------------------------------------------------------------------
# Define a barrier
def barrier():
  if not lparallel: return
  mpi.barrier()

# ---------------------------------------------------------------------------
# Broadcast an object to all procerssors
def broadcast(obj,root=0):
  if not lparallel: return obj
  return mpi.bcast(obj,root)

# ---------------------------------------------------------------------------
# Gather an object from all processors into a list and scatter back to all
def gatherall(obj):
  if not lparallel: return obj
  obj = gather(obj)
  return mpi.bcast(obj)
    
# ---------------------------------------------------------------------------
# General gatherarray which returns an array object combining the
# first dimension.
def gatherarray(a,root=0,othersempty=0,bcast=0):
  if not lparallel: return a
  # --- First check if input can be converted to an array
  isinputok = 1
  try:
    if type(a) in [type(0.),type(0)]:
      a = array([a])
    else:
      a = array(a)
  except:
    isinputok = 0
  # --- Make sure the input is ok on all of the processors
  isinputok = globalmin(isinputok)
  # --- If any returned an error, then all exit (to avoid a deadlock)
  if not isinputok:
    print("Object could not be converted to an array")
    return None
  # --- Now, actually gather the array.
  # --- The check of whether the result is ok may not be needed.
  try:
    result = gather(a,root)
    isinputok = 1
  except:
    isinputok = 0
  # --- Make sure again that the input is ok on all of the processors
  isinputok = globalmin(isinputok)
  if not isinputok:
    print("Error in gather object")
    if type(a) == ArrayType: print("Object has shape ",shape(a))
    return None
  # --- All processors but root simply return either the input argument
  # --- or an empty array unless the result is to be broadcast
  if me != root and not bcast:
    if othersempty: return zeros(len(shape(a))*[0],a.typecode())
    else: return a
  # --- Root processor reshapes the data, removing the first dimension
  # --- Do it bit by bit since the data passed by the other processors may
  # --- not be all the same size.
  if me == root:
    newlen = 0
    for i in range(npes):
      newlen = newlen + shape(result[i])[0]
    newshape = list(shape(result[0]))
    newshape[0] = newlen
    newresult = zeros(newshape,a.typecode())
    i1 = 0
    for i in range(npes):
      i2 = i1 + shape(result[i])[0]
      newresult[i1:i2,...] = result[i]
      i1 = i2
  else:
    newresult = 0
  if bcast: newresult = mpi.bcast(newresult,root)
  return newresult
  ## --- Old way
  ## --- Its easy if all of the arrays passed from the other processors
  ## --- are the same size.
  #result = array(result)
  #ss = list(shape(result))
  #snew = ss[1:]
  #snew[0] = ss[0]*ss[1]
  #result.shape = snew
  ## --- Return the result
  #return result


# ---------------------------------------------------------------------------
# Find the nonzero value of array over all processors. This assumes that the
# non-zero values for each index are the same for all processors.
# Resulting data is broadcast to all processors.
def parallelnonzeroarray(a):
  dmax = parallelmax(a)
  dmin = parallelmin(a)
  result = where(not_equal(dmax,0),dmax,dmin)
  return result

# ---------------------------------------------------------------------------
# Generic global operation on a distributed array.
def globalop(a,localop,mpiop,defaultval):
  if type(a) in [FloatType,IntType]:
    local = a
  elif len(a) > 0:
    try:
      # --- Protect application of localop since the data may not be
      # --- appropriate on all processors, for example it may be an empty array.
      local = localop(a)
    except:
      local = defaultval
  else:
    local = defaultval
  if not lparallel: return local
  return mpi.allreduce(local,getattr(mpi,mpiop))

# ---------------------------------------------------------------------------
# Specific operations on a distributed array.
def globalmax(a):
  return globalop(a,max,"MAX",-1.e36)
def globalmin(a):
  return globalop(a,min,"MIN",+1.e36)
def globalsum(a):
  return globalop(a,sum,"SUM",0.)
def globalave(a):
  s = globalop(a,sum,"SUM",0.)
  if type(a) in [FloatType,IntType]: a = [a]
  n = globalsum(len(a))
  if n > 0: return s/n
  else:     return 0.

# ---------------------------------------------------------------------------
# Generic parallel element-by-element operation on a distributed array.
def parallelop(a,mpiop):
  if not lparallel: return a
  if type(a) == type(array([])):
    a1d = ravel(a) + 0
    for i in range(len(a1d)):
      a1d[i] = mpi.allreduce(a1d[i],getattr(mpi,mpiop))
    a1d.shape = shape(a)
    return a1d
  else:
    return mpi.allreduce(a,getattr(mpi,mpiop))


# ---------------------------------------------------------------------------
# Specific parallel element-by-element operations on a distributed array.
def parallelmax(a):
  return parallelop(a,"MAX")
def parallelmin(a):
  return parallelop(a,"MIN")
def parallelsum(a):
  if not lparallel: return a
  return mpi.allreduce(a,mpi.SUM)


