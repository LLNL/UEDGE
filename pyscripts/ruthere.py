from numpy import *
import time
import signal
import sys
 
#############################################################################
# From Dave Grote:
# --- Setup signal handler to capture Control-C
# --- To use, first call arminterrupt(). Then at the place where the interrupt
# --- is allowed, call ruthere(). This will raise a KeyboardInterrupt if
# --- Control-C had been pressed.
# --- When a interrupt request is received, all this handler does is set a
# --- flag to that effect. Then, a subsequent call to ruthere will check
# --- that flag, and if set, raise an exception. This allows a graceful
# --- stop with the current time step completed.
# --- Set the following two in case ruthere is called before arminterrupt.
_defaultcontrolC = signal.getsignal(signal.SIGINT)
_controlCrecieved = False
savetracebacklimit = 0

def _handlecontrolC(signum, frame):
    global _controlCrecieved
    _controlCrecieved = True
 
def ruthere(reset=True):
    """
Checks if an interrupt was requested (usually control-C). If so, then raise
an exception. If reset is True, restore the original interrupt handler so that the
calling code does not have to, and so that, if there is an exception, it gets
restored (since the calling code is not returned to).
    """
    global _controlCrecieved
    global _defaultcontrolC
    global savetracebacklimit
    if _controlCrecieved:
        if reset:
            signal.signal(signal.SIGINT, _defaultcontrolC)
        _controlCrecieved = False
        raise KeyboardInterrupt("Interrupt requested")
 
def arminterrupt():
    global _controlCrecieved
    global _defaultcontrolC
    global savetracebacklimit
    _controlCrecieved = False
    _defaultcontrolC = signal.getsignal(signal.SIGINT)
    try:
       savetracebacklimit = sys.tracebacklimit
    except:
       savetracebacklimit = None       
    signal.signal(signal.SIGINT, _handlecontrolC)
    sys.tracebacklimit = 0

def disarminterrupt():
    global _defaultcontrolC
    global savetracebacklimit
    signal.signal(signal.SIGINT, _defaultcontrolC)
    sys.tracebacklimit = savetracebacklimit
  
  
#=========================================================================

arminterrupt()


 
