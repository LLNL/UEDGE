      subroutine glbwrlog(ioun)
 
c This is extracted from basis' subroutine glbdone in files.m in order to
c to get the timing info into the log-file. (At some point, they moved the
c line call baspecho(0) in subroutine glbdone before the call glbwrtim.)
C Note that since ostime is re-called here, the numbers in the log may differ
c slightly from what appears on the screen, and that glbwrtim in fact writes
c only cpu and sys.
 
      real cpu,io,sys,mem
      integer ioun
      call ostime(cpu,io,sys,mem)
      call glbwrtim(ioun,cpu,io,sys,mem)
      return
      end
