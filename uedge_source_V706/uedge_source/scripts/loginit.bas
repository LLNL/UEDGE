# Turn on log of basis commands and output.

# See corsica/src/main.cc for an example constructing names from the
# command line, then writing this file and calling basmain.
# See /usr/local/nbasis/include/main.C and /usr/local/nbasis/include/main.m
# for default c++ and Fortran main's.

character*32 codename = "xuedge"  # executable's name, used in version.
#probname             = "uedge"   # basis constructs plot-file name from this
character*64 logname  = trim("uedge")//".log" # probname not yet set fr ex-line
# (could use basnxtsq(nam) to construct next nam in sequence)
integer logunit       = basopen(logname,"w")

baspecho(logunit)

# would be nice to have (requires customized main) ---
# logunit << "> # Command name:"
# logunit << "> #   /...(name of executable goes here)"
# logunit << "> # Options:"
# logunit << "> #   rest of command-line goes here"
# logunit << " "

mdef end()= call glbwrlog(logunit);{end} mend   # write cpu,io,sys,mem to log

#---------------------------------------------------------------------
# There follow two methods for putting version info at top of log and plotfile.
# Comment-out one or both as desired.
#
# 1st method
#---------------------------------------------------------------------
# time,date,vers sent to log & cgm-file (ok if don't mind always getting cgm):
#1## output graphics; version; output tty
#
# 2nd method
#---------------------------------------------------------------------
# Write version info to logfile only:
output .null; version; output tty  # to .null only. not echoed to log.
# Put quotes around each line in file .null and write into file .nullq:
call basisexe("sed -e s/""^#.*""/"" ""/ -e s/""^\(.*\)""/'""'""\1""'""'/ .null > .nullq")

# Write file .nullq to log:
$i=basopen(".nullq","r")
noisy=yes
do $j=1,100;  $i>>$l;  if(eof.eq.1) break;  logunit<<$l;  enddo
noisy=no

# Write version info to plotfile iff cgm-plots are requested:
# ("output graphics;version" goes to log again, though not screen---
# wish could toggle logging off and on in cgm_init below):

# With no devices on, "output graphics;write.." or any plot command init's cgm;
# but with win or ps or ... on, cgm is not init'd until an explicit cgm on.
# The following handles latter case (win;...;cgm on).
# The former requires redef's of all plot commands (!) to begin w. cgm_init.
# This has not been undertaken; there will be no version info in cgm-file in
# this case. Cleanest thing to do is intercept fortran cgm init routine.

define cgm ezcdodev_ command_SSc(SwcS) cgm
#
# Note: the following lines need not be commented out if meth 1 or no version
# info on plots is preferred.
#
integer cgm_on=0
mdef cgm_init=  # wo nf or sf, don't get whole version msg
  output graphics; version; nf; output tty; cgm_on=1
  #output graphics; version; stdplot << return; output tty; cgm_on=1
mend

function ezcdodev_(arg1;arg2,arg3)
  if(exists(".arg2")) then
    if(exists(".arg3")) then;  call ezcdodev(arg1,arg2,arg3)
    else;                      call ezcdodev(arg1,arg2)
    endif
  elseif(exists(".arg3")) then
    << "ezcdodev_: should never be here."; kaboom(0)
  else
    call ezcdodev(arg1)
  endif

  if(exists(".arg2")) then
    if(arg2.eq."close") then; cgm_on=0
    elseif(arg2.eq."on" .and. cgm_on.eq.0) then; cgm_init
    endif
  elseif(cgm_on.eq.0) then; cgm_init # I think default arg2 is "on"
  endif
endf
#---------------------------------------------------------------------
