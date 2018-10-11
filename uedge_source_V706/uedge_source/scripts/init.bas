echo     = no
verbose  = no        # turns off "End of inp...", "Resum..."

# Uncomment the read line to use the features of logfile, executable date, etc
# Will result in plotfiles being uedge.00n.cngm (probname=uedge)
###read loginit.bas

# Uncomment next three lines to redefine real to doubleprec on 32-bit machines
###if(type(albdsi)=="double precision") then
###  read dbprec
###endif

# Chameleons handy for use in functions. Reserve $a-z for user, at prompt.
chameleon $a_, $b_, $c_, $d_, $e_, $f_, $g_, $h_, $i_, $j_, $k_, $l_, $m_,
          $n_, $o_, $p_, $q_, $r_, $s_, $t_, $u_, $v_, $w_, $x_, $y_, $z_

verbose = yes; echo = yes
stdoutnl = 0

