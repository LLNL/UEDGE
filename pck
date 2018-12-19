#!/bin/csh -f
#
# PCK - simple generated script to dispatch to correct architecture make
#

set HostArch = `$BASIS_ROOT/bin/cfgman use`
set Pact     = `awk '($1 == "PACTRoot") { print $3 }' dev/$HostArch/include/configured`

cd $HostArch
$Pact/bin/pact ARCH=$HostArch $argv
exit($status)
