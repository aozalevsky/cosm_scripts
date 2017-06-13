#!/bin/bash

export X3DNA=/home/www-data/web2py/applications/cosm/private/x3dna-v2.1
export PATH=${PATH}:${X3DNA}/bin
T="/tmp/$(uuidgen).pdb"
fiber -b -seq=${1} ${T}
find_pair ${T} stdout | analyze stdin
x3dna_utils cp_std BDNA
# sed 's/ 35.*$/ 34.3/' -i bp_step.par
sed 's/ 35.*$/ 33.75/' -i bp_step.par
sed 's/3\.3../3.200/' -i bp_step.par
rebuild -atomic bp_step.par ${2}
rm -rf ${T}
