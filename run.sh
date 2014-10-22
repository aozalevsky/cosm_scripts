#!/bin/bash

############################
# Required: json file only; creates separate directory for the structure
#           cosm.ff and cosm8.ff in home directory
#   
# Examples: 
#    origami.json (hexagonal lattice, no sequence no mdrun -nt) --> molecular dynamics result
#        bash WF.sh origami h 0
#    origami.json (square lattice, with sequence, mdrun -nt 4) --> molecular dynamics result
#        bash WF.sh origami s 4 seq origami.csv
#
# WARNING:
#   replaces $1 directory!
############################

BASE=/home/domain/silwer/work/cosmo-web/applications/cosm/private/cosm-web
export GMXLIB=${BASE}/static/gromacs
export PATH=${BASE}:${PATH}
export TEMPLATE_PATH=${BASE}/static/templates


function mkmd {
grompp -f ${GMXLIB}/${ccg}.ff/md-vacuum${1}.mdp -c em2 -p topol -o md -maxwarn 1

mdrun -deffnm md$nt -pd -v &
pid=$!
trap "kill $pid 2> /dev/null" EXIT
while kill -0 $pid 2> /dev/null; do
    # Do stuff
    if [[ -n $(grep 'nan ' md.log -m 1) ]]
    then
#    echo 'bang'
    kill $pid
    rm md.log
    if [[ $1 < 5 ]]
    then
        mkmd $(($1+1))
    else
        echo 'no more md.mdp'
        exit 1
    fi
#    exit 1
#    else 
#    echo "ok $1"
    fi
    sleep 5
done

trap - EXIT
}

job=${1%.json}

set -e

if [ $# != 4 ]
then
   if [ $# != 5 ]
    then
        exit 1
    else
    #seq=$(basename $4)
    #cp $4 $seq
    #oligs=$(basename $5)
    #cp $5 $oligs
    json2cosm.py -i ${job}.json -o ${job}.pdb -r ${job}_r -t ${job}_t -l $2 --seq $4 --oligs $5
    fi
else
json2cosm.py -i ${job}.json -o ${job}.pdb -r ${job}_r -t ${job}_t -l $2
fi

### Step 2

if [ "${2::1}" = "h" ]
then 
ccg="cosm002"
else
ccg="cosm002sq"
fi
if [ "$3" == "0" ]
then
nt=""
else
nt=" -nt $3"
fi


pdb2gmx -f ${job}.pdb -o beg.gro -merge all -ff ${ccg} -water none
awk -v name=${job} '/; Include Position restraint file/{print "#include \"" name "_r\""}1' topol.top > topol_tmp
mv topol_tmp topol.top 
grompp -f ${GMXLIB}/${ccg}.ff/minl.mdp -c beg -p topol -o em
mdrun -deffnm em -v
grompp -f ${GMXLIB}/${ccg}.ff/min-implicit.mdp -c em -p topol -o em2
mdrun -deffnm em2 -pd$nt -v

mkmd 0

echo 1 | trjconv -f md -s md -o ${job}_end.pdb -conect -b 20000
echo 1 | trjconv -f md -s md -o ${job}_md.pdb -conect -skip 500

if [ $# == 4 ]
then
    cosm2full.py -i ${job}_end.pdb -t ${job}_t -l $2 -p $4 -o ${job}_end_full.pdb
else
    cosm2full.py -i ${job}_end.pdb -t ${job}_t -l $2 -s $5 -o ${job}_end_full.pdb
fi
appendcnct.py ${job}
