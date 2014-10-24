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

if [[ $(hostname) == 'cosm' ]]
then
    BASE=/home/www-data/web2py/applications/cosm/private/cosm-web
else
    BASE=$HOME/work/cosmo-web/applications/cosm/private/cosm-web
fi
echo $BASE
export GMXLIB=${BASE}/static/gromacs
export PATH=${BASE}:${PATH}
export TEMPLATE_PATH=${BASE}/static/templates


function error_exit {
  echo "<ERROR: ${1}>"
  exit 1
}


function python_wrapper {
  error=$(${1} 2>&1 | grep '^Ex' | sed 's/Exception: //')
  if [[ -n $error ]]
    then 
      echo "<ERROR: ${error}>"
      exit 1
  fi
}

function stage {
  echo "<STAGE: ${1}>"
}


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

#set -e

stage 1

job=${1%.json}

if [[ $# != 4 ]]
then
   if [[ $# != 5 ]]
    then
        exit 1
    else
    #seq=$(basename $4)
    #cp $4 $seq
    #oligs=$(basename $5)
    #cp $5 $oligs
    python_wrapper "json2cosm.py -i ${job}.json -o ${job}.pdb -r ${job}_r -t ${job}_t -l ${2} --seq ${4} --oligs ${5}"
    fi
else
python_wrapper "json2cosm.py -i ${job}.json -o ${job}.pdb -r ${job}_r -t ${job}_t -l ${2}"
fi


### Step 2

stage 2

if [[ "${2::1}" = "h" ]]
then 
ccg="cosm002"
else
ccg="cosm002sq"
fi
if [[ "$3" == "0" ]]
then
nt=""
else
nt=" -nt $3"
fi


pdb2gmx -f ${job}.pdb -o beg.gro -merge all -ff ${ccg} -water none || error_exit  "Unable to create topology file"
awk -v name=${job} '/; Include Position restraint file/{print "#include \"" name "_r\""}1' topol.top > topol_tmp 
mv topol_tmp topol.top 


stage 3
grompp -f ${GMXLIB}/${ccg}.ff/minl.mdp -c beg -p topol -o em || error_exit "Unable to perform energy minimization (run 1)"
mdrun -deffnm em -v || error_exit "Unable to perform energy minimization (run 1)"


stage 4
grompp -f ${GMXLIB}/${ccg}.ff/min-implicit.mdp -c em -p topol -o em2 || error_exit "Unable to perform energy minimization (run 2)"
mdrun -deffnm em2 -pd$nt -v || error_exit "Unable to perform energy minimization (run 2)"


stage 5
mkmd 0 || error_exit "Unable to perform MD"

stage 6
echo 1 | trjconv -f md.gro -s md -o ${job}_end.pdb -conect || error_exit "Unable to prepare output files"
echo 1 | trjconv -f md -s md -o ${job}_md.pdb -conect -skip 500 || error_exit "Unable to prepare output files"


stage 7
if [[ $# == 4 ]]
then
    python_wrapper "cosm2full.py -i ${job}_end.pdb -t ${job}_t -l $2 -p $4 -o ${job}_end_full.pdb"
else
    python_wrapper "cosm2full.py -i ${job}_end.pdb -t ${job}_t -l $2 -s $5 -o ${job}_end_full.pdb"
fi
python_wrapper "appendcnct.py ${job}"

exit $?
