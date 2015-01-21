#!/bin/bash

# for aozalevsky:
#   read comments!

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
    echo $1
    if [[ $1 == "USER"* ]]
    then
        echo 'user'
        exit 1 # show error text to user
    elif [[ $1 == "ADMIN"* ]]
    then
    if [[ $1 == *"PY"* ]]
    then
        echo 'taisniqm'
        exit 1 # mail to taisniqm@gmail.com
    else
        echo 'aoz'
        exit 1 # mail to aozalevsky
    fi
    elif [[ $1 == "Unable"* ]]
    then
        echo 'aoz' # minimization/md probem. ?
        exit 1 #mail to aozalevsky
    else
        echo 'taisniqm'
        exit 1 # mail to taisniqm@gmail.com
    fi        
    exit 1
}


function python_wrapper {
  err=$(${1} 2>&1)
  error=$( echo "${err}" | grep -o "Exception: .*")
 if [[ -n $error ]]
 then
    error=$(echo ${error} | sed 's/Exception: //')
    error_exit "${error}"
  elif [[ -n $error ]]  
    then 
      error_exit "${error}"
 elif [[ -n $err ]]
    then
    error_exit "${err}"
  fi
}

function topol_error {
  error=$(${1} 2>&1 | egrep '^ERROR.*topol|Fatal error' )
    if [[ $error == *"topol"* ]]
    then
        echo 'ffbonded error'
        echo $error
        exit 1 # mail to taisniqm@gmail.com
    elif [[ -n $error ]]
    then
        error_exit "Unable to perform energy minimization (run 1)"
    fi
}

function stage {
  echo "<STAGE: ${1}>"
}


function mkmd {
grompp -f ${GMXLIB}/md-vacuum${1}.mdp -c em2 -p topol -o md -maxwarn 1

mdrun -deffnm md$nt -pd -v &
pid=$!
trap "kill $pid 2> /dev/null" EXIT
while kill -0 $pid 2> /dev/null; do
    if [[ -n $(grep 'nan ' md.log -m 1) ]]
    then
    kill $pid
    rm md.log
    if [[ $1 < 5 ]]
    then
        mkmd $(($1+1))
    else
       error_exit 'BANG'
    fi
    fi
    sleep 5
done

trap - EXIT
}

#set -e

stage 1

job=${1%.json}

lattice=$(python val_lattice.py ${job}.json) || error_exit "USER: Error in input json file"

if [[ $lattice == 'u' ]]
    then
    echo 'Lattice type required' # ask user for lattice type
    exit 1
fi

if [[ $# != 3 ]]
then
    if [[ $# != 4 ]]
    then
        exit 1
    else
    cp ../${3} ../${4} .
    python_wrapper "python json2cosm.py -i ${job}.json -o ${job}_tmp.pdb -r ${job}_r -t ${job}_t -l ${lattice} --seq ${3} --oligs ${4} -m ${job}_map"
    fi
else
python_wrapper "python json2cosm.py -i ${job}.json -o ${job}_tmp.pdb -r ${job}_r -t ${job}_t -l ${lattice} -m ${job}_map"
fi

stage 2

if [[ "${2::1}" = "h" ]]
then 
ccg="calori"
else
ccg="calori_sq"
fi
if [[ "$2" == "0" ]]
then
nt=""
else
nt=" -nt $2"
fi

python_wrapper "python plot.py ${job}_map ${job}_map.pdf"
python movechains.py ${job}_tmp.pdb > ${job}_cosm.pdb || error_exit "Double particles movement error"
if grep -Fq "5' - 3' ends" ${job}_r
then
    sctype=circular
else
    sctype=linear
fi

pdb2gmx -f ${job}_cosm.pdb -o beg.gro -merge all -ff ${ccg} -water none || error_exit  "Unable to create topology file"
awk -v name=${job} '/; Include Position restraint file/{print "#include \"" name "_r\""}1' topol.top > topol_tmp 
mv topol_tmp topol.top 

stage 3

topol_error "grompp -f ${GMXLIB}/minl.mdp -c beg -p topol -o em -maxwarn 1"
mdrun -deffnm em -v || error_exit "Unable to perform energy minimization (run 1)"

stage 4


grompp -f ${GMXLIB}/min-implicit.mdp -c em -p topol -o em2 -maxwarn 1 || error_exit "Unable to perform energy minimization (run 2)"
mdrun -deffnm em2 -pd${nt} -v || error_exit "Unable to perform energy minimization (run 2)"
editconf -f em2.gro -o ${job}_end.pdb -conect
echo 1 | trjconv -f em -s em -o ${job}_multi.pdb -conect -skip 10 || error_exit "Unable to prepare output files"

stage 5
if [[ $# == 3 ]]
then
    python_wrapper "cosm2full.py -i ${job}_end.pdb -t ${job}_t -l ${lattice} -p ${3} -o ${job}_end_full.pdb"
else
    python_wrapper "cosm2full.py -i ${job}_end.pdb -t ${job}_t -l ${lattice} -s ${4} -o ${job}_end_full.pdb"
fi
python_wrapper "python appendcnct.py ${job}"

# ------- continue button (MD) ------------

stage 6

mkmd 0 || error_exit "Unable to perform MD"

stage 7

editconf -f md.gro -o ${job}_mdend.pdb 
echo 1 | trjconv -f md -s md -o ${job}_mdmulti.pdb -conect -skip 500 || error_exit "Unable to prepare output files"

stage 8

echo $sctype
echo $lattice

exit $?
