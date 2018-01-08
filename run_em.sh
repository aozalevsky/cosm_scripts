#!/bin/bash

#set -v -x -e

#set -e -x

BASE=$(dirname $(readlink -f ${0}))
export PATH=${BASE}:${PATH}
export GMXLIB=${BASE}/static/gromacs
export TEMPLATE_PATH=${BASE}/static/templates
export MPLCONFIGDIR=${BASE}/.matplotlib

function error_exit {
    echo -ne "<ERROR: ${1}>"
    
    case ${1} in
    USER*)
        exit 1 # show error text to user
    ;;
    ADMIN*)
    if [[ ${1} == *"PY"* ]]
    then
        exit 2 # mail to taisniqm@gmail.com
    else
        exit 2 # mail to aozalevsky
    fi
    ;;
    Unable*)
        exit 2 #mail to aozalevsky
    ;;
    *)
        exit 2 # mail to taisniqm@gmail.com
    ;;
esac
}

function python_wrapper {
  err=$(${1} 2>&1)
  # search for error message in python outpyt
  # with true fix for grep return code 
  error=$(echo "${err}" | grep -o "Exception: .*")
  if [[ ${error} ]]
    then
      error=$(echo ${error} | sed "s/Exception: //g")
      error_exit "${error}"
  fi
}

function stage {
  echo -ne "<STAGE: ${1}>"
}


GC=50
NT=''
lattice='u'

while [[ $# > 1 ]]
do
key="$1"

case ${key} in
    -j|--json)
    json="$2"
    shift # past argument=value
    ;;
    -l|--lattice)
    lattice="$2"
    shift # past argument=value
    ;;
    -nt|--num_threads)
    ### It's better not to interfere with GROMACS heuristic. It knows better 
    ### how to handle load balancing
    nt="-nt $2"
    shift # past argument=value
    ;;
    -gc|--gcratio)
    gc="$2"
    shift # past argument=value
    ;;
    *)
    # unknown option
    ;;
esac
shift
done




stage 'Preparing files'

job=${json%.json}

python_wrapper "util.py ${job}.json"
lattice=$(util.py ${job}.json)

### Select proper forcefield

case "${lattice}" in

    "h")
    ccg="cosm";;

    "s")
    ccg="cosmsq";;

    "*")
    error_exit "Wrong lattice type" ;;

esac

## Check GC ratio

if [[ -z ${gc} ]]
    then 
        gc=${GC}
fi

python_wrapper "json2cosm.py -i ${job}.json -o ${job}_tmp.pdb -r tmp_r -t ${job}_t -l ${lattice} -m ${job}_map"
python_wrapper "cosm_restr.py -i ${job}_tmp.pdb -t ${job}_t -l ${lattice} -o ${job}_r"

grep -h "; 5' - 3' ends" tmp_r >> ${job}_r
rm -f tmp_r

### Very strange way to detect scaffold type. But ok, let it be ###

if grep -Fq "5' - 3' ends" ${job}_r
then
    sctype="circular"
else
    sctype="linear"
fi

# Pass scaffold type up to scheduler
echo -ne "<SCTYPE: ${sctype}>"

# Plot project map
python_wrapper "plot.py ${job}_map ${job}_map.pdf" || error_exit 'Unable to create project scheme'



python_wrapper "movechains.py ${job}_tmp.pdb ${job}_cosm.pdb" || error_exit "Double particles movement error"


pdb2gmx -f ${job}_cosm.pdb -o beg.gro -merge all -ff ${ccg} -water none || error_exit  "Unable to create topology file"
awk -v name=${job} '/; Include Position restraint file/{print "#include \"" name "_r\""}1' topol.top > topol_tmp 
mv topol_tmp topol.top 

stage 'Modelling stage 1'

grompp -f ${GMXLIB}/em1 -c beg -p topol -o em1 -maxwarn 1 || error_exit "Unable to prepare em run"
mdrun -deffnm em1 -v || error_exit "Unable to perform energy minimization (run 1)"

grompp -f ${GMXLIB}/em2 -c em1 -p topol -o em2 -maxwarn 1 || error_exit "Unable to prepare em run"
mdrun -deffnm em2 -pd -v || error_exit "Unable to perform energy minimization (run 2)"


echo 0 | trjconv -s em1 -f em1.trr -o em1.xtc && rm em.trr

stage 'Modelling stage 2'

grompp -f ${GMXLIB}/md.mdp -c em2 -p topol -o md -maxwarn 1 || error_exit "Unable to prepare md run"
mdrun -deffnm md -v -pd ${NT} || error_exit "Unable to perform md"

stage 'Preparing output'

editconf -f md.gro -o ${job}_amd.pdb
grep '^CONECT' ${job}_cosm.pdb >> ${job}_amd.pdb

python_wrapper "cosm2full.py -i ${job}_amd.pdb -t ${job}_t -l ${lattice} -p ${gc} -o ${job}_amd_full.pdb"

cp ${job}_cosm.pdb ${job}_md_trj.pdb

echo 1 | trjconv -f em1 -s em1 -o ${job}_tmp.pdb -conect -skip 50

cat ${job}_tmp.pdb >> ${job}_md_trj.pdb && rm ${job}_tmp.pdb

echo 1 | trjconv -f md -s md -o ${job}_tmp.pdb -conect -skip 50

cat ${job}_tmp.pdb >> ${job}_md_trj.pdb && rm -f ${job}_tmp.pdb

rm -rf ${job}_dynamics && mkdir ${job}_dynamics
cp ${job}_cosm.pdb ${job}_dynamics/.
cp ${job}_map.pdf ${job}_dynamics/. 
cp ${job}_amd.pdb ${job}_dynamics/.
cp ${job}_md_trj.pdb ${job}_dynamics/.
cp ${job}_amd_full.pdb ${job}_dynamics/.

rm -f ${job}_dynamics.zip && zip -r ${job}_dynamics.zip ${job}_dynamics && rm -rf ${job}_dynamics

exit $?
