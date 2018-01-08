#!/bin/bash

#set -v -x -e

#set -e -x

# for aozalevsky:
#   read comments!

BASE=$(dirname $(readlink -f ${0}))
export PATH=${BASE}:${PATH}
export GMXLIB=${BASE}/static/gromacs
export TEMPLATE_PATH=${BASE}/static/templates
export MPLCONFIGDIR=${BASE}/.matplotlib
export PYTHON_EGG_CACHE=${BASE}


function error_exit {
    case ${1} in
    USER*)
        echo -ne "<ERROR: ${1}>"
        exit 1 # show error text to user
    ;;
    ADMIN*)
    if [[ ${1} == *"PY"* ]]
    then
        echo -ne "<ERROR: ${1}>"
        exit 2 # mail to taisniqm@gmail.com
    else
        echo -ne "<ERROR: ${1}>"
        exit 2 # mail to aozalevsky
    fi
    ;;
    Unable*)
        echo -ne "<ERROR: ${1}>" # minimization/md probem. ?
        exit 2 #mail to aozalevsky
    ;;
    *)
        echo -ne "<ERROR: ${1}>"
        exit 2 # mail to taisniqm@gmail.com
    ;;
esac
}

function stage {
  echo -ne "<STAGE: ${1}>"
}

#set -e

NT=''

while [[ $# > 1 ]]
do
key="$1"

case ${key} in
    -j|--json)
    json="$2"
    shift # past argument=value
    ;;
    -nt|--num_threads)
    ### It's better not to interfere with GROMACS heuristic. It knows better 
    ### how to handle load balancing
    NT="-nt $2"
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


stage 'Modelling'

grompp -f ${GMXLIB}/md.mdp -c em2 -p topol -o md -maxwarn 1 || error_exit "Unable to prepare md run"
mdrun -deffnm md -v -pd ${NT} || error_exit "Unable to perform md"

stage 'Preparing output'

rm -f ${job}_md_trj.pdb && cp ${job}_aem.pdb ${job}_md_trj.pdb

echo 1 | trjconv -f md -s md -o ${job}_tmp.pdb -conect -skip 50

cat ${job}_tmp.pdb >> ${job}_md_trj.pdb && rm -f ${job}_tmp.pdb

rm -rf ${job}_dynamics && mkdir ${job}_dynamics
cp ${job}_cosm.pdb ${job}_dynamics/.
cp ${job}_map.pdf ${job}_dynamics/. 
cp ${job}_aem.pdb ${job}_dynamics/.
cp ${job}_em_trj.pdb ${job}_dynamics/.
cp ${job}_md_trj.pdb ${job}_dynamics/.
cp ${job}_aem_full.pdb ${job}_dynamics/.

rm -f ${job}_dynamics.zip && zip -r ${job}_dynamics.zip ${job}_dynamics && rm -rf ${job}_dynamics

exit $?
