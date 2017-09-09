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
CLDIR='cluster'


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

stage 'Clustering'

rm -rf ${CLDIR} && mkdir ${CLDIR} 

echo 0 | trjconv -f md -s md -o cluster/md.pdb -sep # -skip 5

cd cluster
 
mpirun -n 4 affbio -m matrix.hdf5 -f *.pdb -s ../${job}_amd.pdb  --nopbc --task load_pdb calc_rmsd prepare_matrix
affbio -m matrix.hdf5 -t calc_median set_preference aff_cluster --factor 120


for i in $(cat aff_centers.out  | cut -f 1)
 do
  grep '^CONECT' ../${job}_cosm.pdb >> ${i}
done

affbio -m matrix.hdf5 -t render --moltype origami --draw_nums --merged_labels -o ${job}_cluster_tier1.png --bcolor --noclear

cd ..

cp ${CLDIR}/${job}_${CLDIR}_tier1_color.png .

stage 'Preparing output'

# Prepare output pack

rm -rf ${job}_clustering && mkdir ${job}_clustering
cp ${job}_cosm.pdb ${job}_clustering/.
cp ${job}_map.pdf ${job}_clustering/. 
cp ${job}_amd.pdb ${job}_clustering/.
cp ${job}_md_trj.pdb ${job}_clustering/.
cp ${job}_amd_full.pdb ${job}_clustering/.
cp ${CLDIR}/${job}_${CLDIR}_tier1_color.png ${job}_clustering/.

mkdir ${job}_clustering/allpics/
mv ${CLDIR}/*.png ${job}_clustering/allpics/
# cp ${job}_cluster_tier2.png ${job}_clustering/.

rm -f ${job}_clustering.zip && zip -r ${job}_clustering.zip ${job}_clustering && rm -rf ${job}_clustering

rm -rf ${CLDIR}

exit $?
