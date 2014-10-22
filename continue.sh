#!/bin/bash

# Continue MD from last successfull tpr file
# append - append to last trajectory
# -cpi - use checkpoint
# -v - be verbose

BASE=/home/domain/silwer/work/cosmo-web/applications/test/private/cosm-web
export GMXLIB=${BASE}/static/gromacs
export PATH=${BASE}:${PATH}
export TEMPLATE_PATH=${BASE}/static/templates

job=${1%.json}


mdrun -s md.tpr -append -cpi md.cpi -v

echo 1 | trjconv -f md.gro -s md -o ${job}_end.pdb -conect
echo 1 | trjconv -f md -s md -o ${job}_md.pdb -conect -skip 500

if [ $# == 4 ]
then
    cosm2full.py -i ${job}_end.pdb -t ${job}_t -l $2 -p $4 -o ${job}_end_full.pdb
else
    cosm2full.py -i ${job}_end.pdb -t ${job}_t -l $2 -s $5 -o ${job}_end_full.pdb
fi
appendcnct.py ${job}
