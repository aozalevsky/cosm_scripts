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

set -e

if [ -d "$1" ]; then
    rm -f -r $1_swp
    mv $1 $1_swp
fi
mkdir $1
cp $1.json $1
cd $1

if [ $# != 3 ]
then
   if [ $# != 5 ]
    then
        exit 1
    else
    cp ../$4 .
    cp ../$5 .
    python ../json2cosm.py -i $1.json -o $1.pdb -r $1_r -t $1_t -l $2 --seq $4 --oligs $5
    fi
else
python ../json2cosm.py -i $1.json -o $1.pdb -r $1_r -t $1_t -l $2
fi

### Step 2

rm -f cosm002.ff
rm -f residuetypes.dat
if [ "${2::1}" = "h" ]
then 
ccg="cosm002"
else
ccg="cosm002sq"
fi
ln -s ../${ccg}.ff .
ln -s ${ccg}.ff/residuetypes.dat .
if [ "$3" == "0" ]
then
nt=""
else
nt=" -nt $3"
fi


pdb2gmx -f $1.pdb -o beg.gro -merge all -ff ${ccg} -water none
awk -v name=$1 '/; Include Position restraint file/{print "#include \"" name "_r\""}1' topol.top > topol_tmp
mv topol_tmp topol.top 
grompp -f ${ccg}.ff/minl.mdp -c beg -p topol -o em
mdrun -deffnm em -v
grompp -f ${ccg}.ff/min-implicit.mdp -c em -p topol -o em2
mdrun -deffnm em2 -pd$nt -v
grompp -f ${ccg}.ff/md-vacuum.mdp -c em2 -p topol -o md
mdrun -deffnm md$nt -pd -v
echo 1 | trjconv -f md -s md -o $1_end.pdb -conect -b 20000
echo 1 | trjconv -f md -s md -o $1_md.pdb -conect -skip 500
cd ../
cp $1/$1_end.pdb .
if [ $# == 3 ]
then
    python cosm2full.py -i $1_end.pdb -t $1/$1_t -l $2 -p 50 -o $1_end_full.pdb
else
    python cosm2full.py -i $1_end.pdb -t $1/$1_t -l $2 -s $5 -o $1_end_full.pdb
fi
#sed -i '$ d' $1_end.pdb
python appendcnct.py $1
