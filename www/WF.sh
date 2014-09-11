############################
# Required: json file only; creates separate directory for the structure
#           cosm.ff and cosm8.ff in home directory
#   
# Examples: 
#    origami.json (hexagonal lattice, no sequence no mdrun -nt) --> molecular dynamics result
#        bash WF.sh origami h 0
#    origami.json (square lattice, with sequence, mdrun -nt 4) --> molecular dynamics result
#        bash WF.sh origami s 4
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

### Step 1

python ../json2cosm.py -i $1.json -o $1.pdb -r $1_r -t $1_t -l $2 

### Step 2

rm -f cosm.ff
rm -f residuetypes.dat
if [ "${2::1}" = "h" ]
then 
ccg="cosm"
else
ccg="cosm8"
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
mdrun -deffnm em
grompp -f ${ccg}.ff/min-implicit.mdp -c em -p topol -o em2
mdrun -deffnm em2 -pd$nt
grompp -f ${ccg}.ff/md-vacuum.mdp -c em2 -p topol -o md
mdrun -deffnm md$nt -pd
#trjconv -f md -s md -o md.pdb -skip 500
