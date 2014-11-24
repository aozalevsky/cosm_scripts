name=$1

python json2coords.py -i $name.json -l h -t 1 > ${name}_coords1
python json2coords.py -i $name.json -l h -t 2 > ${name}_coords2
python plot.py ${name}_coords1 ${name}_coords2 ${name}_map.pdf
