#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
NGM_ROOT="$(realpath "$ROOT_DIR/extern/nano_geo_matrix")"
CXXFLAGS="${CXXFLAGS:-} -I$NGM_ROOT/include -I$NGM_ROOT/modules -I$NGM_ROOT/modules/cup"

echo "> Compiling code..."
g++ $CXXFLAGS -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib ../src/nanoshell_num.cxx -o ../bin/nsh -lgsl -lgslcblas -lm -larmadillo
echo "> Done!"
omemi=2.615
omema=3.02

N=100

dome=`echo "($omema - $omemi)/$N"| bc -l`
rm ../data/output/spe.dat
for (( i=0; i<=$N; i++ ))
    do 
        ome=`echo "$omemi+$dome*$i"|bc -l`
        output=`./../bin/nsh $ome`
        echo $output 
        echo $output >> ../data/output/spe.dat
done
