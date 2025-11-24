#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
echo
echo "> Compiling codes..."
g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/rho2ome_sp.cxx -o ../bin/rho2ome_sp -lgsl -lgslcblas -lm -larmadillo
g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/nanoshell_ome_al_p3.cxx -o ../bin/oap -lgsl -lgslcblas -lm -larmadillo
echo "> ...Done!"
echo
none=1.e-30

read  a  Dome  ome21  G  omemi  omema  metal  model  gain_model  solvent E0 rap host < ../data/input/nanosphere_eV.dat
if [ $metal == "gold" ]; then
    ex1=1.6;
    ex2=2.8;
    omi=1.6;
    oma=2.8;
elif [ $metal == "silver" ]; then
    ex1=2.0;
    ex2=4.5;
    omi=2.2;
    oma=3.4;    
fi

echo "> Removing existing directory ../data/output/rho"
rm -fr ../data/output/rho
echo "> Creating directory ../data/output/rho"
mkdir ../data/output/rho

for rho in 0.4 0.5 0.6 0.7 0.8; do
    echo "> rho = $rho"
    fro=(`./../bin/rho2ome_sp $rho`)
    omeG=${fro[0]}
    Gth=${fro[1]} 
    GG=`echo "1.2*$Gth"|bc -l`
    echo "> Setting the gain at the frohlich frequency ome_sp = "$omeG" eV"
    echo "> and setting the gain at 1.2*Gth = "$G"..."
    echo $a $Dome $omeG $GG $omi $oma $metal $model $gain_model $solvent $E0 $rho $host> ../data/input/nanosphere_eV.dat
    ./../bin/oap
    cp "../data/output/oGp/ome_p3.dat" "../data/output/rho/$rho.dat"
    echo $omeG > "../data/output/rho/omeB-$rho.dat"
    bash nano-shell-sketch.bash $rho
    mv "../img/nanoshell.png" "../data/output/rho/$rho.png"
done

echo "> Resetting the original input file..."
echo $a $Dome $ome21 $G $omemi $omema $metal $model $gain_model $solvent $E0 $rap $host> ../data/input/nanosphere_eV.dat
echo "> ...Done!"

echo "> Producing the image..."
gnuplot "[nanoph-2024-0491]_rho_emi.gp"
echo "image ready in img/output/"
exit

