#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
alifol="../img/output/alp_p3_emi/"
kaifol="../img/output/kap_p3_emi/"
compile_requested=false
while getopts ":c" opt; do
    case "$opt" in
        c) compile_requested=true ;;
    esac
done
shift $((OPTIND-1))

if [ "$compile_requested" = true ]; then
    echo
    echo "> Compiling codes..."
    g++ -Wall -I/usr/local/include -L/usr/local/lib ../src/frohlich.cxx -o ../bin/fro -lgsl -lgslcblas -lm -larmadillo
    g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/nanoshell_ome_al_p3.cxx -o ../bin/oap -lgsl -lgslcblas -lm -larmadillo
    g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/nanoshell_omeG_p3.cxx -o ../bin/oGp -lgsl -lgslcblas -lm -larmadillo
    echo "> ...Done!"
    echo
fi
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
fro=(`./../bin/fro`)
omeG=${fro[0]}
Gth=${fro[1]}
echo "> Setting the gain at the frohlich frequency"
echo "> and setting the probe field to "$none"..."
echo $a $Dome $omeG $G $omi $oma $metal $model $gain_model $solvent $none $rap $host> ../data/input/nanosphere_eV.dat
focurr="rho"$rap"-"$host"-at-"$metal"-in-"$solvent"/"
echo "> ...Done!"
echo
echo "> the emission threshold is Gth = "$Gth
echo
echo "> Producing HEAT MAPS"
echo $a $Dome $omeG $G $ex1 $ex2 $metal $model $gain_model $solvent $none $rap $host> ../data/input/nanosphere_eV.dat
./../bin/oGp
echo $a $Dome $omeG $G $omi $oma $metal $model $gain_model $solvent $none $rap $host> ../data/input/nanosphere_eV.dat
ogpfol="../img/output/oGp/"
mnewdir=$ogpfol$focurr
echo "> Removing old files from the image folder: "$ogpfol
rm -fr $mnewdir
echo "> ...Done!"
echo "> Creating image folder: "$mnewdir
mkdir -p $mnewdir
gnuplot emi_maps.gp
mv $ogpfol*.png $mnewdir

echo
echo
echo " Producing SPECTRA"
echo
anewdir=$alifol$focurr
knewdir=$kaifol$focurr
echo "> Removing old files from the image folder: "$alifol
rm -fr $anewdir
echo "> Removing old files from the image folder: "$kaifol
rm -fr $knewdir
echo "> ...Done!"
echo "> Creating image folder: "$anewdir
mkdir -p $anewdir
echo "> Creating image folder: "$knewdir
mkdir -p $knewdir
echo "> Creating image folders: map_isa and map_isk"
mkdir -p $mnewdir"map_isa"
mkdir -p $mnewdir"map_isk"
echo "> ...Done!"
echo 
focurr="rho"$rap"-"$host"-at-"$metal"-in-"$solvent"/"
echo "> Removing old files from the image folder: GIF"
gifdir="../img/output/GIF/"$focurr
rm -fr $gifdir
echo "> Creating image folders: GIF/"$focurr
echo $gifdir
mkdir -p $gifdir
mkdir -p $gifdir"kappa/"
mkdir -p $gifdir"alpha/"
for ii in {0..20..1}; do
    iG=`echo $ii*0.1|bc -l`
    iG=`printf '%.3f\n' $iG`
    iname=`printf "%03d" $ii`
    mG=`echo $iG*$Gth|bc -l`
    echo -e "   > Calculating the emission spectrum for G = "$iG"*Gth\t=\t"$mG
    echo $a $Dome $omeG $mG $ex1 $ex2 $metal $model $gain_model $solvent $E0 $rap $host> ../data/input/nanosphere_eV.dat
    ./../bin/oap
    echo $a $Dome $omeG $mG $omi $oma $metal $model $gain_model $solvent $E0 $rap $host> ../data/input/nanosphere_eV.dat
    gnuplot emi_wh.gp
    mv "../img/output/emi_wh.pdf" $anewdir"/emi_"$iG"Gth.pdf"
    convert -density 300 $anewdir"/emi_"$iG"Gth.pdf" $anewdir"/emi_"$iname".png"
    rm $anewdir"/emi_"$iG"Gth.pdf"
    mv "../img/output/kmi_wh.pdf" $knewdir"/emi_"$iG"Gth.pdf"
    convert -density 300 $knewdir"/emi_"$iG"Gth.pdf" $knewdir"/emi_"$iname".png"
    rm $knewdir"/emi_"$iG"Gth.pdf"
    gnuplot -c emi_maps_step.gp $mG
    mv $ogpfol"map_isa_step.png" $mnewdir"map_isa/emi_"$iname".png"
    mv $ogpfol"map_isk_step.png" $mnewdir"map_isk/emi_"$iname".png"
    convert -resize x940 $anewdir"/emi_"$iname".png" "../img/output/temp0.png"
    convert -resize x940 $mnewdir"map_isa/emi_"$iname".png" "../img/output/temp1.png"
    convert +append "../img/output/temp0.png" "../img/output/temp1.png" $gifdir"alpha/emi_"$iname".png"
    convert -resize x940 $knewdir"/emi_"$iname".png" "../img/output/temp0.png"
    convert -resize x940 $mnewdir"map_isk/emi_"$iname".png" "../img/output/temp1.png"
    convert +append "../img/output/temp0.png" "../img/output/temp1.png" $gifdir"kappa/emi_"$iname".png"
    rm "../img/output/temp0.png" 
    rm "../img/output/temp1.png"
done
convert -delay 20 -loop 0 -dispose previous $gifdir"alpha/*.png"  $gifdir"alpha.gif"
convert -delay 20 -loop 0 -dispose previous $gifdir"kappa/*.png"  $gifdir"kappa.gif"
echo "> Resetting the original input file..."
echo $a $Dome $ome21 $G $omemi $omema $metal $model $gain_model $solvent $E0 $rap $host> ../data/input/nanosphere_eV.dat
echo "> ...Done!"
