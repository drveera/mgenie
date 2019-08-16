
plinkscript=$2
pgenbase=$3
outprefix=$4

cp $pgenbase.bim $1
r2=0.08
lcount=$(< $1 wc -l)
while [ $lcount -gt 40000 ]
do
    r2=$(echo "($r2 - 0.005)" | bc -l)
    echo $r2
    $plinkscript --bfile $pgenbase \
                 --indep-pairwise 5000 500 $r2 \
                 --out $outprefix
    lcount=$(< $1 wc -l)
done

