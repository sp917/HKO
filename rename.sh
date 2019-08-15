#!/bin/basj

cd ../../mnt/c/Users/hko/Desktop/plots/Results/nclplots/

for sn in BHD BR1 CCB CCH CP1 CPH CT2 CT3 CT4
do
    echo ${sn}
    cd TS_${sn}
    mv psfc_modified_run_${sn}d01_tendency.png psfc_modified_UVmu0.75_Tmu0.5_QVmu0.5_${sn}d01_tendency.png
    for mu in 0.1 0.2 0.3 0.4 0.5 0.6
    do 
        echo ${mu}
        mv psfc_modified_UVmu${mu}_Tmu0.5_QVmu0.5_${sn}d01_tendency.png psfc_modified_UVmu${mu}0_Tmu0.5_QVmu0.5_${sn}d01_tendency.png
    done
    cd ..
done

cd
