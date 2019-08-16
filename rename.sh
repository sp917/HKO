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
        for muT in 0 0.25 0.3 0.4 0.5 0.75
        do
            for muQV in -0.2 -0.1 0 0.1 0.25 -0.25 0.5
            do
                mv psfc_modified_UVmu${mu}_Tmu${muT}_QVmu${muQV}_${sn}d01_tendency.png psfc_modified_UVmu${mu}0_Tmu${muT}_QVmu${muQV}_${sn}d01_tendency.png
            done
        done
    done
    cd ..
done

cd
