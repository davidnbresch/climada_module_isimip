#!/bin/bash
# Template bash script

# define paths
run_path=$HOME/euler_jobs/job_isimip2a_runALL


for ghm in CLM DBH H08 JULES-TUC JULES-UoE LPJmL MATSIRO MPI-HM ORCHIDEE PCR-GLOBWB VEGAS VIC WaterGAP ; do
    for forcing in gswp3 princeton watch wfdei ; do
        if [ "${ghm}_${forcing}" == "CLM_gswp3" ]
        then
#            echo Skipping run for ${ghm}_${forcing}
            continue
        fi
        filename=job2a_${ghm//-/_}_${forcing//-/_}
        if [ -e "${filename}.m" ]
        then
#            echo Skipping run for ${ghm}_${forcing} because file exists
            continue
        fi

        echo Lauching run for ${ghm}_${forcing} , file ${filename}.m
        cp job_isimip2a_template.m ${filename}.m
        sed -i -e "s/GHMNAME/$ghm/g" "${filename}.m"
       	sed -i -e "s/FORCINGNAME/$forcing/g" "${filename}.m"
       	sed -i -e "s/MATSIROSUBTRACTION/0/g" "${filename}.m"
        bsub -W 12:00 -R "rusage[mem=20000]" -n 1 matlab -nodisplay -singleCompThread -r "${filename}"

        # Same thing but subtracting MATSIRO flood fraction
        filename=${filename}_mFRCmastiro
        cp job_isimip2a_template.m ${filename}.m
        echo Lauching run for ${ghm}_${forcing} , file ${filename}.m
        sed -i -e "s/GHMNAME/$ghm/g" "${filename}.m"
       	sed -i -e "s/FORCINGNAME/$forcing/g" "${filename}.m"
       	sed -i -e "s/MATSIROSUBTRACTION/1/g" "${filename}.m"
        bsub -W 12:00 -R "rusage[mem=20000]" -n 1 matlab -nodisplay -singleCompThread -r "${filename}"
    done
done

