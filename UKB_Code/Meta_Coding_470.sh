# load bigsnpr image
docker load -i r_with_plink.tar.gz

if [ $1 -lt 41 ]
then
       CHR=1
 elif [ $1 -lt 66 ]
then
       CHR=2
 elif [ $1 -lt 87 ]
then
       CHR=3
 elif [ $1 -lt 102 ]
then 
       CHR=4
 elif [ $1 -lt 120 ]
then 
       CHR=5
 elif [ $1 -lt 141 ]
then 
       CHR=6
 elif [ $1 -lt 159 ]
then 
       CHR=7
 elif [ $1 -lt 173 ]
then 
       CHR=8
 elif [ $1 -lt 189 ]
then 
       CHR=9
 elif [ $1 -lt 204 ]
then 
       CHR=10
 elif [ $1 -lt 230 ]
then 
       CHR=11
 elif [ $1 -lt 250 ]
then 
       CHR=12
 elif [ $1 -lt 257 ]
then 
       CHR=13
 elif [ $1 -lt 269 ]
then 
       CHR=14
 elif [ $1 -lt 281 ]
then 
       CHR=15
 elif [ $1 -lt 298 ]
then 
       CHR=16
 elif [ $1 -lt 321 ]
then 
       CHR=17
 elif [ $1 -lt 327 ]
then 
       CHR=18
 elif [ $1 -lt 355 ]
then 
       CHR=19
 elif [ $1 -lt 366 ]
then 
       CHR=20
 elif [ $1 -lt 371 ]
then 
       CHR=21
else
       CHR=22
fi

dx download UKB_PRS:MetaSTAAR/Null_Models/470K_WES/obj.STAAR.UKB.TC.20240707.Rdata
dx download UKB_PRS:UKB_470K_WES_AGDS_uncompressed/ukb.470k.wes.chr${CHR}.pass.annotated.gds
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/Annotation_name_catalog.csv

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript Meta_Coding_470.R ${1}"