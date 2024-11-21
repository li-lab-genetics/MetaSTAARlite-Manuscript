# load bigsnpr image
docker load -i r_with_plink.tar.gz

if [ $1 -lt 19 ]
then
       CHR=1
 elif [ $1 -lt 35 ]
then
       CHR=2
 elif [ $1 -lt 46 ]
then
       CHR=3
 elif [ $1 -lt 56 ]
then 
       CHR=4
 elif [ $1 -lt 69 ]
then 
       CHR=5
 elif [ $1 -lt 79 ]
then 
       CHR=6
 elif [ $1 -lt 89 ]
then 
       CHR=7
 elif [ $1 -lt 100 ]
then 
       CHR=8
 elif [ $1 -lt 108 ]
then 
       CHR=9
 elif [ $1 -lt 117 ]
then 
       CHR=10
 elif [ $1 -lt 128 ]
then 
       CHR=11
 elif [ $1 -lt 141 ]
then 
       CHR=12
 elif [ $1 -lt 147 ]
then 
       CHR=13
 elif [ $1 -lt 156 ]
then 
       CHR=14
 elif [ $1 -lt 166 ]
then 
       CHR=15
 elif [ $1 -lt 178 ]
then 
       CHR=16
 elif [ $1 -lt 191 ]
then 
       CHR=17
 elif [ $1 -lt 198 ]
then 
       CHR=18
 elif [ $1 -lt 208 ]
then 
       CHR=19
 elif [ $1 -lt 214 ]
then 
       CHR=20
 elif [ $1 -lt 218 ]
then 
       CHR=21
else
       CHR=22
fi

dx download UKB_PRS:MetaSTAAR/Null_Models/obj.STAAR.UKB.TC.20240530.95055.Rdata
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/ukb.200k.wgs.chr${CHR}.pass.annotated.gds
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/Annotation_name_catalog.csv

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript Meta_ncRNA_95055.R ${1}"