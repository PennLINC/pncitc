---
layout: default
title: PNC ITC
nav_order: 1
---

Analyses were run using [Pehlivanova et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5858592/)'s n427 sample and then running restQA exclusions on them (all this information was from Pehlivanova .csvs). The .csvs for this are available in the `samplerecreation` folder in this project's repo.

All final analyses were run in `/cbica/projects/pncitc/ignore`.


### Sample replication

The code for sample replication is: 

```
# RAN THIS LOCALLY DUE TO CLUSTER ISSUES WITH R
setwd("/Users/kahinim/Desktop")

# read the subject demographics
restdatapnc=read.csv('n2416_RestQAData_20170714.csv') # pnc QA for resting-state data
nmel=read.csv('n452_pnc_itc_whole_sample_20160825.csv') # Marieta final subject QA  
z = read.csv('n427_fsSubcortVol.csv')
pncitc=merge(nmel,restdatapnc, by=c('bblid','scanid')) # merge by Ids  
pncitc=merge(pncitc,z, by=c('bblid','scanid')) # merge by Ids  
# select the neccessary variable for screening and further processing
# age, logk, sex, rest exclusion  variables: voxelwise and motion
pncit1 <- data.frame(
  pncitc$bblid,
  pncitc$healthExclude,
  pncitc$scanid,
  pncitc$logk,
  pncitc$ageAtScan,
  pncitc$logAlpha,pncitc$sex,pncitc$race,pncitc$race2,pncitc$restExclude,pncitc$restExcludeVoxelwise,
  pncitc$restNoDataExclude,pncitc$relMeanRMSmotion,pncitc$restNSpikesMotion,pncitc$restNSpikesMotionExclude,pncitc$restRelMeanRMSMotionExclude, pncitc$restNoDataExclude, pncitc$restVoxelwiseCoverageExclude, pncitc$meduCnbGo1, pncitc$feduCnbGo1
)
colnames(pncit1)=c('bblid','healthExclude',
                   'scanid','logk','ageAtScan','logAlpha','sex','race','race2','restExclude','restExcludeVoxelwise',
                   'restNoDataExclude','relMeanRMSmotion','restNSpikesMotion','restNSpikesMotionExclude','restRelMeanRMSMotionExclude','restNoDataExclude', 'restVoxelwiseCoverageExclude', 'Medu', 'Fedu')

pncit1=pncit1[which(pncit1$restExcludeVoxelwise==0),]
pncit1=pncit1[which(pncit1$restNoDataExclude==0),]
pncit1=pncit1[which(pncit1$restRelMeanRMSMotionExclude==0),]
pncit1=pncit1[which(pncit1$restVoxelwiseCoverageExclude==0),]
pncit1=pncit1[which(pncit1$healthExclude==0),]
pncit1=pncit1[which(pncit1$restNSpikesMotionExclude==0),]
pncit1=pncit1[-which(is.na(pncit1$relMeanRMSmotion)),]


# during manual checking, one subject (id:96832) has data points 90 less than 120 expected 
pncit1=pncit1[-which(pncit1$bblid==96832),]

#get the ids of final subjects
ids=data.frame(pncit1$bblid,pncit1$scanid) # get bblid and scanid for futher analyses 

# write out demographics and bblid and scanid
write.csv(ids,'n293_blbid_scanid.csv',row.names = FALSE,quote = FALSE)
pncit1$age=pncit1$ageAtScan/12
write.csv(pncit1,'n293_demographics.csv',row.names = FALSE,quote = FALSE)


```

### 1. CWAS-MDMR

The computation of  CWASMDMR was  done with  `cwasmdr` singularity image (`/cbica/projects/pncitc/cwasmdmr.simg`). We used packages from the connectir project at [https://github.com/czarrar/connectir](https://github.com/czarrar/connectir)

Distance matrix was first computed with the following script, using masks under the `subjectData` folder in this project's repo: 

```
#!/bin/bash
#$ -l h_vmem=320G #QSUB, can take some hours
#$ -l tmpfree=200G
singimage=/cbica/projects/pncitc/cwasmdmr.simg 
scriptdir=/usr/local/bin
mdmrouput=/cbica/projects/pncitc/ignore/cwas293 #output directory
brainmask=/cbica/projects/pncitc/subjectData/PNCgrey.nii.gz # greymatter mask from pnc   
bgim=/cbica/projects/pncitc/subjectData/PNCbrain.nii.gz # pnc template from pnc
imagelist=/cbica/projects/pncitc/subjectData/imageinput_rest3.csv #list of image in nifti # HAD TO RE-GENERATE THIS LIST AS THE FILE PATHS HAD CHANGED
rm  -rf $mdmrouput # remove previous run if exist 
metric=pearson # pearson correlation 
# compute distance matrix
singularity exec -e -B /cbica/projects/pncitc $singimage $scriptdir/Rscript $scriptdir/connectir_subdist.R $mdmrouput --infuncs1=$imagelist --ztransform --automask1  -c 3 -t 4 --brainmask1=$brainmask --method="$metric" --bg=$bgim  --overwrite --memlimit=200

```

The output of distance matrix: `/cbica/projects/pncitc/ignore/cwas293`
   
The  distance matrix  was used for mdmr computation with `logk` as the main factor.
other covariates used are `sex`, `age`,  and `relative rms`:

 ```math  
 distancematrix = f(logk)+relMeanRMSmotion+sex+age
 ```
   
The script used for mdmr computation is as below: 
```
#!/bin/bash
#$ -l h_vmem=300G
#$ -l tmpfree=300G
singularity exec -e -B /cbica/projects/pncitc  \
/cbica/projects/pncitc/cwasmdmr.simg \
/usr/local/bin/Rscript /usr/local/bin/connectir_mdmr.R -i /cbica/projects/pncitc/ignore/cwas293 -f 'logk+relMeanRMSmotion+sex+age' -m /cbica/projects/pncitc/ignore/samplerecreation/n293_demographics.csv --factors2perm='logk' --save-perms -c 5 -t 5  --ignoreprocerror --memlimit=300 logk_motion_sex_age
```

The output is at: `/cbica/projects/pncitc/ignore/cwas293/logk_motion_sex_age`

### 2. Significant clusters from mdmr
The cluster analysis was computed  with the script `scripts/grf_fslcluster.sh`, written based on  [FSL cluster analysis](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Cluster) with  Gaussian Random Field (GRF) theory.

The script `cluster.sh` called `grf_fslcluster.sh`, with `z=3.09`.  


`cluster.sh` reads as below:
```
 #!/bin/bash # NO NEED TO QSUB
dir=/cbica/projects/pncitc
bash grf_fslcluster.sh -i ${dir}/ignore/cwas293/logk_motion_sex_age/zstats_logk.nii.gz  -m ${dir}/ignore/cwas293/mask.nii.gz -t 3.09 -o ${dir}/ignore/cluster_output 
```

while `grf_fslcluster.sh`reads as: 
```
#!/usr/bin/env bash

###################################################################
###################################################################

###################################################################
# Combine all text file output
###################################################################

###################################################################
# Usage function
###################################################################
Usage(){
  echo ""; echo ""; echo ""
  echo "Usage: `basename $0`  grf_fslcluster.sh -i zstat -m mask -t threshold -o output"
  echo ""
  echo "Compulsory arguments:"
  echo "  -i : zstats: compulsory"
  echo "  -m: mask"
  echo "  -o : Output file name"
  echo "       "
  exit 2
}

###################################################################
# Parse arguments
###################################################################
while getopts "i:t:m:o:h" OPTION ; do
  case ${OPTION} in
    i)
      zstat=${OPTARG}
      ;;
    t)
      thresh=${OPTARG}
      ;;
    m)
      mask=${OPTARG}
      ;;
    o)
      outdir=${OPTARG}
      ;;
    h)
      Usage
      ;;
    *)
      Usage
      ;;
    esac
done

###################################################################
# Ensure that all compulsory arguments have been defined
###################################################################
[[ -z ${outdir} ]] && Usage
[[ -z ${zstat} ]] && Usage
[[ -z ${mask} ]] && Usage

###################################################################
# Now run through each file that we find and append it to the output file
###################################################################
 
if [[ -z ${thresh} ]]; then 
   thresh=2.3
   echo "voxel threshold is 2.3 (default)"
fi 

echo " find d and v " 
dv=$(smoothest -z ${zstat} -m ${mask})

id0=$(echo $dv |cut -d' ' -f2)
id1=$(echo $dv |cut -d' ' -f4)
echo " the dlh is ${id0}"
echo "                  "
echo " the number of volume: ${id1}"
echo $thresh
echo $dv
mkdir -p ${outdir}/cluster_Z${thresh}

cluster -i ${zstat} -d ${id0} --volume=${id1} -t ${thresh} -p 0.05 \
   -o  ${outdir}/cluster_Z${thresh}/cluster_Z${thresh} >  \
    ${outdir}/cluster_Z${thresh}/cluster_Z${thresh}.csv

echo "done"
```


The output of cluster masks is at: `/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09`. One cluster was found:

| Cluster Index | Voxels | P | -log10(P) | MAX | MAX X (vox) | MAX Y (vox) | MAX Z (vox) | COG X (vox) | COG Y (vox) | COG Z (vox) |
| ---- | ---- | ---- | ------ | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|1|	12|	0.000103|	3.99|	3.54|	30|	44|	30|	30.9|	43.8|	30.3|


However, the output image of clusters produced was an .img and .hdr which needed to be turned into a nifti format, I used: 

 ```fslchfiletype NIFTI_GZ cluster_Z3.09.img cluster_Z3.09.nii```



### 3. Seed-based correlation 
A mask was generated from the cluster_z3_09.nii, using fslmath at `/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09/mask1` :

(See [https://mandymejia.com/fsl-maths-commands/](https://mandymejia.com/fsl-maths-commands/))

```
fslmaths cluster_Z3.09.nii.gz -thr 1 -uthr 2 mask1/mask1.nii.gz
```

The mask generated was again in .hdr and .img format, so I used `fslchfiletype` as before to turn them into niftis.The mask was upsampled from 4mm to 2mm and was used as a seed for seed-based correlation. This upsampling was done using the pnc_2mm template, found under `subjectData` in this project's repo. 

The code used for this was: 

```
3dresample -master pnc_template_brain_2mm.nii.gz -input mask1.nii.gz -prefix mask1_2mm.nii.gz
```

The path to the resampled seed is: 
`/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09/mask1/mask1_2mm.nii.gz`

The seed-based correlation was computed with the following script, and the `xcpengine.simg` file under `/cbica/projects/pncitc/ignore`:

```
#!/bin/bash #DEFINITELY QSUB THIS
cd /cbica/projects/pncitc/ignore
XCPEDIR=xcpEngine
seedpoint1=/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09/mask1/mask1_2mm.nii.gz

bblid=/cbica/projects/pncitc/demographics/n293_bblid_scandid.csv
image=/cbica/projects/pncitc/subjectData/rest293/
outputdir=/cbica/projects/pncitc/ignore/seedcorrmaps

mkdir -p ${outputdir}

cat $bblid | while IFS="," read -r a b ; 

do
     img=$(ls -f $image/${a}_${b}_rest.nii.gz)
     singularity exec --cleanenv -B /cbica/projects/pncitc/ignore /cbica/projects/pncitc/ignore/xcpengine.simg /xcpEngine/utils/seedconnectivity -i $img -s $seedpoint1 -o $outputdir -p ${a},${b} -k 6 -n mask1\
     #rm -rf $outputdir/seed/mask1/${a}_${b}_connectivity_mask1_seed.nii.gz

done
```

### 4. Linear regression with FSL `flameo` 

Flameo regression computation requires `design`,`contrast` and `group` text files. The script `makeflameodesig.R` was used to make these:

```
# script to make the design matrices for flameo

library(pracma)
demogr=read.csv('/cbica/projects/pncitc/demographics/n293_demographics.csv') 
#logk+relMeanRMSmotion+age+sex 
desigmatlogkonly=cbind(rep(1,293),demogr$logk,demogr$sex,demogr$relMeanRMSmotion,demogr$age)

grp=ones(293,1) # only one group

contrast4=zeros(5,5); 
diag(contrast4)=1; 


write.table(desigmatlogkonly,'/cbica/projects/pncitc/demographics//desigmatlogkonly.txt',sep=' ',quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(contrast4,'/cbica/projects/pncitc/demographics//contrast4.txt',sep=' ',quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(grp,'/cbica/projects/pncitc/demographics//grp.txt',sep=' ',quote = FALSE,row.names = FALSE,col.names = FALSE)
```

The needed files were then copied over into `/cbica/projects/pncitc/ignore/seedcorrmaps`.


I also converted the designlogkmat from .txt to .mat and did the same for grp.txt, using: 

```
Text2Vest desigmatlogkonly.txt desigmatlogkonly.mat
Text2Vest grp.txt grp.grp
Text2Vest contrast4.txt contrast4.con
```
I also created the flameo csvs pointing to the mask1Z_sm6.nii.gz niftis generated in previous step under the same directory, 
called `mask1.csv`  (`/cbica/projects/pncitc/ignore/seedcorrmaps`). 

The flameo linear regression was computed with this script:

```
#!/bin/bash
#$ -l h_vmem=200G
#$ -l tmpfree=200G # this will throw an error about ": $PATH does not agree with $PATH_modshare counter.", but it shouldn't stop the process. 
bblid=/cbica/projects/pncitc/demographics/n293_bblid_scandid.csv #FIRST HALF RUNS QUICKLY, WOULD QSUB BECAUSE THE SECOND PART TAKES A BIT LONGER
imagedir=/cbica/projects/pncitc/ignore/seedcorrmaps/seed
scriptdir=/cbica/projects/pncitc/ignore/seedcorrmaps
outputdir=/cbica/projects/pncitc/ignore/seedcorrmaps/seed
demogdir=/cbica/projects/pncitc/ignore/seedcorrmaps


imagelist1=$scriptdir/mask1.csv

rm -rf $imagelist1


cat $bblid | while IFS="," read -r a b ; 

do 
     img1=$(ls -f $imagedir/mask1/${a}_${b}_connectivity_mask1Z_sm6.nii.gz)
     
     echo $img1 >> $imagelist1

done 


mask=/cbica/projects/pncitc/subjectData/PNCgrey2mm.nii.gz

fslmerge -t ${outputdir}/4Dcopeseed1.nii.gz $(cat $imagelist1)

flameo --copefile=${outputdir}/4Dcopeseed1.nii.gz   --mask=${mask}   --dm=${demogdir}/desigmatlogkonly.mat  --tc=${demogdir}/contrast4.con  --cs=${demogdir}/grp.grp --runmode=flame1 --ld=$outputdir/mask1/logk #SECOND PART, WHICH TAKES LONGER

```



The outputs of flameo regression: 

`/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk`


In the directory, there are zvalues: 

  `zstat1 : average `

  `zstat2 : logk `

  `zstat3 : sex`

  `zstat4 : motion`

  `zstat5 : age`


### 5. Vizualisation of Results - iPython in CBICA

All computations were done in PNC template. For vizualisation, all the nifti files were tranformed to MNI before as below. Transform
files are available in this project's repo under `PNC_transforms`. 

#### Transform to MNI
```
# import all the requirements and hide warnings
import warnings
warnings.filterwarnings("ignore")


import nilearn.plotting as plott
import nilearn.image as img
from nilearn import datasets,surface
import matplotlib.pyplot as plt
from nipype.interfaces.ants import ApplyTransforms

big_fsaverage = datasets.fetch_surf_fsaverage('fsaverage') # for viz 

#registration paramteters

ref='/cbica/projects/pncitc/subjectData/PNC_transforms/MNI152_T1_2mm_brain.nii.gz'
transform1='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_0Warp.nii.gz'
transform2='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_1Affine.mat'
at = ApplyTransforms()
at.inputs.dimension = 3
at.inputs.reference_image = ref
at.inputs.interpolation = 'NearestNeighbor'
at.inputs.default_value = 0
at.inputs.transforms = [transform1, transform2]
at.inputs.invert_transform_flags = [False, False]
flame1dir='/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/'
zstats=['zstat1','zstat2']
viewim=[]
for i in range(len(zstats)):
    at.inputs.input_image = flame1dir + zstats[i]+'.nii.gz'
    at.inputs.output_image = flame1dir + zstats[i]+'MNI.nii.gz'
    at.run()

    
clusterdirectory = '/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09'
zstats=['/mask1/mask1_2mm','/mask1/mask1']
for i in range(len(zstats)):
    at.inputs.input_image = clusterdirectory + zstats[i]+'.nii.gz'
    at.inputs.output_image = clusterdirectory + zstats[i]+'MNI.nii.gz'
    at.run()
```

 a. for clusters and mean of seed-based correlation: 
_Note: had to use `flchfiletype` on the `copeseed` images before running the script, and move the `.hdr` and `.img` files out of the directory/ remove them altogether - having the nifti and img in the same directory can cause an error. Another tip: if you are using the os.system functionality, running the code from the same directory in which the files are in helps avoid "can't open/read file" errors_
```
# import all the requirements and hide warnings
import warnings
warnings.filterwarnings("ignore")


import nilearn.plotting as plott
import nilearn.image as img
from nilearn import datasets,surface
import matplotlib.pyplot as plt
from nipype.interfaces.ants import ApplyTransforms
from nipype.interfaces.fsl import MultiImageMaths,maths, MeanImage
import os

big_fsaverage = datasets.fetch_surf_fsaverage('fsaverage') # for viz 

#registration paramteters

ref='/cbica/projects/pncitc/subjectData/PNC_transforms/MNI152_T1_2mm_brain.nii.gz'
transform1='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_0Warp.nii.gz'
transform2='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_1Affine.mat'
at = ApplyTransforms()
at.inputs.dimension = 3
at.inputs.reference_image = ref
at.inputs.interpolation = 'NearestNeighbor'
at.inputs.default_value = 0
at.inputs.transforms = [transform1, transform2]
at.inputs.invert_transform_flags = [False, False]

output_image = '/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09/mask1/mask1_2mmMNI.nii.gz'
#put it on surface 
img1=img.load_img(output_image)

# plot of mean of seed-based correlation 

# average of all subject 
seedbasedir='/cbica/projects/pncitc/ignore/seedcorrmaps/seed/'
corrtm=['4Dcopeseed1'] # make sure to change to nifti and remove .img
label = corrtm
viewim=[]
meanimage=MeanImage()
for i in range(len(corrtm)):
    meanimage.inputs.in_file=seedbasedir + corrtm[i]+ '.nii.gz'
    meanimage.inputs.dimension='T'
    meanimage.inputs.out_file=seedbasedir +corrtm[i] + 'mean.nii.gz' 
    meanimage.run()
    os.system('fslchfiletype NIFTI_GZ 4Dcopeseed1mean.img 4Dcopeseed1mean.nii.gz')
    os.system('rm -rf 4Dcopeseed1mean.hdr  4Dcopeseed1mean.img')
    at.inputs.input_image = meanimage.inputs.out_file
    at.inputs.output_image = seedbasedir +corrtm[i] + 'meanMNI.nii.gz'
    at.run()
    img1=img.load_img(at.inputs.output_image)
    v= plott.view_img_on_surf(img1, surf_mesh='fsaverage',threshold=0.1,title='meanseedcorr :'+label[i],cmap = 'coolwarm', symmetric_cmap=True) # for mean seed corr

    viewim.append(v)
    
ii = 0    
for x in viewim:
  ii += 1
  x.save_as_html("/cbica/projects/pncitc/ignore/meanseedbasedcorr" + str(ii) + ".html")
 

```

  b. for mask1 (ran in iPython on CBICA): `notebook/flameomask1.ipynb`
  

  ```
# import all the requirements and hide warnings
import warnings
warnings.filterwarnings("ignore")


import nilearn.plotting as plott
import nilearn.image as img
from nilearn import datasets,surface
import matplotlib.pyplot as plt
from nipype.interfaces.ants import ApplyTransforms

big_fsaverage = datasets.fetch_surf_fsaverage('fsaverage') # for viz 

#registration paramteters

ref='/cbica/projects/pncitc/subjectData/PNC_transforms/MNI152_T1_2mm_brain.nii.gz'
transform1='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_0Warp.nii.gz'
transform2='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_1Affine.mat'
at = ApplyTransforms()
at.inputs.dimension = 3
at.inputs.reference_image = ref
at.inputs.interpolation = 'NearestNeighbor'
at.inputs.default_value = 0
at.inputs.transforms = [transform1, transform2]
at.inputs.invert_transform_flags = [False, False]

flame1dir='/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/'
zstats=['zstat1','zstat2']
label=['mean','logk']
viewim=[]
for i in range(len(zstats)):
    output_image = flame1dir + zstats[i]+'MNI.nii.gz' 
    img1=img.load_img(output_image)
    v= plott.view_img_on_surf(img1, surf_mesh='fsaverage',threshold=3.09,vmax=5,title='zstat :'+label[i],cmap = 'coolwarm', symmetric_cmap=True) 
    viewim.append(v)

ii = 0
for x in viewim:
    ii+=1
    x.save_as_html("/cbica/projects/pncitc/ignore/cluster"+str(ii)+".html")
  ```


Images I generated were saved in `/cbica/projects/pncitc/ignore` in the .html format. I ended up using only the logk visualization here in the manuscript. 

### 6. Regional plot of significant regions of logk
 
 
I used R for this: 


```
library(RNifti, lib.loc = '/cbica/projects/pncitc/mehtareplicate')
library(pracma, lib.loc = '/cbica/projects/pncitc/mehtareplicate')
library(ggplot2, lib.loc = '/cbica/projects/pncitc/mehtareplicate')
library(nlme, lib.loc = '/cbica/projects/pncitc/mehtareplicate')
library(visreg, lib.loc = '/cbica/projects/pncitc/mehtareplicate')
library(Matrix, lib.loc = '/cbica/projects/pncitc/mehtareplicate') # and really any other packages that give you issues - as this can happen
mask1=readNifti('/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/zstat2.nii.gz') 

#get the  postive masks
p_m1=mask1; p_m1[p_m1<3.09]=0
mask1=readNifti('/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/zstat2.nii.gz') 

#get the negative masks
n_m1=mask1; n_m1[n_m1>-3.09]=0


writeNifti(p_m1, '/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/p_m1.nii.gz', template = NULL, datatype = "auto", version = 1)
writeNifti(n_m1, '/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/n_m1.nii.gz', template = NULL, datatype = "auto", version = 1)

b=read.csv('/cbica/projects/pncitc/demographics/n293_bblid_scandid.csv',header=FALSE)

#make table 

corrdata=zeros(293,4)

for (i in 1:293) {
  img1=readNifti(paste0('/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/',b[i,1],'_',b[i,2],'_connectivity_mask1Z_sm6.nii.gz')) # flameo output
  datap1=img1[p_m1!=0]
  datam1=img1[n_m1!=0]
  corrdata[i,]=c(b[i,1],b[i,2],mean(datap1),mean(datam1))
}

colnames(corrdata)=c('bblid','scanid','mask1pos','mask1neg')

write.csv(corrdata,'/cbica/projects/pncitc/demographics/n293_meanseedcorr.csv',quote = FALSE,row.names = FALSE)

# merge CSV
x = read.csv('/cbica/projects/pncitc/demographics/n293_meanseedcorr.csv')
y = read.csv('/cbica/projects/pncitc/demographics/n307_demographics.csv') # demographics are right, when merged the n307 will become n293
z = read.csv('/cbica/projects/pncitc/demographics/n2416_RestQAData_20170714.csv')

final1=merge(x,y, by=c('bblid','scanid')) # merge by Ids  
final2=merge(final1,z, by=c('bblid','scanid')) # merge by Ids 
write.csv(final2,'n293_data.csv',quote = FALSE,row.names = FALSE)

# write as .rds
saveRDS(final2, file = "/cbica/projects/pncitc/demographics/my_data.RDS") 

#start plotting
ddata=readRDS('/cbica/projects/pncitc/demographics/my_data.RDS')
poscluster1mask_nologk=lm(mask1pos~age+sex+relMeanRMSmotion,data=ddata)
negcluster1mask_nologk=lm(mask1neg~age+sex+relMeanRMSmotion,data=ddata)

ylab<-"Correlation (z(r))"

ddata$poscluster1resid<-poscluster1mask_nologk$residuals+mean(ddata$mask1pos)
ggplot(ddata,aes(x=logk,y=poscluster1resid)) + geom_smooth(method = 'lm', colour=('#b40101'), fill = "#ef1212",size=2,alpha=.8)+xlim(c(-8.75,-1))+ geom_point() + xlab("Discount Rate (logK)") +ylab(ylab) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=20), axis.title = element_text(size=26), axis.title.y = element_text(margin = margin(t = 0, r = 27, b = 0, l = 0)))
ggsave('/cbica/projects/pncitc/ignore/cluster1pos.png')

ddata$negcluster1resid<-negcluster1mask_nologk$residuals+mean(ddata$mask1neg)

ggplot(ddata,aes(x=logk,y=negcluster1resid)) + geom_smooth(method = 'lm', colour=('#0c3e6d'), fill = "#69abde",size=2,alpha=1)+xlim(c(-8.75,-1))+ geom_point() + xlab("Discount Rate (logK)") +ylab(ylab) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=20), axis.title = element_text(size=26), axis.title.y = element_text(margin = margin(t = 0, r = 27, b = 0, l = 0)))

ggsave('/cbica/projects/pncitc/ignore/cluster1neg.png')
```
 
Finally, I generated the insets seen in the manuscript (these were the niftis written out in the script below, I then projected them to the surface in `ignore` using similar code as in the previous visualization steps:

```
# import all the requirements and hide warnings
import warnings
warnings.filterwarnings("ignore")
import os

import nilearn.plotting as plott
import nilearn.image as img
from nilearn import datasets,surface
import matplotlib.pyplot as plt
from nipype.interfaces.ants import ApplyTransforms

big_fsaverage = datasets.fetch_surf_fsaverage('fsaverage') # for viz 

#registration paramteters

ref='/cbica/projects/pncitc/subjectData/PNC_transforms/MNI152_T1_2mm_brain.nii.gz'
transform1='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_0Warp.nii.gz'
transform2='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_1Affine.mat'
at = ApplyTransforms()
at.inputs.dimension = 3
at.inputs.reference_image = ref
at.inputs.interpolation = 'NearestNeighbor'
at.inputs.default_value = 0
at.inputs.transforms = [transform1, transform2]
at.inputs.invert_transform_flags = [False, False]

dir='/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/'
masks=['n_m1', 'p_m1']
mask = masks
for i in range(len(masks)):
    at.inputs.input_image = dir + masks[i]+'.nii.gz'
    at.inputs.output_image = dir + masks[i]+'MNI.nii.gz'
    img2 = mask[i]+'.img'
    nii = mask[i]+'.nii.gz'
    os.system('fslchfiletype NIFTI_GZ ' + img2 + ' ' + nii)
    os.system('rm -rf *img')
    at.run()

viewim=[]
label = ['neg','pos']

for i in range(len(masks)):
  output_image = dir + masks[i]+'MNI.nii.gz'
  img1=img.load_img(output_image)
  v= plott.view_img_on_surf(img1, surf_mesh='fsaverage',vmax=5,title='inset :'+label[i],cmap = 'coolwarm', symmetric_cmap=True, threshold = 0.09) 
  viewim.append(v)

ii = 0
for x in viewim:
  ii+=1
  x.save_as_html('/cbica/projects/pncitc/ignore/inset'+str(ii)+'.html')

```

Finally, for the glass brain cluster image, I locally ran this on the `/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09/mask1/mask1_2mmMNI.nii.gz` file: 

```
from nilearn import plotting
import nibabel as nib

stat_map_img = '/Users/kahinim/Desktop/mask1_2mmMNI.nii.gz'
plotting.plot_glass_brain(stat_map_img, output_file='/Users/kahinim/Desktop/test.png', display_mode='ortho', colorbar=False, figure=None, axes=None, title=None, threshold='auto', annotate=True, black_bg=False, cmap=None, alpha=0.7, vmin=None, vmax=None, plot_abs=True, symmetric_cbar='auto', resampling_interpolation='continuous')
```
### 7. logk by sex/age interaction

1. I created the directories `cwas293logkbyge` and `cwas293logkbysex` for the interaction results to be written into...
2. I re-ran the mdmr script using the formula `logk*sex` or `logk*age`, keeping other parameters the same.  Here are the scripts:

```
# sex

#!/bin/bash
#$ -l h_vmem=300G
#$ -l tmpfree=300G
singularity exec -e -B /cbica/projects/pncitc  \
/cbica/projects/pncitc/cwasmdmr.simg \
/usr/local/bin/Rscript /usr/local/bin/connectir_mdmr.R -i /cbica/projects/pncitc/ignore/cwas293sex -f 'logk*sex+age+relMeanRMSmotion' -m /cbica/projects/pncitc/samplerecreation/n293_demographics.csv --factors2perm='logk:sex' --save-perms -c 5 -t 5  --ignoreprocerror --memlimit=300 logk_motion_sex_age

# age
#!/bin/bash
#$ -l h_vmem=300G
#$ -l tmpfree=300G
singularity exec -e -B /cbica/projects/pncitc  \
/cbica/projects/pncitc/cwasmdmr.simg \
/usr/local/bin/Rscript /usr/local/bin/connectir_mdmr.R -i /cbica/projects/pncitc/ignore/cwas293age -f 'logk*age+sex+relMeanRMSmotion' -m /cbica/projects/pncitc/samplerecreation/n293_demographics.csv --factors2perm='logk:age' --save-perms -c 5 -t 5  --ignoreprocerror --memlimit=300 logk_motion_sex_age
```
3. For clustering, I created two directories: `cluster_output_sex` and `cluster_output_age`. The scripts I used for clustering were as below: 

```
# for sex
#!/bin/bash # NO NEED TO QSUB
dir=/cbica/projects/pncitc
bash grf_fslcluster.sh -i ${dir}/ignore/cwas293sex/logk_motion_sex_age/zstats_logk:sex.nii.gz  -m ${dir}/ignore/cwas293sex/mask.nii.gz -t 3.09 -o ${dir}/ignore/cluster_output_sex

#for age
#!/bin/bash # NO NEED TO QSUB
dir=/cbica/projects/pncitc
bash grf_fslcluster.sh -i ${dir}/ignore/cwas293age/logk_motion_sex_age/zstats_logk:age.nii.gz  -m ${dir}/ignore/cwas293age/mask.nii.gz -t 3.09 -o ${dir}/ignore/cluster_output_age
```

No significant clusters were found. 

### 8. SES and age effects on behavioral data
1. I obtained maternal level of education from  `/cbica/projects/pncitc/dropbox/pehlivanovaPncItc/subjectData/demoBehavData/n452_pnc_itc_whole_sample_201608256.csv`
2. I merged it with `n293_data.csv` using `R` and saved `n293_full_demo_data.csv` to the `demographics` directory on CBICA: 

```
setwd("/Users/kahinim/Desktop")

# read the subject demographics
demo=read.csv('n452_pnc_itc_whole_sample_20160825.csv')
data=read.csv('n293_data.csv') 
pncitc=merge(demo,data, by=c('bblid','scanid')) # merge by Ids  
write.csv(pncitc,'n293_full_demo_data.csv',quote = FALSE,row.names = FALSE)
```
3. After obtaining this data, I used R to run a bivariate analysis: 

```
age = pncitc$age
medu1 = pncitc$meduCnbGo1
fedu1 = pncitc$feduCnbGo1
edu = (medu1+fedu1)/2
ses = edu
logk =  pncitc$logk.x
reg_age = lm(logk~age)
plot(logk~age)
reg_ses = lm(logk~ses)
plot(logk~ses)
```
I also used the mgcv package to check for nonlinear age effects, as in Pehlivanova et al., 2018, using  logk ~ s(age). This yielded a GCV score of 2.246636.


4. There appeared to be a correlation of -0.155 for SES
5. We conducted sensitivity analysis with SES but the results did not change. The code for this is here: 

#### Sample replication

The code for sample replication is: 

```
setwd("/Users/kahinim/Desktop")

# read the subject demographics
restdatapnc=read.csv('n2416_RestQAData_20170714.csv') # pnc QA for resting-state data
nmel=read.csv('n452_pnc_itc_whole_sample_20160825.csv') # Marieta final subject QA  
z = read.csv('n427_fsSubcortVol.csv')
pncitc=merge(nmel,restdatapnc, by=c('bblid','scanid')) # merge by Ids  
pncitc=merge(pncitc,z, by=c('bblid','scanid')) # merge by Ids  
# select the neccessary variable for screening and further processing
# age, logk, sex, rest exclusion  variables: voxelwise and motion
pncit1 <- data.frame(
  pncitc$bblid,
  pncitc$healthExclude,
  pncitc$scanid,
  pncitc$logk,
  pncitc$ageAtScan,
  pncitc$logAlpha,pncitc$sex,pncitc$race,pncitc$race2,pncitc$restExclude,pncitc$restExcludeVoxelwise,
  pncitc$restNoDataExclude,pncitc$relMeanRMSmotion,pncitc$restNSpikesMotion,pncitc$restNSpikesMotionExclude,pncitc$restRelMeanRMSMotionExclude, pncitc$restNoDataExclude, pncitc$restVoxelwiseCoverageExclude, pncitc$meduCnbGo1, pncitc$feduCnbGo1
)
colnames(pncit1)=c('bblid','healthExclude',
                   'scanid','logk','ageAtScan','logAlpha','sex','race','race2','restExclude','restExcludeVoxelwise',
                   'restNoDataExclude','relMeanRMSmotion','restNSpikesMotion','restNSpikesMotionExclude','restRelMeanRMSMotionExclude','restNoDataExclude', 'restVoxelwiseCoverageExclude', 'Medu', 'Fedu')

pncit1=pncit1[which(pncit1$restExcludeVoxelwise==0),]
pncit1=pncit1[which(pncit1$restNoDataExclude==0),]
pncit1=pncit1[which(pncit1$restRelMeanRMSMotionExclude==0),]
pncit1=pncit1[which(pncit1$restVoxelwiseCoverageExclude==0),]
pncit1=pncit1[which(pncit1$healthExclude==0),]
pncit1=pncit1[which(pncit1$restNSpikesMotionExclude==0),]
pncit1=pncit1[-which(is.na(pncit1$relMeanRMSmotion)),]


# during manual checking, one subject (id:96832) has data points 90 less than 120 expected 
pncit1=pncit1[-which(pncit1$bblid==96832),]

#get the ids of final subjects
ids=data.frame(pncit1$bblid,pncit1$scanid) # get bblid and scanid for futher analyses 

# write out demographics and bblid and scanid
write.csv(ids,'n293_blbid_scanid.csv',row.names = FALSE,quote = FALSE)
pncit1$age=pncit1$ageAtScan/12
write.csv(pncit1,'n293_demographics.csv',row.names = FALSE,quote = FALSE)

pncit1=read.csv('n293_demographics.csv')
medu = pncit1$Medu # average of mat and pat edu
fedu = pncit1$Fedu
edu = (medu+fedu)/2
pncit1$edu = edu
pncit1=pncit1[which(pncit1$edu>0.0),] # remove NaN values
write.csv(pncit1,'n282_demographics.csv',row.names = FALSE,quote = FALSE)
```
I also used the mgcv package to check for nonlinear age effects, as in Pehlivanova et al., 2018, using  logk ~ s(age). This yielded a GCV score of 2.246636.

#### 1. CWAS-MDMR

The computation of  CWASMDMR was  done with  `cwasmdr` singularity image (`/cbica/projects/pncitc/cwasmdmr.simg`). We used packages from the connectir project at [https://github.com/czarrar/connectir](https://github.com/czarrar/connectir)

Distance matrix was first computed with the following script: 

```
#!/bin/bash
#$ -l h_vmem=320G #QSUB, can take some hours
#$ -l tmpfree=200G
singimage=/cbica/projects/pncitc/cwasmdmr.simg 
scriptdir=/usr/local/bin
mdmrouput=/cbica/projects/pncitc/finalreplication/cwas282 #output directory
brainmask=/cbica/projects/pncitc/subjectData/PNCgrey.nii.gz # greymatter mask from pnc   
bgim=/cbica/projects/pncitc/subjectData/PNCbrain.nii.gz # pnc template from pnc
imagelist=/cbica/projects/pncitc/finalreplication/imageinput_rest.csv #list of image in nifti # HAD TO RE-GENERATE THIS LIST AS THE FILE PATHS HAD CHANGED AS WELL AS THE SAMPLE
rm  -rf $mdmrouput # remove previous run if exist 
metric=pearson # pearson correlation 
# compute distance matrix
singularity exec -e -B /cbica/projects/pncitc $singimage $scriptdir/Rscript $scriptdir/connectir_subdist.R $mdmrouput --infuncs1=$imagelist --ztransform --automask1  -c 3 -t 4 --brainmask1=$brainmask --method="$metric" --bg=$bgim  --overwrite --memlimit=200
```

The output of distance matrix: `/cbica/projects/pncitc/finalreplication/cwas282`
   
The  distance matrix  was used for mdmr computation with `logk` as the main factor.
other covariates used are `sex`, `age`, `edu` and `relative rms`:

 ```math  
 distancematrix = f(logk)+relMeanRMSmotion+sex+age+edu
 ```
   
The script used for mdmr computation is as below: 
```
#!/bin/bash
#$ -l h_vmem=300G
#$ -l tmpfree=300G
singularity exec -e -B /cbica/projects/pncitc  \
/cbica/projects/pncitc/cwasmdmr.simg \
/usr/local/bin/Rscript /usr/local/bin/connectir_mdmr.R -i /cbica/projects/pncitc/finalreplication/cwas282 -f 'logk+relMeanRMSmotion+sex+age+edu' -m /cbica/projects/pncitc/finalreplication/samplerecreation/n282_demographics.csv --factors2perm='logk' --save-perms -c 5 -t 5  --ignoreprocerror --memlimit=300 logk_motion_sex_age_edu
```


#### 2. Significant clusters from mdmr

As before,
`cluster.sh` reads as below and calls `grf_fslcluster.sh`:
```
#!/bin/bash # NO NEED TO QSUB
dir=/cbica/projects/pncitc
bash grf_fslcluster.sh -i ${dir}/finalreplication/cwas282/logk_motion_sex_age_edu/zstats_logk.nii.gz  -m ${dir}/finalreplication/cwas282/mask.nii.gz -t 3.09 -o ${dir}/finalreplication/cluster_output
```

while `grf_fslcluster.sh` reads as: 
```
#!/usr/bin/env bash

###################################################################
###################################################################

###################################################################
# Combine all text file output
###################################################################

###################################################################
# Usage function
###################################################################
Usage(){
  echo ""; echo ""; echo ""
  echo "Usage: `basename $0`  grf_fslcluster.sh -i zstat -m mask -t threshold -o output"
  echo ""
  echo "Compulsory arguments:"
  echo "  -i : zstats: compulsory"
  echo "  -m: mask"
  echo "  -o : Output file name"
  echo "       "
  exit 2
}

###################################################################
# Parse arguments
###################################################################
while getopts "i:t:m:o:h" OPTION ; do
  case ${OPTION} in
    i)
      zstat=${OPTARG}
      ;;
    t)
      thresh=${OPTARG}
      ;;
    m)
      mask=${OPTARG}
      ;;
    o)
      outdir=${OPTARG}
      ;;
    h)
      Usage
      ;;
    *)
      Usage
      ;;
    esac
done

###################################################################
# Ensure that all compulsory arguments have been defined
###################################################################
[[ -z ${outdir} ]] && Usage
[[ -z ${zstat} ]] && Usage
[[ -z ${mask} ]] && Usage

###################################################################
# Now run through each file that we find and append it to the output file
###################################################################
 
if [[ -z ${thresh} ]]; then 
   thresh=2.3
   echo "voxel threshold is 2.3 (default)"
fi 

echo " find d and v " 
dv=$(smoothest -z ${zstat} -m ${mask})

id0=$(echo $dv |cut -d' ' -f2)
id1=$(echo $dv |cut -d' ' -f4)
echo " the dlh is ${id0}"
echo "                  "
echo " the number of volume: ${id1}"
echo $thresh
echo $dv
mkdir -p ${outdir}/cluster_Z${thresh}

cluster -i ${zstat} -d ${id0} --volume=${id1} -t ${thresh} -p 0.05 \
   -o  ${outdir}/cluster_Z${thresh}/cluster_Z${thresh} >  \
    ${outdir}/cluster_Z${thresh}/cluster_Z${thresh}.csv

echo "done"
```

The output of cluster masks is at: `/cbica/projects/pncitc/finalreplication/cluster_output/cluster_Z3.09`. We found the same cluster to be stable.


### 7. Spin testing

First, we generated fsaverage5 surfaces from the MNI statmaps we used for visualization locally, using the code below. We also added manually visualized/compared
the images produced here to the images already projected on the surface as a sanity check:

```

library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench/')

library(rgl) #to use ciftiTools graphics
library(rmarkdown)
library(gifti) #to read in your surface giftis 

command=sprintf("-volume-to-surface-mapping /Users/kahinim/Desktop/ITC/finalIms/4Dcopeseed1meanMNI.nii.gz /Users/kahinim/Desktop/fsaverage5_lh.pial_avg.surf.gii /Users/kahinim/Desktop/mean_l.shape.gii -trilinear")
ciftiTools::run_wb_cmd(command, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

command=sprintf("-volume-to-surface-mapping /Users/kahinim/Desktop/ITC/finalIms/4Dcopeseed1meanMNI.nii.gz /Users/kahinim/Desktop/fsaverage5_rh.pial_avg.surf.gii /Users/kahinim/Desktop/mean_r.shape.gii -trilinear")
ciftiTools::run_wb_cmd(command, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

test_l <- read_gifti("/Users/kahinim/Desktop/mean_l.shape.gii")
test_r <- read_gifti("/Users/kahinim/Desktop/mean_r.shape.gii")
test.ciftimap <- as_cifti(cortexL = test_l$data$normal, cortexR = test_r$data$normal) 
write_cifti(test.ciftimap, "/Users/kahinim/Desktop/mean.dscalar.nii")

surfL_name <- read_surf("/Users/kahinim/Desktop/fsaverage5_lh.inflated_avg.surf.gii")
surfR_name <- read_surf("/Users/kahinim/Desktop/fsaverage5_rh.inflated_avg.surf.gii")

test.ciftimap <- add_surf(test.ciftimap, surfL=surfL_name, surfR=surfR_name)
view_cifti(test.ciftimap, widget=TRUE, colors=c("blue","gray","red"), zlim=c(-0.5,0.5))

command=sprintf("-volume-to-surface-mapping /Users/kahinim/Desktop/ITC/finalIms/zstat2MNI.nii.gz /Users/kahinim/Desktop/fsaverage5_lh.pial_avg.surf.gii /Users/kahinim/Desktop/logk_l.shape.gii -trilinear")
ciftiTools::run_wb_cmd(command, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

command=sprintf("-volume-to-surface-mapping /Users/kahinim/Desktop/ITC/finalIms/zstat2MNI.nii.gz /Users/kahinim/Desktop/fsaverage5_rh.pial_avg.surf.gii /Users/kahinim/Desktop/logk_r.shape.gii -trilinear")
ciftiTools::run_wb_cmd(command, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

test_l <- read_gifti("/Users/kahinim/Desktop/logk_l.shape.gii")
test_r <- read_gifti("/Users/kahinim/Desktop/logk_r.shape.gii")
test.ciftimap <- as_cifti(cortexL = test_l$data$normal, cortexR = test_r$data$normal) 
write_cifti(test.ciftimap, "/Users/kahinim/Desktop/logk.dscalar.nii")

surfL_name <- read_surf("/Users/kahinim/Desktop/fsaverage5_lh.inflated_avg.surf.gii")
surfR_name <- read_surf("/Users/kahinim/Desktop/fsaverage5_rh.inflated_avg.surf.gii")

test.ciftimap <- add_surf(test.ciftimap, surfL=surfL_name, surfR=surfR_name)
view_cifti(test.ciftimap, widget=TRUE, colors=c("blue","gray"), zlim=c(-5,-3))
view_cifti(test.ciftimap, widget=TRUE, colors=c("gray","red"), zlim=c(3,5))


```

Surface files used for this projection are available [here](https://github.com/PennLINC/PennLINC.github.io/tree/main/docs/Tutorials/fsaverage5).

Spin testing was performed on the data column in the gifti, which was read in using the library 'gifti' in R: 
```
library(gifti)
col_file = '/Users/kahinim/Desktop/logk_r.shape.gii'
data_r = readgii(col_file)$data
write.table(data_r, file = "/Users/kahinim/Desktop/R_logk_data.txt", sep = "", row.names = FALSE)

col_file = '/Users/kahinim/Desktop/logk_l.shape.gii'
data_r = readgii(col_file)$data
write.table(data_r, file = "/Users/kahinim/Desktop/L_logk_data.txt", sep = "", row.names = FALSE)

col_file = '/Users/kahinim/Desktop/mean_r.shape.gii'
data_r = readgii(col_file)$data
write.table(data_r, file = "/Users/kahinim/Desktop/R_mean_data.txt", sep = "", row.names = FALSE)

col_file = '/Users/kahinim/Desktop/mean_l.shape.gii'
data_r = readgii(col_file)$data
write.table(data_r, file = "/Users/kahinim/Desktop/L_mean_data.txt", sep = "", row.names = FALSE)

```

Then, using the code detailed in the `spin_test` folder in this project's repo, Dr. Erica Baller performed the actual spin testing.  

### 8. Code added during revision

1. The histogram was created using R using `hist`
----------------------------------------------------------------------------------------------------------------------------------------------------------
