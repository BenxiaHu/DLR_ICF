# DLR_ICF  
This tutorial will introduce how to run DLR_ICF to analyze DLR and ICF.

### DLR_ICF can be used to analyze DLR and ICF based on normlazed Hi-C contact matrix.  
DLR: distal-to-local ratio  
ICF: inter-chromosomal fraction  

#### calculate the DIR and ICF for individual Hi-C contact matrix
#### usage: 
```DLR_ICF_main [-h] [-b] -I INPUTPATH -f FILENAME -d DISTANCE -r RESOLUTION -O OUTPATH -c CHRSIZE -o OUTFILE``` 

the format of input file is cool format.  
the output file is bedgraph which can be visualized in IGV/UCSC or other genome browsers. 

example:  
1. used the iced normalized contact matrix to calculate DLR/ICF.  
```DLR_ICF_main **-b** -I "the path of matrix" -f "test" -d 1000000 -r 10000 -o 1Mb_test_10kb -O "the path of output" -c hg38.chrom.sizes```

2. used the balanced contact matrix (mcool) to calculate DLR/ICF.  
```DLR_ICF_main -I "the path of matrix" -f "test" -d 1000000 -r 10000 -o 1Mb_test_10kb -O "the path of output" -c hg38.chrom.sizes```

                     
optional arguments:  
|  |   |    |   |   |
|:----:|:-----:|:----:|:------:|:------:|  
| -h |  |--help|| show this help message and exit |
| -b ||  --balanced |   | contact matrix is iced or balanced |
| -I | INPUTPATH  | --inputpath | INPUTPATH |path of input file  |  
| -f | FILENAME   | --filename    | FILENAME |name of input file |
| -d | DISTANCE  | --distance |DISTANCE|the distance of distal chromation interactions|
| -r    |   RESOLUTION| --resolution | RESOLUTION| resolution of contact matrix  | 
| -O | OUTPATH    | --outpath |  OUTPATH |path of output file  |  
| -c | CHRSIZE    | --chrsize |  CHRSIZE |chromosome size file  |
| -o | OUTFILE    | --outfile |  OUTFILE |name of output file  |


#### divide the DIR and ICF into different groups based on compartment definitions.  
```DLR_ICF_separation [-h] [-b] -I INPUTPATH -f FILENAME -d DISTANCE -p Compartment -r RESOLUTION -O OUTPATH -c CHRSIZE -o OUTFILE``` 

optional arguments:  
|  |   |    |   |   |
|:----:|:-----:|:----:|:------:|:------:|  
| -h |  |--help|| show this help message and exit |
| -b ||  --balanced |   | contact matrix is iced or balanced |
| -I | INPUTPATH  | --inputpath | INPUTPATH |path of input file  |  
| -f | FILENAME   | --filename    | FILENAME |name of input file |
| -d | DISTANCE  | --distance |DISTANCE|the distance of distal chromation interactions|
| -p | compartment  | --compartment |compartment|compartment of chromaitn|
| -r    |   RESOLUTION| --resolution | RESOLUTION| resolution of contact matrix  | 
| -O | OUTPATH    | --outpath |  OUTPATH |path of output file  |  
| -c | CHRSIZE    | --chrsize |  CHRSIZE |chromosome size file  |
| -o | OUTFILE    | --outfile |  OUTFILE |name of output file  |

#### output ICF groups.  
```ICF_chromatin [-h] [-b] -I INPUTPATH -f FILENAME -p Compartment -r RESOLUTION -O OUTPATH -c CHRSIZE -o OUTFILE``` 

optional arguments:  
|  |   |    |   |   |
|:----:|:-----:|:----:|:------:|:------:|  
| -h |  |--help|| show this help message and exit |
| -b ||  --balanced |   | contact matrix is iced or balanced |
| -I | INPUTPATH  | --inputpath | INPUTPATH |path of input file  |  
| -f | FILENAME   | --filename    | FILENAME |name of input file |
| -p | compartment  | --compartment |compartment|compartment of chromaitn|
| -r    |   RESOLUTION| --resolution | RESOLUTION| resolution of contact matrix  | 
| -O | OUTPATH    | --outpath |  OUTPATH |path of output file  |  
| -c | CHRSIZE    | --chrsize |  CHRSIZE |chromosome size file  |
| -o | OUTFILE    | --outfile |  OUTFILE |name of output file  |

#### compare the difference in the DIR and ICF for 2 Hi-C contact matrices.  
#### usage: 
```DLR_ICF_comparison [-h] -i INPUTPATH -t TREATMENT -c CONTROL -r RESOLUTION -O OUTPATH -o OUTFILE```  

example:  
```DLR_ICF_comparison -i "the path of matrix" t "caseid" -c "controlid" -r 10000 -o 1Mb_test_10kb -O "the path of output"```  

Here we convert the DLR/ICF ratio into z-score, and calculate p value and fdr.  

the output file is bedgraph which can be visualized in IGV/UCSC or other genome browsers. 

optional arguments:  
|  |   |    |   |   |
|:----:|:-----:|:----:|:------:|:------:|  
|  -h    |      |   --help   |           | show this help message and exit         |
|  -i    |   INPUTPATH | --inputpath  | INPUTPATH | path of input file             |
|  -t    |   TREATMENT | --treatment  | TREATMENT | name of treatment file         |
|  -c    |   CONTROL   | --control    | CONTROL   | name of control file             |
|  -r    |   RESOLUTION| --resolution | RESOLUTION| resolution of contact matrix  | 
|  -O    |   OUTPATH   | --outpath    | OUTPATH   | path of output file              |
|  -o   |    OUTFILE   | --outfile    | OUTFILE   | name of output file              |

### Installation 
#### requirement for installation
python>=3.8  
numpy  
pandas  
argparse  
cooler   
h5py  
scipy.stats   
statsmodels.stats.multitest  

#### pip install DLR-ICF==1.0.4
https://pypi.org/project/DLR-ICF/1.0.4/

#### reference
Transcription Elongation Can Affect Genome 3D Structure (https://doi.org/10.1016/j.cell.2018.07.047)
