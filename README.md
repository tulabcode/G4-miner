# G4-miner

G4-miner is a genome-wide approach to directly identify G4 structures from standard sequencing data. G4-miner inspects sequencing quality in the whole genome, seizes unexpected fluctuations in local, and ascribes some of the quality fluctuations to the unstable formation of G4.



## Installation
Using Git:
```shell
git clone https://github.com/tulabcode/G4-miner
```



## Prerequisities
G4-miner is implemented in perl and shell, a command line software.

Required:
1. Linux system
2. [g4predict](https://github.com/mparker2/g4predict), install it following the installtion before run G4-miner.
3. perl,  recommended version 5.18.2 or higher
4. It is recommended to install in advance three perl moudles: Getopt::Long, File::Basename, File::Spec.
5. python2,  recommended version 2.7.16


## Running G4-miner for predictting
Follow the instructions:
```shell
# enter the G4-miner directory
cd path_to_G4-miner_directory

# add rights to execute
chmod +x G4-miner.pl

# run G4-miner
./G4-Miner.pl or perl /path_to_G4-miner_directory/G4-miner.pl
Usage:
        perl G4-miner.pl [options]
Options:
                -i <file> <input bam/sam file>
                -d <path> <the root directory for saving output>
                -g <file> <the fasta file of the whole genome reference>
                -rm <int> <rm the directory of $hdir/$samp, default: 1, yes>
                -ch <string> <specify some specific chromosomes for analysing, 
                             -ch chr1 -ch chr2 ... -ch chrX> <Optional, Default genome wide>
                -py <file> <the absolute path of python2> (default: `which python2`)
                -g4 <file> <the absolute path of g4predict> (default: `which g4predict`)
                -rl <int> <the longest length of reads, default = 150 nt>
                -dl <int> <the length of detection region, defalut = 75 nt>
                -sp <float> <Specie for selecting cutoff of G ratio, here have six species:
                             human/mouse/fruitfly/thaliana/elegans/yeast,
                             default: human with 0.28, if the specie did not include in this
                             paper, can calculate by one more step>
                -n <string> <the name of the sample>
```

The default model is set to a model trained on the human dataset, as described in our paper.

**NOTE**:
The following table shows the corresponding G-ratio of different species:

|           **Species**            | **G-ratio** |
| :------------------------------: | :---------: |
|       **Homo  sapiens**      |    0.28     |
|       **Mus  musculus**       |    0.28     |
| **Drosophila  melanogaster**  |    0.27     |
|   **Arabidopsis  Thaliana**   |    0.24     |
|  **Caenorhabditis  elegans**  |    0.25     |
| **Saccharomyces  cerevisiae** |    0.29     |


## Reference genomes and data sources
|           **Species**           | **Reference genome source**                                  | **Data  source (SRR code)** |
| :-----------------------------: | :----------------------------------------------------------- | :-------------------------: |
|       **Homo sapiens**       | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/             |         SRR9644818          |
|       **Mus musculus**       | https://hgdownload.soe.ucsc.edu/goldenPath/mm10/             |         SRR13179566         |
| **Drosophila melanogaster**  | https://hgdownload.soe.ucsc.edu/goldenPath/dm6/              |         SRR12822760         |
|   **Arabidopsis Thaliana**   | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/ |         SRR11608990         |
|  **Caenorhabditis elegans**  | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/ |         SRR8816429          |
| **Saccharomyces cerevisiae** | https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/          |         SRR13747318         |


## sacCer Test data
### Download
Here, we provide a dataset from Saccharomyces cerevisiae for testing G4-miner. Contains the reference genome file of Saccharomyces cerevisiae (sacCer3.fa) and indexed (sacCer3.fa.fai), and there is also a bam file that is aligned to the genome, sorted (SRR13747318.sort.bam) and indexed (SRR13747318.sort.bam.bai).
Could download from [Southeast University Document Cloud](https://pan.seu.edu.cn:443/link/970CF07A79536277306E76D9A1BFD1BA).
Download test datasets, create a test directory, upload the downloaded files by using `rz` or `FileZilla`, and then run test.

### run Test
Just following the instruction for running:
```shell
# Download test datasets

# create a test directory
cd ~ && mkdir -p work/sacCer && cd ~/work/sacCer

# upload the downloaded files by using `rz` or `FileZilla`

# run test
perl /path_to_G4-miner_directory/G4-miner.pl -i ~/work/sacCer/SRR13747318.sort.bam -d ~/work/sacCer/G4-miner -n sacCer -g ~/work/sacCer/sacCer3.fa -py /usr/bin/python2 -g4 /usr/local/bin/g4predict -sp yeast
## here "~"" represents the personal home directory of Linux
```
### Result analysis
This project provides a case analysis of the species **Saccharomyces cerevisiae** (*Yeast*), and the output results are shown in the folder "~/work/sacCer/G4-miner".

In the output directory, each chromosome generated a corresponding experimental result folder. According to the experimental results of each chromosome, this tool counted the number of all positive regions, negative regions and false-positive regions, and calculated the number of regions that meet the conditions of different M and N parameters and found the parameter values corresponding to the most representative low-quality regions. In addition, this tool also distinguished between forward and reverse strand. The suffix with ‘f’ corresponds to the forward strand, and the suffix with ‘r’ corresponds to the reverse strand.

In each chromosome directory, the forward and reverse strand respectively generated 11 result files, and each file contained the following content:

**<Median.gz>**: The median of the quality value corresponding to each base on the chromosome

**<PQ_g4predict>**: The PQ file on the whole genome predicted by G4predict, named "PQ_g4predict_minus" on the forward strand, and "PQ_g4predict_plus" on the reverse strand

**<positive_array>**: The number of positive regions corresponding to different M and N parameter conditions

**<negative_sloci.gz>**: Position information of the negative region corresponding to different M and N parameter conditions

**<false_positive_sloci.gz>**: Position information of the false-positive region corresponding to different M and N parameter conditions

**<negative_array>**: The number of negative regions corresponding to different M and N parameter conditions.

**<false_positive_array>**: The number of false-positive regions corresponding to different M and N parameter conditions.

**<OQ.gz>**: Low-quality region sequence collection, that is, the observed quadruplex (OQ), there is overlap between the sequences

**<G4_OQ_unmerged>**: A collection of sequences related to the G-quadruplex with overlap screened by the characteristics of the sequence

**<G4_OQ_merged>**: A collection of sequences related to the G-quadruplex without overlap screened by the characteristics of the sequence

**<G4_PQinOQ>**: PQ set included by OQ


## Contact information
Jing Tu,  jtu@seu.edu.cn


## Citing G4-miner Article
Tu, J., Duan, M., Liu, W. et al. Direct genome-wide identification of G-quadruplex structures by whole-genome resequencing. Nat Commun 12, 6014 (2021). https://doi.org/10.1038/s41467-021-26312-w

