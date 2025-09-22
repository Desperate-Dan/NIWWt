# NIWWt
Setting up the pipeline for SARS-CoV-2 WW analysis of NI samples.

## How to run
Currently works assuming you already have [Nextflow](https://www.nextflow.io/docs/latest/install.html) installed. Next create the environment using the following command:
```
git clone https://github.com/Desperate-Dan/NIWWt.git
cd NIWWt
conda env create -f ./environment.yml
```
The code itself assumes that your conda environments are located at ```${HOME}/miniconda3/envs/NIWWt```. If this is not the case you'll need to change that bit of the code in [main.nf](https://github.com/Desperate-Dan/NIWWt/blob/main/main.nf).

To run the tool itself you can run the below command:
```
nextflow -profile conda /path/to/NIWWt/main.nf --fastq /path/to/fastq/directory/
```
The pipeline expects a directory of pairs of fastq files. There are default reference and primer bed files included in ```resources``` but you can use your own with the ```--ref /path/to/ref``` and ```--bedfile /path/to/bedfile``` flags.
