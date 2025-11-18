# NIWWt
Setting up the pipeline for SARS-CoV-2 WW analysis of NI samples. Broadly it uses [TrimGalore](https://github.com/FelixKrueger/TrimGalore) for read trimming, [bwa mem](https://github.com/lh3/bwa) for mapping, [iVar](https://github.com/andersen-lab/ivar) for primer-trimming and [freyja](https://github.com/andersen-lab/Freyja) for SARS_CoV-2 variant demultiplexing.

## How to run
Currently works assuming you already have [Nextflow](https://www.nextflow.io/docs/latest/install.html) installed. Next create the environment using the following command:
```
git clone https://github.com/Desperate-Dan/NIWWt.git
cd NIWWt
conda env create -f ./environment.yml
```
The code itself assumes that your conda environments are located at ```${HOME}/miniconda3/envs/NIWWt```. If this is not the case you'll need to change those bits of the code in [main_pipeline.nf](https://github.com/Desperate-Dan/NIWWt/blob/main/modules/main_pipeline.nf).

To run the tool itself you can run the below command:
```
nextflow run -profile conda /path/to/NIWWt/main.nf --fastq /path/to/fastq/directory/
```
The pipeline expects a directory of pairs of fastq files and will produce a ```results``` directory in the directory where you launched your command. There are default reference and primer bed files included in ```resources``` but you can use your own with the ```--ref /path/to/ref``` and ```--bedfile /path/to/bedfile``` flags.
