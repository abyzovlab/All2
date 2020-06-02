# All2
A tool for filtering variants from all2all comparison of multiple clones or single cells

## **Prerequisite:**
1. Python 3.6 with the following packages
	1. Pandas
	1. Matplotlib
	1. Seaborn
	1. Numpy
1. Other dependencies:
    1. Samtools
    
## **Setup:** 
### Download
   ```
      git clone https://github.com/abyzovlab/Scellector.git    
   ```
### Configuration setup
1. Samtools (download and install):
     ```
     wget https://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2
    
     tar -xvf samtools-1.9.tar.bz2
    
    cd samtools-1.9/
    
    /.configure
    
    make
     ```
1. Reference (download and update "REFERENCE"  in config_file/config.txt):

    ```
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz
    
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.fai.gz
    
    gunzip human_g1k_v37_decoy.fasta.gz
    
    gunzip human_g1k_v37_decoy.fasta.fai.gz
   ```