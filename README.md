# RNA-NRD 

## RNA-NRD: a Non-redundant RNA Structural Dataset for Benchmarking and Functional Analys  

### Install Instructions 
RNA-NRD source code is implemented using Python 3.8.10 and can be executed in 64-bit Linux Machine. It uses the tool STAR3D for 3D structure alignment which is provided here. STAR3D is implemented by using java 1.7 and requires JRE to run.

#### Install JRE:  
```
Debian/Ubuntu: apt install default-jre
Fedora/CentOS: dnf install default-jre 
```
#### Install python3:
```
Debian/Ubuntu: apt install python3.8  
Fedora/CentOS: dnf install python3.8 
```
#### Install pip3: 
```
Debian/Ubuntu: apt install python3-pip  
Fedora/CentOS: dnf install python3-pip  
```
#### Install required Python libraries:  
It is required to install several python libraries to run RNA-NRD pipeline. These libraries are included in the [requirements.txt](requirements.txt) file. To install all required python libraries, please navigate to the RNA-NRD home directory in the terminal and execute the following command.

```
pip install -r requirements.txt
``` 

#### Python packages that should already exist:  
os, sys, shutil, math, random, subprocess, glob, time, argparse, logging, requests  
  
*** If any of the above mentioed package doesn't exist, then please install with command 'pip3 install package-name' ***



### Collect RNA chains from PDB
*** Steps:
1. Download RNA chains from PDB in csv format. Go to homepage button -> “^[0-9] Nucleic Acid Containing Structures” -> “Tabular Report” -> “Create CustomReport” -> select attributes -> “Run Report”.
2. Select these 12 attributes: Experimental Method, Release Date, PDB ID, Number of Distinct RNA Entities, Resolution (Å), Sequence, Entity Polymer Type, Polymer Entity Sequence Length, Source Organism, Macromolecule Name, Chain ID (Asym ID), Entry Id (Polymer Entity Identifier).
3. Convertcsv files to txt file format and merge them. See the example code below and change the file names accordingly:
		
		****************************************************************************************
		csvformat -T rcsb_pdb_custom_report_0001-2500.csv > rcsb_pdb_custom_report_0001-2500.txt
		sed -i 's/\ /_/g' rcsb_pdb_custom_report_0001-2500.txt

		csvformat -T rcsb_pdb_custom_report_2501-5000.csv > rcsb_pdb_custom_report_2501-5000.txt
		sed -i 's/\ /_/g' rcsb_pdb_custom_report_2501-5000.txt

		csvformat -T rcsb_pdb_custom_report_5001-5329.csv > rcsb_pdb_custom_report_5001-5329.txt
		sed -i 's/\ /_/g' rcsb_pdb_custom_report_5001-5329.txt

		cat *.txt > Data.txt
		****************************************************************************************
		
*** Note: 'Merged_data.txt' already exists inside folder 'Data' which contains all the RNA chains collected from PDB on 03/17/21.



### Collect RNA chain family information from Rfam website   
*** Requirement: Install chromedrive and update the chromedrive path in the code Rfam_parser.py on line 26   
*** Instructions: Run the code inside the folder 'Rfam_parser.py'. If chromedriver not found, then update the 'DRIVER_PATH' in the code.     
*** Run command: python3 Rfam_parser.py  
*** Output: generates file 'PDB_family_name.txt' contaiing Rfam family information insdie the folder 'Rfam_family'.   
*** Note: 'PDB_family_name.txt' already exists inside folder 'Rfam_family' which contains all the Rfam family information collected from Rfam website on 01/23/23.     


### Run Instructions
  
#### Generate Non-redundant Dataset  
  
*** Run command: python3 Run.py [-o 'Output/Nonredundat_datalist'] [-r 4] [-a 80] [-fr 50]  
*** Help command: python3 Run5_NR_Dataset.py -h  
*** Optional arguments:  
  -h, --help  show this help message and exit  
  -i [I]      Input file name containing RNA chains. Default: 'Data.txt'.  
  -o1 [O1]    Output file name with organism division. Default: 'RNA-NRD_Dataset'.  
  -o2 [O2]    Output file name without organism division. Default: 'RNA-NRD_without_Organism_Division_Dataset'.  
  -t [T]      Provide sequence identity threshold. Default: 80.  
  -r [R]      Provide RMSD threshold. Default: 4.  
  -a [A]      Provide structural alignment ratio threshold. Default: 80.  
  -org [ORG]  Generate RNA-NRD dataset without organism based division. Default: True.  
*** Output: Generates final nonredundant detaset output file inside user defined location (Default: 'Output/Nonredundat_datalist')  

            
### Important Notes
*** Please make sure the file 'All_Organism_list' is inside the folder 'Organism_list'. It contains list of all organisms currently present in PDB. In case the input file has any new organisms, then the new organisms will be added in the file 'All_Organism_list'.                


### Terms  
Where appropriate, please cite the following RNA-NRD paper:  
Nabila et al. "RNA-NRD: a Non-redundant RNA Structural Dataset for Benchmarking and Functional Analys."  
