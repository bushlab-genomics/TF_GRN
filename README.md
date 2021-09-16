# TF GRN Project
The code and files in this repo correspond to the paper titled **Modeling transcriptional regulation using gene regulatory networks based on multi-omics data sources**. We have provided test files here for creating gene regulatory networks (GRNs) and building ElasticNet(ENET) prediction models. 

## Requirements(software)
python3\
R(>=3.4.2)\
bedtools(https://github.com/arq5x/bedtools2/releases)

### python packages
gseapy(0.9.13)\
numpy(1.17.2)\
pandas(0.25.1)\
scipy(1.4.1)\
scikit-learn(0.21.2)

### R packages
pandaR(1.18.0)\
biomaRt(2.42.1)\
dplyr(0.8.5)

### Tutorial
The main script for running the code is "run_tf_grn.py", which was written to create PANDA GRNs. This script mainly requires a file containing list of genes and a TF peak file for inputs. The ouput is a weighted PANDA GRN built using different regulatory regions for each gene.\

Following options are mainly required to run the script\
* -i, --inputs: The inputs required for the script in the following sequence (gene list, TF bed peak file, CTCF bed peak file, PANDA PPI file, PANDA expression file. User can write "none" instead of the name of CTCF peak file if the user wants to define the cis-regulatory regions of genes without using CTCF peaks. Additionally, user can replace PANDA PPI or expression files  if they want to create PANDA GRN without these datasets. 
* -w, --window: The bp distance that the user needs to mention to define the cis-regulatory regions in order to find TFBS.(default=50000bp)
* -m, --multiprocessin: Select this flag in order to run the script as a parallel process. 
* -n, --cpu : Use this option along with the -m flag to specify the number of cores/cpus to be used for parallel processing.
* -d, --hammingdistance: This option specifies the hamming distance to be used for the convergence of the PANDA GRN(default = 0.001). 
* -o, --output: Specify the output prefix for various output files as well as for the output directory, where they will be stored. 
* -t, --type: This option is used to specify the kind of analysis that will be done to produce the GRNs. The following options are available:
    * "grn" : This option creates a GRN using cis-regulatory TFBS. 
    * "promoter": This option creates a GRN based on promoter TFBS and requires the following additional options to be used:
        * -p, --promoterwindow: The bp distance upstream of the TSS of the gene to be used to define promoter TFBS (default = 5000bp). 
    * "distal": This option creates a GRN based on distal TFBS and requires the -p option to be used.  
    * "intron" : This option creates a GRN based on intronic TFBS. 
    * "hic" : This option creates a Hi-C contacts normalized GRN and requires the -p option to be used in addition to the following options:
        * -c, --contacts: Specify the file containing normalized Hi-C contacts.
        * -u, --upweighting: This flag is to be used if the user wants to upweight the promoter TFBS not containing any Hi-C contacts. 
        * -s, --scalingfactor: This option specifies the scaling factor in the range [0,1] to scale the Hi-C normalized TF-TG interactions(default = 1.0). 
        
  ### Use-cases
  The following command will create a PANDA GRN with TFBS defined using a 50Kbp cis-regulatory window defined by CTCF peak boundaries using parallel processing and hamming distance threshold of 0.001. It will store the results in the "toy_grn" directory.\
  ```
  python run_tf_grn.py -i gene_list.txt toy_peaks.bed ctcf_peaks.bed toy_panda_ppi.csv toy_panda_expression.csv -t grn -w -m -n 20 -d -o toy_grn
  ```
  
  The following command will create a PANDA GRN with TFBS defined using a 50Kbp cis-regulatory window without the CTCF boundaries and the PPI data. It also won't use parallel processing.\
  ```
  python run_tf_grn.py -i gene_list.txt toy_peaks.bed none none toy_panda_expression.csv -t grn -w -d -o toy_grn
  ```
  The following command will create a PANDA GRN using intronic TFBS.\
  ```
  python run_tf_grn.py -i gene_list.txt toy_peaks.bed ctcf_peaks.bed toy_panda_ppi.csv toy_panda_expression.csv -t intron -w -m -n 20 -d -o toy_intron
  ```
  The following commands will create a PANDA GRN using promoter and distal TFBS respectively based on promoter window of 5000bp:
  ```
  python run_tf_grn.py -i gene_list.txt toy_peaks.bed ctcf_peaks.bed toy_panda_ppi.csv toy_panda_expression.csv -t promoter -w -m -n 20 -p -d -o toy_promoter
  python run_tf_grn.py -i gene_list.txt toy_peaks.bed ctcf_peaks.bed toy_panda_ppi.csv toy_panda_expression.csv -t distal -w -m -n 20 -p -d -o toy_distal
  ```
  The following command will create a Hi-C normalized PANDA GRN using a scaling factor of 1.0 as well as by upweighting the promoter TFBS:
  ```
 python run_tf_grn.py -i gene_list.txt toy_peaks.bed ctcf_peaks.bed toy_panda_ppi.csv toy_panda_expression.csv -t hic -c toy_hic_contacts.csv -u -w -m -n 20 -p -s -d -o toy_hic
  ```
  The other script that is of use is "Code/get_pred_eval.py". The output of this script is a dataframe , produced in the "pred_results" directory, containing Pearson's Correlation Coefficient and MSE produced after running the ENET prediction models. This script has the following options:\n
  ```
  * -i, --input: The name of the input file containing the matrices/GRNs to be used for features in order to build the ENET prediction models. 
  * -e, --expression: The name of the expression file containing the outcome FPKM values to be used for training the testing the ENET models. Place this file in the Data/ folder. 
  * -o, --output: The output prefix for the prediction results (PCC and MSE) that will be produced from the trained ENET models within the "pred_results" directoy. 
  * -r, --random: The number of random iterations for which to train and test the ENET models(default = 20). 
  * -m, --multiprocess: Select this flag if the users wants to run the script in a parallel manner. 
  * -n, --cpu: The number of cores/cpu to be specified along with the -m flag. 
  ```
  The following set of scripts will build the prediction models using the normal and the Hi-C normalized GRNs for 20 iterations for a common set of genes.\
  ```
  python run_tf_grn.py -i gene_list.txt toy_peaks.bed ctcf_peaks.bed toy_panda_ppi.csv toy_panda_expression.csv -t grn -w -m -n 20 -d -o toy_grn**\
  python run_tf_grn.py -i gene_list.txt toy_peaks.bed ctcf_peaks.bed toy_panda_ppi.csv toy_panda_expression.csv -t hic -c toy_hic_contacts.csv -u -w -m -n 20 -p -s -d -o toy_hic\
  echo "toy_hic_results/toy_hic_reg_net.tsv" >> toy_inputs.txt
  echo "toy_grn_results/toy_grn_reg_net.tsv" >> toy_inputs.txt
  python Code/get_pred_eval.py -i toy_inputs.txt -e toy_pred_expression.csv -m -n 20 -r -o toy
  ```
  
  
  
  
