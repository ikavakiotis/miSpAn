# miSpAn


## Running miSpAn
In order to run miSpAn download the miSpAn folder and run the mispan.py. The input parameters and paths should be declared in cfg.txt file which can be found in cfg folder


## Congifuaration file - cfg.txt
The following parameters should be set:
- **mirConfidence_file** This file contains the information about the miRNA confidence from miRBase_V22. miSpAn provides the current annotation file, so users should change it only if they would like to provide a file with updated annotations.  e.g. _mirConfidence_file=mature_precursors_confidence_miRBase22_human_DB.txt_  
- **expression_repository_path:** The path to the folder containing all expression data. User should prvide the **exact path** to the folder, e.g: 
  - Windows: expression_repository_path=C:\data\expression\samples\
  - linux: expression_repository_path=C:\data\expression\samples\
- **metaData_file:** The path to the metadata file. User should provide the **exact path** to the file, e.g: 
  - Windows: expression_repository_path=C:\data\expression\metadata.txt
  - linux: expression_repository_path=C:\data\expression\metadata.txt
- **mirna_confidence:** Available options: high, low, all (only one). The analysis will only proceed with the selected miRNAs depending on the annotation from mirConfidence_file. 
- **Conditions:** The conditions that will be considered in the study. The conditions should be separated with (#) and should be written exatly as they appear in metadata file (e.g. Conditions=skin#brain)
- **status:** Status comparison. The status values should be separated with (#) and should be written exatly as they appear in metadata file (e.g. status=skin#brain) If no status comparison leave it blank (e.g. status= ).
- **tsi_threshold:**. The TSI threshold for the presented results (e.g tsi_threshold=0.98)


## Results - Output files

All generated output files belong to one of the following categories (also stated in file's name first letters): Datasets (Ds) in .tsv format, Reports (Rp), Results (Rs) and Figures (Fig) in png format.  

Below the detailed structure of the output folder is presented. Each time an new run/analysis is preformed, a new folter (named with current date and time) inside the output folder is created and all files are generated in this folder. In case users are not utilizing "statuses" only "all" folders are created. If users utilize statuses then new folders named after the statuses are created to containg the corresponding results. For simplicity reasons only the structure of "all" folder is depicted.

- Folder: _ouput_/_current datetime/_
  - File: _Rp_logfile.txt_ - **Report** file containing input parameters for tracking reasons
  - Folder: _1_msp_data_integration_
    - File: _Ds_samplesExpression_rc.txt_ - A **dataset** containing the expression values (read counts) of each miRNA (rows) across all samples (columns) and it is offered for conducting differential expression analysis with dedicated software, such as DESeq2 or limma.  
    - File: _Ds_samplesExpression_rpm.txt_ - A **dataset** in which samples serve as rows, whereas miRNAs as features/columns, containing the expression of miRNAs in reads per million (RPM). This dataset is utilized for the downstream analysis. 
    
  - Folder: _2_msp_condition_expression_
    - File: _Rp_samplesPerCondition.txt_ - **Report** file containing the number of samples per condition
    - Folder: _all_
      - File: _Ds_all_ConditionExpression_log2rpm.txt_ A **dataset** in which conditions serve as rows, whereas miRNAs as features/columns, containing the expression of miRNAs in lor2(RPM). 
      - File: _Ds_all_ConditionExpression_rpm.txt_ - A **dataset** in which conditions serve as rows, whereas miRNAs as features/columns, containing the expression of miRNAs in reads per million (RPM).
      - File: _Fig_all_Condition_Clustermap.png_ - A **figure** (clustermap) depicting the expression of all miRNAs in different conditions. 
    - Folder: status_1
      - _The same structure as Folder: all_ - Contains results for the analysis conducted utilizing only samples belonging to status_1
    - Folder: status_2
      - _The same structure as Folder: all_ - Contains results for the analysis conducted utilizing only samples belonging to status_2
 - Folder: _3_msp_condition_specificity_
    - File: _Rp_samplesPerCondition.txt_ - **Report** file containing the number of samples per condition
    - Folder: _all_
      - File: _Rs_all_specificMiRNAs.txt_ - **Result** file in tsv format containing the most specific miRNAs sorted by their TSI value. Columns include: (A)TSI score, (B)miRNA, (C)The Condition in which miRNA has the highest expression (to which is specific), and (D) Expression in log2RPM
      - File: _Rs_all_tauTSI.txt_ - **Result** file containing the sorted miRNAs according to their TSI value.
      - File: _Rs_all_spMiRNAs_PerCondition.txt_ - **Result** file in tsv format containing the specific miRNAs in each condition
      - File: _Rs_all_spMiRNAs_ConditionExpression.txt_ **Result** file in tsv format. In this file each column is dedicated to one specific miRNA. The first row contains the tsi value and the miRNA name (e.g. 0.998_hsa-miR-449b-3p). The rest columns contain the conditions along with their expression of the specific miRNAs sorted from the highest to the lowest.
      - File: _Fig_all_Condition_Clustermap.png_ - A **figure** (clustermap) depicting the expression of all miRNAs in different conditions. 
      - File: A **dataset** containing the data from which the _Fig_all_Condition_Clustermap.png_ has occured, offering the opportunity to the user to produce their own figures. 
    - Folder: status_1
      - _The same structure as Folder: all_ - Contains results for the analysis conducted utilizing only samples belonging to status_1
    - Folder: status_2
      - _The same structure as Folder: all_ - Contains results for the analysis conducted utilizing only samples belonging to status_2
    
  


