

import msp_data_integration as din
import msp_tissue_expression as tex

import msp_utils as utl
import msp_condition_datasets as cd
import msp_condition_specificity_metrics as csm
import msp_condition_specificity as cs
import msp_visualizations as viz
import msp_config_parser as cfgp


import datetime
import time
import sys
import os

def pipeline(working_dir):
    
    '''Case 1: If not considering status, run the analysis with all data'''
    
    '''Case 2: If considering status, 
        a) run the analysis as a whole to provide results with "extended" conditions e.g brain --> hlty_brain dis_brain
        b) run the analysis as different datasets, i.e healthy dataset and diseased dataset
    
    
    Case 1 and 2a are perfomred by the same code. If i_statusComparison is true the analysis continues
    
    '''
    
    tempFolder=time.strftime("%Y%m%d-%H%M%S")
    
    
    
    
    df_pathsVariables=cfgp.parse_config_file(working_dir)
    
    
    iPath_samples=df_pathsVariables['expression_repository_path']
    ifile_mirConfidence="lib"+os.sep+df_pathsVariables['mirConfidence_file']
    ifile_metaData=df_pathsVariables['metaData_file']
    ls_conditions=df_pathsVariables['conditions']
    i_confidence=df_pathsVariables['mirna_confidence']
    i_statusComparison=df_pathsVariables['status']
    i_thresholdToPrint=df_pathsVariables['tsi_threshold']
    
    opath_main=working_dir+"output"+os.sep+tempFolder+os.sep
    #print opath_main
    #print iPath_samples
    #sys.exit()
    #print ifile_mirConfidence
    #print ifile_metaData
    #print ls_conditions
    #print i_confidence
    #print i_statusComparison
    #print i_thresholdToPrint
    #sys.exit()
    
    
    if i_statusComparison:
        firstBinary=i_statusComparison[0]
        secondBinary=i_statusComparison[1]
    else:
        firstBinary=""
        secondBinary=""
    
    folder_1_di="1_msp_data_integration"
    folder_2_texpr="2_msp_condition_expression"
    folder_3_tspec="3_msp_condition_specificity"
    
    folder_all="all"
    folder_firstBinary=firstBinary
    folder_secondBinary=secondBinary
    
    opath_1_di=os.path.join(opath_main,folder_1_di)
    
    opath_2_texpr_all=os.path.join(opath_main,folder_2_texpr, folder_all)
    opath_2_texpr_firstBinary=os.path.join(opath_main,folder_2_texpr, folder_firstBinary)
    opath_2_texpr_secondBinary=os.path.join(opath_main,folder_2_texpr, folder_secondBinary)
    
    opath_2_texpr_bare=os.path.join(opath_main,folder_2_texpr)
    
        
    opath_3_tspec_all=os.path.join(opath_main,folder_3_tspec, folder_all)
    opath_3_tspec_firstBinary=os.path.join(opath_main,folder_3_tspec, folder_firstBinary)
    opath_3_tspec_secondBinary=os.path.join(opath_main,folder_3_tspec, folder_secondBinary)
    
    if not os.path.exists(opath_1_di):
        os.makedirs(opath_1_di)
    if not os.path.exists(opath_2_texpr_all):
        os.makedirs(opath_2_texpr_all)
    if not os.path.exists(opath_3_tspec_all):
        os.makedirs(opath_3_tspec_all)
    
    if i_statusComparison:
        if not os.path.exists(opath_2_texpr_firstBinary):
            os.makedirs(opath_2_texpr_firstBinary)
        if not os.path.exists(opath_2_texpr_secondBinary):
            os.makedirs(opath_2_texpr_secondBinary)
        if not os.path.exists(opath_3_tspec_firstBinary):
            os.makedirs(opath_3_tspec_firstBinary)
        if not os.path.exists(opath_3_tspec_secondBinary):
            os.makedirs(opath_3_tspec_secondBinary)
    
    
    #sys.exit()
    
    '''OUTPUT FILES'''
    
    '''Output folder 1 (Table 1 - 2 files): 1_msp_data_integration'''
    
    '''logfile: Parameters'''
    ofile_logfile=os.path.join(opath_main,"Rp_logfile.txt")
    
    '''T1.1 - 1 file'''
    ofile_rc_dataset=os.path.join(opath_1_di,"Ds_samplesExpression_rc.txt")
    '''T1.2 - 1 file'''
    ofile_rpm_dataset=os.path.join(opath_1_di,"Ds_samplesExpression_rpm.txt")
    
    
    '''Output folder component 2 (Table 2 - max 8 files): 2_msp_condition_expression'''
    
    '''T2.1 - Max 6 files'''
    ofile_allConditionDataset_rpm=os.path.join(opath_2_texpr_all,"Ds_all_ConditionExpression_rpm.txt")
    ofile_firstBinaryConditionDataset_rpm=os.path.join(opath_2_texpr_firstBinary,"Ds_"+firstBinary+"_ConditionExpression_rpm.txt")
    ofile_secondBinaryConditionDataset_rpm=os.path.join(opath_2_texpr_secondBinary,"Ds_"+secondBinary+"_ConditionExpression_rpm.txt")
  
    ofile_allConditionDataset_log2rpm=os.path.join(opath_2_texpr_all,"Ds_all_ConditionExpression_log2rpm.txt")
    ofile_firstBinaryConditionDataset_log2rpm=os.path.join(opath_2_texpr_firstBinary,"Ds_"+firstBinary+"_ConditionExpression_log2rpm.txt")
    ofile_secondBinaryConditionDataset_log2rpm=os.path.join(opath_2_texpr_secondBinary,"Ds_"+secondBinary+"_ConditionExpression_log2rpm.txt")
    
    '''T2.2 - 3 file'''
    ofile_all_clustermap=os.path.join(opath_2_texpr_all,"Fig_all_Condition_Clustermap.png")
    ofile_firstBinary_clustermap=os.path.join(opath_2_texpr_firstBinary,"Fig_"+firstBinary+"_Condition_Clustermap.png")
    ofile_secondBinary_clustermap=os.path.join(opath_2_texpr_secondBinary,"Fig_"+secondBinary+"_Condition_Clustermap.png")
    
    '''T2.3 - 1 file'''
    ofile_numberofSamplesCondition=os.path.join(opath_2_texpr_bare,"Rp_samplesPerCondition.txt")
    
    
    '''Output folder 3 (Table 3 - files: 3_msp_condition_specificity'''
    
    '''T3.1 - Max 3 files'''
    ofile_AllTSI=os.path.join(opath_3_tspec_all,"Rs_all_tauTSI.txt")
    ofile_firstBinaryTSI=os.path.join(opath_3_tspec_firstBinary,"Rs_"+firstBinary+"_tauTSI_dis.txt")
    ofile_secondBinaryTSI=os.path.join(opath_3_tspec_secondBinary,"Rs_"+secondBinary+"_tauTSI.txt")
    
    '''TT3.2 - Max 3 files'''
    ofile_all_specificMirnas=os.path.join(opath_3_tspec_all,"Rs_all_specificMiRNAs.txt")
    ofile_firstBinary_specificMirnas=os.path.join(opath_3_tspec_firstBinary,"Rs_"+firstBinary+"_specificMiRNAs.txt")
    ofile_secondBinary_specificMirnas=os.path.join(opath_3_tspec_secondBinary,"Rs_"+secondBinary+"_specificMiRNAs.txt")
    
    '''TT3.2 - Max 3 files'''
    ofile_all_specificMirnasperTissue=os.path.join(opath_3_tspec_all,"Rs_all_spMiRNAs_PerCondition.txt")
    ofile_firstBinary_specificMirnasperTissue=os.path.join(opath_3_tspec_firstBinary,"Rs_"+firstBinary+"_spMiRNAs_PerCondition.txt")
    ofile_secondBinary_specificMirnasperTissue=os.path.join(opath_3_tspec_secondBinary,"Rs_"+secondBinary+"_spMiRNAs_PerCondition.txt")
   
    '''T3.3 - Max 3 files'''  
    ofile_all_specificMirnasperTissueExtended=os.path.join(opath_3_tspec_all,"Rs_all_spMiRNAs_ConditionExpression.txt")
    ofile_firstBinary_specificMirnasperTissueExtended=os.path.join(opath_3_tspec_firstBinary,"Rs_"+firstBinary+"_spMiRNAs_ConditionExpression.txt")
    ofile_secondBinary_specificMirnasperTissueExtended=os.path.join(opath_3_tspec_secondBinary,"Rs_"+secondBinary+"_spMiRNAs_ConditionExpression.txt")
    #ofile_topTIScompareIF="C:/Users/IK/Dropbox/data/miSPan/outputFiles/3_msp_specificity/Compare_tsi_hlthyDis.txt"
    
    '''T3.4 - Max 3 files'''
    ofile_all_specificMirnasHeadmapData=os.path.join(opath_3_tspec_all,"Ds_all_spMiRNAs_HeatmapData.txt")
    ofile_firstBinary_specificMirnasHeadmapData=os.path.join(opath_3_tspec_firstBinary,"Ds_"+firstBinary+"_spMiRNAs_HeatmapData.txt")
    ofile_secondBinary_specificMirnasHeadmapData=os.path.join(opath_3_tspec_secondBinary,"Ds_"+secondBinary+"_spMiRNAs_HeatmapData.txt")
    
    '''T3.5 - Max 3 files'''
    ofile_all_specificMirnasHeadmap=os.path.join(opath_3_tspec_all,"Fig_all_spMiRNAs_Heatmap.png")
    ofile_firstBinary_specificMirnasHeadmap=os.path.join(opath_3_tspec_firstBinary,"Fig_"+firstBinary+"_spMiRNAs_Heatmap.png")
    ofile_secondBinary_specificMirnasHeadmap=os.path.join(opath_3_tspec_secondBinary,"Fig_"+secondBinary+"_spMiRNAs_Heatmap.png")
    
    
    fw=open(ofile_logfile, 'w')
    
    fw.write("Condition Specificty Analysis\n\n")
    
    fw.write("Date: "+tempFolder+"\n\n")
    #fw.write("Date: "+str(datetime.datetime.now())+"\n\n")
    #nowtime=
    
    fw.write("Input parameters and paths:\n\n")
    
    
    fw.write("Path to samples: "+iPath_samples+"\n")
    fw.write("Metadata file: "+ifile_metaData+"\n")
    fw.write("MicroRNA confidence file: "+ifile_mirConfidence+"\n")
    fw.write("Conditions: "+str(ls_conditions)+"\n")
    fw.write("Statuses: "+str(i_statusComparison)+"\n")
    fw.write("MicroRNA confidence: "+i_confidence+"\n")
    fw.write("Threshold for specific miRNAs: "+str(i_thresholdToPrint)+"\n")
    fw.write("Main output path: "+opath_main+"\n")
    
    fw.close()
    
    '''Code for cases 1, 2a, 2b''' 
    
    '''1. Data Integration component: 1_msp_data_integration. RPM dataset generation, DESeq dataset'''
    
    
    
    tup_confidence=din.dataset_integration(iPath_samples,ifile_mirConfidence, i_confidence, ifile_metaData,ls_conditions,ofile_rc_dataset, ofile_rpm_dataset)    
    '''Data Integration component END'''
    
    #sys.exit()
    
    '''2. Tissue Expression component'''
    tup_samples_to_categories=tex.samples_to_categories(ofile_rpm_dataset,ifile_metaData, i_statusComparison)
    
    '''Condition datasets: '''
    i_expressionValue='rpm'
    conditionData_rpm=cd.condition_dataset_construction(tup_samples_to_categories, i_expressionValue, i_statusComparison, ofile_allConditionDataset_rpm, ofile_firstBinaryConditionDataset_rpm, ofile_secondBinaryConditionDataset_rpm)
    
    '''(dict_mixedConditionsToAVG,sortedAtlMixedCondKeysToAVG,microRNAs, groupSizes)'''
    i_expressionValue='log2rpm'
    conditionData_log2rpm=cd.condition_dataset_construction(tup_samples_to_categories, i_expressionValue, i_statusComparison, ofile_allConditionDataset_log2rpm, ofile_firstBinaryConditionDataset_log2rpm, ofile_secondBinaryConditionDataset_log2rpm)
    
    
    
    #sys.exit()
    
    viz.clustermap_conditions_log2rpm(ofile_allConditionDataset_log2rpm, ofile_all_clustermap)
    
    if i_statusComparison:
        viz.clustermap_conditions_log2rpm(ofile_firstBinaryConditionDataset_log2rpm, ofile_firstBinary_clustermap)
        viz.clustermap_conditions_log2rpm(ofile_secondBinaryConditionDataset_log2rpm, ofile_secondBinary_clustermap)
    
    
    
    '''Analysis continues with log2 normalized data'''
    dict_mixedConditionsToAVG=conditionData_log2rpm[0]
    
        
    ls_forPrint=[]
    ls_forPrint.append(dict_mixedConditionsToAVG)
    groupSizes=conditionData_log2rpm[3]
    utl.print_group_sizes(groupSizes,ls_forPrint, ofile_numberofSamplesCondition)
    
    '''Tissue Expression component END'''
    
    
    
    
    dict_mirnaCondExpr=utl.reverse_dict(dict_mixedConditionsToAVG)  
    
    dict_mirnasTsi=csm.tau_tsi(dict_mirnaCondExpr, ofile_AllTSI)
    
    #dict_mirnasTsi=csm.tsi(dict_mirnaCondExpr, ofile_AllTSI)
    
    #dict_mirnasTsi=csm.shannon_entropy_hg(dict_mirnaCondExpr, ofile_AllTSI)
    
    
    cs.mirna_tissue_mapping(dict_mirnaCondExpr, dict_mirnasTsi,i_thresholdToPrint,ofile_all_specificMirnasperTissueExtended, ofile_all_specificMirnas)
    
    tuple_All_cond_specific_mirnas=cs.condition_specific_mirnas(ofile_all_specificMirnasperTissueExtended)
 
    viz.specific_mirna_OneDataset_Heatmap(dict_mirnaCondExpr, dict_mirnasTsi, ofile_all_specificMirnasperTissueExtended,ofile_all_specificMirnasHeadmapData, ofile_all_specificMirnasHeadmap, ofile_all_specificMirnasperTissue)
    #viz.clustermap_conditions_log2rpm(ofile_all_HeadmapData, ofile_specificMirnasClustermap)
    
    
    
    
    '''Start: Steps in order to calculate TSI separately for Healthy and Diseased'''
    if i_statusComparison:
        
        condition=firstBinary
        dict_firstBinary_mixedConditionsToAVG=utl.conditionSelection(dict_mixedConditionsToAVG, condition)
        
        condition=secondBinary
        dict_secondBinary_mixedConditionsToAVG=utl.conditionSelection(dict_mixedConditionsToAVG, condition)
        
        dict_firstBinary_mirnaCondExpr=utl.reverse_dict(dict_firstBinary_mixedConditionsToAVG)    
        dict_secondBinary__mirnaCondExpr=utl.reverse_dict(dict_secondBinary_mixedConditionsToAVG)
        
        dict_firstBinary_mirnasTsi=csm.tau_tsi(dict_firstBinary_mirnaCondExpr, ofile_firstBinaryTSI)
        dict_secondBinary_mirnasTsi=csm.tau_tsi(dict_secondBinary__mirnaCondExpr, ofile_secondBinaryTSI)
        
        cs.mirna_tissue_mapping(dict_firstBinary_mirnaCondExpr, dict_firstBinary_mirnasTsi,i_thresholdToPrint,ofile_firstBinary_specificMirnasperTissueExtended, ofile_firstBinary_specificMirnas)
        cs.mirna_tissue_mapping(dict_secondBinary__mirnaCondExpr, dict_secondBinary_mirnasTsi,i_thresholdToPrint,ofile_secondBinary_specificMirnasperTissueExtended, ofile_secondBinary_specificMirnas)
        
        tuple_firstBinary_cond_specific_mirnas=cs.condition_specific_mirnas(ofile_firstBinary_specificMirnasperTissueExtended)
        tuple_secondBinary_cond_specific_mirnas=cs.condition_specific_mirnas(ofile_secondBinary_specificMirnasperTissueExtended)
    
        viz.specific_mirna_OneDataset_Heatmap(dict_firstBinary_mirnaCondExpr, dict_firstBinary_mirnasTsi, ofile_firstBinary_specificMirnasperTissueExtended,ofile_firstBinary_specificMirnasHeadmapData, ofile_firstBinary_specificMirnasHeadmap, ofile_firstBinary_specificMirnasperTissue)
        viz.specific_mirna_OneDataset_Heatmap(dict_secondBinary__mirnaCondExpr, dict_secondBinary_mirnasTsi, ofile_secondBinary_specificMirnasperTissueExtended,ofile_secondBinary_specificMirnasHeadmapData, ofile_secondBinary_specificMirnasHeadmap, ofile_secondBinary_specificMirnasperTissue)
       
        
        
    
    
    
