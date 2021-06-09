

import pandas as pd

import time

from os import listdir
from os.path import isfile, join

import numpy as np

import sys


'''Main function running all data integration'''
def dataset_integration(iPath_samples,ifile_mirConfidence, i_confidence, ifile_metadata, ls_conditions, ofile_RC_dataset, ofile_RPM_dataset):
    '''msp'''
    
    tup_conf_sel=confidence_selected_mirnas(ifile_mirConfidence, i_confidence)
    
    conf_sel_mir=tup_conf_sel[1]
    
    tup_conf=tup_conf_sel[0]
    
    expression_dataset_construction(iPath_samples, conf_sel_mir, ifile_metadata, ls_conditions,ofile_RC_dataset, ofile_RPM_dataset)
    
    return tup_conf



'''Funtion for selecting mirnas based on confidence'''
def confidence_selected_mirnas(ifile_mirConfidence, i_confidence):
    '''msp'''
    
    ituple_confidence=parse_miRNA_confidence_file(ifile_mirConfidence)
    
    confindenceSelectedmiRNAs=[]
    if i_confidence == "high":
        confindenceSelectedmiRNAs=ituple_confidence[1]
    elif i_confidence =="low":
        confindenceSelectedmiRNAs=ituple_confidence[2]
    else:
        confindenceSelectedmiRNAs=ituple_confidence[3]
  
    tup_data=(ituple_confidence,confindenceSelectedmiRNAs)
    
    return tup_data
    
'''Auxiliary function for confidence_selected_mirnas'''
def parse_miRNA_confidence_file(ifile_mirConfidence):
    '''msp'''
    
    #start_time = time.time()
    dict_mirConfidence={}
    ls_high=[]
    ls_low=[]
    ls_all=[]
    
    f = open (ifile_mirConfidence)
    
    rr=-1
    
    for l in f.readlines():
        rr=rr+1
        if rr==0:
            continue
        
        col=l.strip().split("\t")
        
        matureMirName=col[1]	
        matureConfidence=col[2]
        
        if matureConfidence == "high":
            ls_high.append(matureMirName)
        elif matureConfidence == "low":
            ls_low.append(matureMirName)
        else:
            ls_all.append(matureMirName)
            
  
        dict_mirConfidence[matureMirName]=matureConfidence
        
  
    dataTuple=(dict_mirConfidence,ls_high,ls_low,ls_all)
    
    return dataTuple




def expression_dataset_construction(iPath_samples, confidenceSelectedMirnas, ifile_metadata, ls_conditions, ofile_rc_dataset, ofile_rpm_dataset):
    '''msp'''
    
    ls_miRNAs=[]
    ls_samples=[]
    #dict_mirnaExpression_rc={}
    dict_sampleExpression_rc={}
    
    #dict_mirnaExpression={}
    dict_sampleExpression_rpm={}
    
    df_metadata=pd.read_csv(ifile_metadata, sep="\t")
    
    #df_metadata.columns=['ID', 'condition']
    if ls_conditions:
        df_metadata_selected = df_metadata[df_metadata['Condition'].isin(ls_conditions)]
    else:
        df_metadata_selected=df_metadata
    #print df_metadata_selected.head(3)
    
    samplesMetadata = df_metadata_selected.iloc[:,0].tolist()
    
    #print "expression_dataset_construction. length of metadata: "+ str(len(samplesMetadata))
    
    #print samplesMetadata
    
    '''read files form path- each file is a single case with 2656 mirna expressions'''
    onlyfiles = [f for f in listdir(iPath_samples) if isfile(join(iPath_samples, f))]
    
    #print "expression_dataset_construction. num files: "+ str(len(onlyfiles))
    #print onlyfiles
    
    '''Counter for files'''
    numfile=0
    
    for sampleFile in onlyfiles:
        
        
        
        '''Get sampleID from file name'''
        sampleID=sampleFile.split("_")[0]
        sampleID=sampleID.split(".")[0]
        
        if sampleID not in samplesMetadata:
            continue
        
        
        '''Counter for files'''
        numfile=numfile+1
       
        '''Create the list with all sampleIDs'''
        if sampleID not in ls_samples:
            ls_samples.append(sampleID)
    
        '''Read each file to extract expression values'''
        dataFile=iPath_samples+sampleFile
        #print dataFile
        
        
        f = open (dataFile)
        
        '''Line counter'''
        rr=0
    
        for l in f.readlines():
            
            rr=rr+1
            
            '''Skip the header'''
            if rr==1:
                continue
           
            col=l.strip().split("\t")
            
            miRNA=col[0]
            
            '''If the selection mode is false, or if its TRUE but the mirna is in the selected list, 
            exclude itfrom the analysis'''
            if  confidenceSelectedMirnas and (miRNA not in confidenceSelectedMirnas):
                continue
            
            readCount=col[1]
            rpm=col[2]
            
            
            if miRNA not in ls_miRNAs:
                ls_miRNAs.append(miRNA)
            
            
            if sampleID in dict_sampleExpression_rc:
                dict_mirnaExpression_rc=dict_sampleExpression_rc[sampleID]
                dict_mirnaExpression_rc[miRNA]=readCount   
            else:
                dict_mirnaExpression_rc={}
                dict_mirnaExpression_rc[miRNA]=readCount
                    
            dict_sampleExpression_rc[sampleID]=dict_mirnaExpression_rc
            
            
            if sampleID in dict_sampleExpression_rpm:
                dict_mirnaExpression_rpm=dict_sampleExpression_rpm[sampleID]
                dict_mirnaExpression_rpm[miRNA]=rpm   
            else:
                dict_mirnaExpression_rpm={}
                dict_mirnaExpression_rpm[miRNA]=rpm
                    
            dict_sampleExpression_rpm[sampleID]=dict_mirnaExpression_rpm
                 
    #print "expression_dataset_construction. Minras in dataset: "+str(len(ls_miRNAs))
    
    '''Dataset file genaration'''
    fw = open (ofile_rc_dataset, "w")
        
    fw.write('SampleID\t')
    fw.write('\t'.join(str(p) for p in ls_miRNAs)) 
    fw.write('\n')
        
    for sample in ls_samples:
            
        toPrint=sample+"\t"    
        dict_mirnaExpression= dict_sampleExpression_rc[sample]
            
        for mirna in ls_miRNAs:
                
            value=dict_mirnaExpression[mirna]
            toPrint=toPrint+value+"\t"
            
        fw.write(toPrint.strip(' \t\n\r')+"\n") 
        
    fw.close()
    
    
    df_data=pd.read_csv(ofile_rc_dataset, sep="\t", header=0)
    
    #print df_data.info()
    
    #df_data.set_index('SampleID',inplace=True)
    
    df_datadeseq=df_data.T
    
    #df_datadeseq = df_datadeseq.iloc[1:]
    
    
    #print df_datadeseq.head(3)
    
    df_datadeseq.to_csv(ofile_rc_dataset,sep="\t", index=True, header=False)
    
    
    
    fw = open (ofile_rpm_dataset, "w")
        
    fw.write('SampleID\t')
    fw.write('\t'.join(str(p) for p in ls_miRNAs)) 
    fw.write('\n')
        
    for sample in ls_samples:
            
        toPrint=sample+"\t"    
        dict_mirnaExpression= dict_sampleExpression_rpm[sample]
            
        for mirna in ls_miRNAs:
                
            value=dict_mirnaExpression[mirna]
            toPrint=toPrint+value+"\t"
            
        fw.write(toPrint.strip(' \t\n\r')+"\n") 
        
    fw.close()
    #print("--- %s seconds ---" % (time.time() - start_time))
         

