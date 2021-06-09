import os
import shutil
import sys
def parse_config_file(working_dir):
    
    
    
    ifile=working_dir+"cfg"+os.sep+"cfg.txt"
    #print ifile
    
    data={}
    
    f = open (ifile)
    
    for l in f.readlines():
        
        firstChar=l[0]
            
        if firstChar == '%':
            continue
        
        key=""
        value=""
        
        #print l
        col=l.strip().split("=")
        
        if col[0] == 'mirConfidence_file':
            key=col[0]
            value=col[1]
        
        if col[0] == 'expression_repository_path':
            key=col[0]
            value=col[1]
            
        if col[0] == 'metaData_file':
            key=col[0]
            value=col[1]
            
        if col[0] == 'mirna_confidence':
            key=col[0]
            value=col[1]
            
        if col[0] == 'conditions':
            key=col[0]
            value=col[1].strip().split("#")    
        
        if col[0] == 'status':
            key=col[0]
            if col[1].strip().split("#")[0]:
                value=col[1].strip().split("#")
            else:
                value=[]
            #value=col[1].strip().split("#")  
            
        if col[0] == 'tsi_threshold':
            key=col[0]
            value=col[1]           
        
        
        data[key]=value
        
    return data
       
     
