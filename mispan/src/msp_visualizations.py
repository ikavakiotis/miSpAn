

import pandas as pd
import numpy as np
import seaborn as sns; 
import matplotlib as plt



def specific_mirna_OneDataset_Heatmap(dict_mirCondAVGLog2Expr, dict_mirnaTSI, ifile_tissueSpecificity,ofile_specificMirnasHeatmapData, ofile_specificMirnasHeatmap, ofile_tissuespecificMirnas):

    
    
    f = open (ifile_tissueSpecificity)
   
    
    rr=-1
    
    microRNAsTSI=""
    tissuesExpr=""
    
    for l in f.readlines():
        rr=rr+1
        
        if rr==0:
            microRNAsTSI=l
        
        if rr==1:
            tissuesExpr=l
        
        if rr==2:
            break
        
    
    
    ls_microRNAsTSI=microRNAsTSI.strip().split("\t")
    ls_tissuesExpr=tissuesExpr.strip().split("\t")
    
    ls_mirnas=[]
    
    ls_mirnasHeader=[]
    
    dict_ContLsMir={}
    
    dict_ContLsMir_dis={}
    dict_ContLsMir_hlthy={}
    
    dict_mirnaMirnaHeader={}
    
    dict_tissueSpMirnas={}
    
    
    for i in range(0,(len(ls_microRNAsTSI))):
        
        tsi_mirna=ls_microRNAsTSI[i]
        tissueExpr=ls_tissuesExpr[i]
        
        mirna=tsi_mirna.split("_",1)[1]
        tissue=tissueExpr.split("_",1)[1]
        
        ls_mirnas.append(mirna)
        
        #disState=tissue.split("_",1)[0]
        
        #mirnaHeader=mirna+"_"+disState
        
        
        if tissue in dict_tissueSpMirnas:
            ls_spMir=dict_tissueSpMirnas[tissue]
            ls_spMir.append(mirna)
            dict_tissueSpMirnas[tissue]=ls_spMir
        else:
            ls_spMir=[]
            ls_spMir.append(mirna)
            dict_tissueSpMirnas[tissue]=ls_spMir
        
        #print mirnaHeader
        
        #dict_mirnaMirnaHeader[mirna]=mirnaHeader
        
        
        #ls_mirnasHeader.append(mirnaHeader)
        
        #if disState == "Dis":
        if tissue in dict_ContLsMir:
            ls_mir=dict_ContLsMir[tissue]
            ls_mir.append(mirna)
            dict_ContLsMir[tissue]=ls_mir
        else:
            ls_mir=[]
            ls_mir.append(mirna)
            dict_ContLsMir[tissue]=ls_mir
        '''
        elif disState == "Hlthy":
            if tissue in dict_ContLsMir_hlthy:
                ls_mir=dict_ContLsMir_hlthy[tissue]
                ls_mir.append(mirna)
                dict_ContLsMir_hlthy[tissue]=ls_mir
            else:
                ls_mir=[]
                ls_mir.append(mirna)
                dict_ContLsMir_hlthy[tissue]=ls_mir
        '''
    #or tsi_mirna in ls_microRNAsTSI
    '''
    for i, k in dict_ContLsMir_hlthy.iteritems():
        print i+" "+str(len(k))+" "+str(k)
    '''
    
    df_data=pd.DataFrame();
    r=0
    ls_cond=[]
    
    
    mir=ls_mirnas[0]
    #print mir
    
    #for mir in ls_mirnas:
         
    dict_condAVGL2=dict_mirCondAVGLog2Expr[mir]
        
       
        
    
    for mixedCond, avglog2 in dict_condAVGL2.iteritems():
        
        ls_cond.append(mixedCond)
        
    #print str(ls_cond)
            #Conditions=ls_disCond+ls_hlthyCond
    df_data['Conditions']=ls_cond
            
       
    for mirna in ls_mirnas:
        
        #print mirna
        ls_acgLog2Column=[]
        dict_condAVGL2=dict_mirCondAVGLog2Expr[mirna]
        
        for condit in ls_cond:
            
            expr=dict_condAVGL2[condit]
            ls_acgLog2Column.append(expr)
        
        
        #mirHead=dict_mirnaMirnaHeader[mirna]
        df_data[mirna]=ls_acgLog2Column 
    
    
    #if r==2:
    #sys.exit()
        #ls_avgLog2Column=ls_disAvgLog2+ls_hlthyAvgLog2
        
        
    
        
    #print df_data.info()
    #print df_data.head(5)
        #df_data.set_index('Tissues')
    df_data.to_csv(ofile_specificMirnasHeatmapData,sep="\t",index=False)
    
    df_condExpression=pd.read_csv(ofile_specificMirnasHeatmapData,sep= "\t", index_col=0) 
    sns.set(color_codes=True)
    
    df_Transposed = df_condExpression.T
    
    cmap = sns.cm.rocket_r
    
    #species = iris.pop("species")
    #g = sns.clustermap(df_Transposed,row_cluster=False, col_cluster=False, figsize=(50,200),cmap=cmap)
    g = sns.clustermap(df_Transposed,row_cluster=False, col_cluster=False, cmap=cmap)
    
    #plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    #plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    
    g.savefig(ofile_specificMirnasHeatmap)
    
    
    
    fw = open(ofile_tissuespecificMirnas, "w")
    
    for tis, ls_m in dict_tissueSpMirnas.iteritems():
        fw.write(tis+"\t"+str(len(ls_m))+"\t"+"\t".join(ls_m)+"\n")
       
    fw.close()
    
    '''
    fw = open(ofile_CounttissuespecificMirnas, "w")
    
    header="Tissue\tDiseased\tHealthy\n"
    fw.write(header)
    
    ls_bareTis=[]
    for tis, ls_m in dict_tissueSpMirnas.iteritems():

        baretis=tis.split("_",1)[1]
        
        if baretis in ls_bareTis:
            continue
        
        dis_tis="Dis_"+baretis
        hlthy_tis="Hlthy_"+baretis
        
        
        if  dis_tis in dict_tissueSpMirnas:
            dis_count=len(dict_tissueSpMirnas[dis_tis])
        else:
            dis_count=0
            
        if  hlthy_tis in dict_tissueSpMirnas:
            hlthy_count=len(dict_tissueSpMirnas[hlthy_tis])
        else:
            hlthy_count=0
        
        ls_bareTis.append(baretis)
        
        fw.write(baretis+"\t"+str(dis_count)+"\t"+str(hlthy_count)+"\n")
        
    fw.close()
    '''
    



def clustermap_conditions_log2rpm(ifile_condAVGLog2RPM, ofile_figure):
    
    df_condExpression=pd.read_csv(ifile_condAVGLog2RPM,sep= "\t", index_col=0) 
    sns.set(color_codes=True)
    
    df_Transposed = df_condExpression.T
    
    
    #print df_Transposed.info
    
    cmap = sns.cm.rocket_r
    
    #species = iris.pop("species")
    g = sns.clustermap(df_Transposed,row_cluster=True, col_cluster=False, figsize=(50,200),cmap=cmap)
    
    #g = sns.clustermap(df_Transposed,row_cluster=False, col_cluster=False, cmap=cmap)
    
    g.savefig(ofile_figure)
    
    '''
    
    df_2=df_Transposed.iloc[0:30, :]
    
    df_3=df_Transposed.iloc[0:80, :]
    
    g2 = sns.clustermap(df_2,row_cluster=False, col_cluster=False, figsize=(50,50),cmap=cmap)
    
    #g = sns.clustermap(df_Transposed,row_cluster=False, col_cluster=False, cmap=cmap)
    
    g2.savefig(ofile_figure2)
    
    
    g3 = sns.clustermap(df_3,row_cluster=False, col_cluster=False, figsize=(50,50),cmap=cmap)
    
    #g = sns.clustermap(df_Transposed,row_cluster=False, col_cluster=False, cmap=cmap)
    
    g3.savefig(ofile_fig_3)
    '''


def circos_dataset(ituple_condition_specific_mirnas, maxMir, ofile_circosData,ofile_circosDataRound):
    '''old specificMicroRNAs to change se print circos'''
    
    dict_tissuemirExpr=ituple_condition_specific_mirnas[0]
    ls_microRNAs=ituple_condition_specific_mirnas[1]
    
    
    fw=open(ofile_circosData,"w")
    
    fwabs=open(ofile_circosDataRound,"w")
    
    
    header="Tissues\t"+"\t".join(ls_microRNAs)          
    
    fw.write(header+"\n")
    fwabs.write(header+"\n")
    
    maxMir=int(maxMir)
    
    if maxMir==0:
        maxMir=len(ls_microRNAs)+1
    
    
    for tissue, dict_mirExpr in dict_tissuemirExpr.iteritems():
        #print tissue
        numMir=0
        line=tissue.replace(" ","_")
        line2=tissue.replace(" ","_")
        #print line
        for mir in ls_microRNAs:
            
            numMir=numMir+1
            if numMir==maxMir:
                continue
            
            if mir in dict_mirExpr:
                expr=dict_mirExpr[mir]
                expr2=dict_mirExpr[mir]
                expr2=round(float(expr2))
                expr2=int(expr2)
            else:
                expr="-"
                expr2="-"
            
            line=line+"\t"+expr
            line2=line2+"\t"+str(expr2)
    
        line=line.strip()+"\n"
        line2=line2.strip()+"\n"
        fw.write(line)  
        fwabs.write(line2)
    #print len(ls_microRNAsTSI)
    
    fwabs.close()

    dataTuple=(dict_tissuemirExpr,)
    
    return dataTuple


