
import time

from pandas import *
import numpy as np
import sys



def condition_specific_mirnas(ifile_tissueSpecificity):
    '''msp'''
    '''old: tissuesSpMirs'''
    
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
    
    
    dict_tissuemirExpr={}
    
    
    ls_microRNAs=[]
    ls_tissues=[]
    
    for i in range(0,len(ls_microRNAsTSI)):
        #print ls_microRNAsTSI[i]
        mirTSI=ls_microRNAsTSI[i]
        tissueExpr=ls_tissuesExpr[i]
        
        #print mirTSI+" "+tissueExpr
        
        mir=mirTSI.split("_")[1]
        tissue=tissueExpr.split('_', 1)[1]
        expr=tissueExpr.split('_', 1)[0]
        
        ls_microRNAs.append(mir)
        ls_tissues.append(tissue)
        #print mir+" "+tissue+" "+expr
        
        if tissue in dict_tissuemirExpr:
            
            dict_mirExpr=dict_tissuemirExpr[tissue]
            dict_mirExpr[mir]=expr
            dict_tissuemirExpr[tissue]=dict_mirExpr
            
        else:
            dict_mirExpr={}
            dict_mirExpr[mir]=expr
            dict_tissuemirExpr[tissue]=dict_mirExpr
            
    '''
    for tisss, dict_mirEx in dict_tissuemirExpr.iteritems():
        
        for mir, expr in dict_mirEx.iteritems():
            
            print tisss+" "+mir+" "+expr
    '''
    
    #print len(dict_tissuemirExpr)
    
    dataTuple=(dict_tissuemirExpr,ls_microRNAs)
    
    return dataTuple




def mirna_tis_comparison(dict_mircondavg_dis,dict_mircondavg_hlthy,dict_mirTSI_dis,dict_mirTSI_hlthy,threshold,ituple_confidence,ofile_topTIScomparison):
    
    '''Used only when we have aclculated separate tsi for two conditions. Is it right? I don think so.'''
    
    
    dict_mirConfidence=ituple_confidence[0]    

    threshold=float(threshold)
    
    
    ls_miRNAhighTSI=[]
    
    
    df_data=DataFrame();
    
    d_view = [ (tsi,miRNA) for miRNA,tsi in dict_mirTSI_dis.iteritems() ]
    d_view.sort(reverse=True) # natively sort tuples by first element
    for tsi,miRNA in d_view:
        
        if tsi>threshold:
            ls_miRNAhighTSI.append(miRNA)
        else:
            break
    
    ls_Cond=[]
    ls_hlthyCond=[]
    
    r=0
    
    for mir in ls_miRNAhighTSI:
         
        dict_condAVGL2_dis=dict_mircondavg_dis[mir]
        
        dict_condAVGL2_hlthy=dict_mircondavg_hlthy[mir]
        
        
        if mir in dict_mirConfidence:
            mirConf=dict_mirConfidence[mir]
        else:
            mirConf="unknown"
        
        
        #if mirConf != "high":
            #continue
        r=r+1
        #print "here"
        #ls_disAvgLog2=[]
        #ls_hlthyAvgLog2=[]
        ls_avgLog2Column=[]
        
        for dis_mixedCond, dis_avglog2 in dict_condAVGL2_dis.iteritems():
            
            #print dis_mixedCond
            
            ls_avgLog2Column.append(dis_avglog2)
            
            
            
            #htly_mixedCond="Hlthy_"+dis_mixedCond.split("_")[1]
            htly_mixedCond=dis_mixedCond.replace("Dis_","Hlthy_")
            
            #print htly_mixedCond
            
        
            if r==1:
                ls_Cond.append(dis_mixedCond)
                ls_Cond.append(htly_mixedCond)
            
            
            hlthy_avglog2=dict_condAVGL2_hlthy[htly_mixedCond]
            ls_avgLog2Column.append(hlthy_avglog2)
            
        if r==1:
            #Conditions=ls_disCond+ls_hlthyCond
            df_data['Tissues']=ls_Cond
            
            
        #ls_avgLog2Column=ls_disAvgLog2+ls_hlthyAvgLog2
        
        
        mirAndCond=mir+"_"+mirConf
        df_data[mirAndCond]=ls_avgLog2Column 
        
        
    
    #print df_data.info()
    #print df_data.head(5)
    #df_data.set_index('Tissues')
    df_data.to_csv(ofile_topTIScomparison,sep="\t",index=False)


def mirna_tissue_mapping(dict_mircondavg, dict_mirTSI,thresholdS,ofile_mirtsiConditions, ofile_specificMirnas):
    
        
    threshold=float(thresholdS)
    
    df_data=DataFrame();
    
    d_view = [ (tsi,miRNA) for miRNA,tsi in dict_mirTSI.iteritems() ]
    d_view.sort(reverse=True)
    # natively sort tuples by first element
    
    f=open(ofile_specificMirnas, 'w')
    f.write("miRNA\tTSI\tMost epxressed condition\tExpression\n")
    
    for tsi,miRNA in d_view:
        strTowrite=""
        #print tsi
        if tsi<threshold:
            break
        
        
        listForColumn=[]
        fTsi=float(tsi)
        rtsi=round(fTsi,3)
        tsi_mirna=str(rtsi)+"_"+miRNA
        #listForColumn.append(tsi_mirna)
        
        strToWrite=miRNA+"\t"+str(fTsi)+"\t"
        
        dict_condAVGExpr=dict_mircondavg[miRNA]
        
        d_viewCond = [(expr,cond) for cond,expr in dict_condAVGExpr.iteritems() ]
        d_viewCond.sort(reverse=True)
        
        rr=0
        for expr,cond in d_viewCond:
            fexpr=float(expr)
            texpr=round(fexpr,3)
            expr_cond=str(texpr)+"_"+cond
            
            listForColumn.append(expr_cond)
            if rr==0:
               strToWrite=strToWrite+cond+"\t"+str(fexpr)+"\n"
               f.write(strToWrite)
               
            rr=rr+1
            
        df_data[tsi_mirna]=listForColumn
        
        
    f.close()
        
        #print df_data.info()
        
        
        #fw.write(str(miRNA)+"\t"+str(tsi)+"\n")
    
    
    df_data.to_csv(ofile_mirtsiConditions,sep="\t",index=False)
    
    