

import numpy as np
import sys


def print_group_sizes(dict_groupSizes,ls_conditionDicts, ofile_groupsizes):
    
    
    fw=open(ofile_groupsizes,"w")
    
    
    fw.write("Number of samples for each condition\n\n")
    
    for condict in ls_conditionDicts:
        
        for cont in condict.keys():
        
            if cont in dict_groupSizes:
                size=dict_groupSizes[cont]
                fw.write(cont+"\t"+str(size)+"\n")
                
    fw.close()


def reverse_dict(dict_condMiRNAExpr):
    '''msp'''
    '''input: 
    a dict: 
        key:condition/tissue
        value: dict: mirna:avglo2Expr

    output:
       a dict:
           key: mirna
           value: dict: condition:avgLog2Expr
    '''

    dict_miRNACondExpr={}
   
    for condition, dict_mirExpr in dict_condMiRNAExpr.iteritems():
        
        
        for mirna, expression in dict_mirExpr.iteritems():
            
            if mirna in dict_miRNACondExpr:
                dict_contExpr=dict_miRNACondExpr[mirna]
                dict_contExpr[condition]=expression
                dict_miRNACondExpr[mirna]=dict_contExpr
            else:
                dict_contExpr={}
                dict_contExpr[condition]=expression
                dict_miRNACondExpr[mirna]=dict_contExpr    
    
    return dict_miRNACondExpr
  
   


def same_conditions(dict_dis,dict_hlthy):
    '''mep'''
    dict_newdis={}
    dict_newhlthy={}
    
    
    #print str(len(dict_dis.keys()))+" "+str(dict_dis.keys())
    #print str(len(dict_hlthy.keys()))+" "+str(dict_hlthy.keys())
   
    
    ls_hlthyKeys=dict_hlthy.keys()
    ls_disKeys=dict_dis.keys()
    
    #ls_hlthyKeys=["Hlthy_a","Hlthy_b","Hlthy_c","Hlthy_d","Hlthy_w"]
    #ls_disKeys=["Dis_b", "Dis_d","Dis_f","Dis_e"]
    
    ls_hlthyBareKeys=[]
    ls_disBareKeys=[]
    
    for k in ls_hlthyKeys:
        barek=k.split("_")[1]
        ls_hlthyBareKeys.append(barek)
        
    for k in ls_disKeys:
        barek=k.split("_")[1]
        ls_disBareKeys.append(barek)
        
    commonKeys=[value for value in ls_disBareKeys if value in ls_hlthyBareKeys] 
    
    #print str(commonKeys)
    
    
    for bareKey in commonKeys:
        hlthyKey='Hlthy_'+bareKey
        disKey='Dis_'+bareKey
        
        hlthyVal=dict_hlthy[hlthyKey]
        dict_newhlthy[hlthyKey]=hlthyVal
        
        disVal=dict_dis[disKey]
        dict_newdis[disKey]=disVal
        
 
    
    dataTuple=(dict_newdis,dict_newhlthy)
    
    #print str(len(dict_seldis.keys()))+" "+str(dict_seldis.keys())
    #print str(len(dict_hlthy.keys()))+" "+str(dict_hlthy.keys())
    
    #sys.exit()
    return dataTuple

def conditionSelection(dict_condmirExpr, condition):
    
    dict_selectedData={}
    
    for cond, dict_mirExpr in  dict_condmirExpr.iteritems():
        prefix=cond.split("_")[0]
        if prefix == condition :
            dict_selectedData[cond]=dict_mirExpr

    return dict_selectedData

def sort_AlternatelyFirstSecondBinary(dict_mixedConditionsToAVG, i_statusComparison):
    '''different from the intitial'''
    
    if i_statusComparison:
        firstBin=i_statusComparison[0]
        secondBin=i_statusComparison[1]
    else:
        firstBin=""
        secondBin=""
    
    finalList=[]
    singleStatus=[]
    
    set_bareConditions=set()
    
    ls_mixedConditions=dict_mixedConditionsToAVG.keys() 
    
    if not i_statusComparison:
        #print "sort_AlternatelyFirstSecondBinary: status false" + str(ls_mixedConditions)
        return ls_mixedConditions
    #ls_mixedConditions=["Hlthy_liver", "Dis_atom", "Dis_brain","Dis_liver", "Hlthy_two","Hlthy_brain"]
    for condition in ls_mixedConditions:
        
        #h_sample="Hlthy_"+sample
        
        bareCondition=condition.split("_")[1]
        
        if bareCondition not in set_bareConditions:
            #print "sort_AlternatelyFirstSecondBinary: bare conditions "+bareCondition
            seconB=secondBin+"_"+bareCondition
            firstB=firstBin+"_"+bareCondition
            set_bareConditions.add(bareCondition)
        
            if (seconB in ls_mixedConditions) and (firstB in ls_mixedConditions): 
                finalList.append(firstB)
                finalList.append(seconB)
            elif (seconB in ls_mixedConditions) and (firstB not in ls_mixedConditions):
                singleStatus.append(seconB)
            elif (seconB not in ls_mixedConditions) and (firstB in ls_mixedConditions):
                singleStatus.append(firstB)
        
    finalList=finalList+singleStatus
    
    '''Debug: print list'''    
    #print str(finalList)
    
    return finalList


def from_expressionlist_to_avg(dict_mixedConditions):
    '''msp'''
    
    '''Input: dict with key: mixed condition and valueu a dict with mirna an list of expressions
        Output: the same dict replaced expresion list with average expression
    '''
    dict_groupSizes={}
    dict_mixedConditionsmiRNAToAVG={}
    
    for condition, dict_miRNAToListExpressions in dict_mixedConditions.iteritems():
        #print condition
        dict_miRNAAVG={}
        for miRNA, ls_epxr in dict_miRNAToListExpressions.iteritems():
            
            if condition not in dict_groupSizes:
                dict_groupSizes[condition]=len(ls_epxr)
            
            avg=sum(ls_epxr)/float(len(ls_epxr))
            dict_miRNAAVG[miRNA]=avg
            
            '''Debug: avg'''
            #print condition+" "+miRNA+" "+str(avg)
        
        dict_mixedConditionsmiRNAToAVG[condition]=dict_miRNAAVG
            
    
    dataTuple=(dict_mixedConditionsmiRNAToAVG, dict_groupSizes)    
      
    return dataTuple


def from_expressionlist_to_log2avg(dict_mixedConditions):
    
    '''msp'''
    '''Input: dict with key: mixed condition and valueu a dict with mirna an list of expressions
        Output: the same dict replaced expresion list with average log2 expression
    '''
    
    
    #dict_mixedConditions=iTuple_data[0]
    
    lowCutoff=1.0
    
    #print len(dict_mixedConditions)
    #print dict_mixedConditions.keys()
    dict_groupSizes={}
    dict_mixedConditionmiRNALog2AVG={}
    
    for mixedCondition, dict_miRNAListExpressions in dict_mixedConditions.iteritems():
        i=0
        #print "\n\n\n"+mixedCondition
        
        dict_miRNAAVGLog2Expression={}
        
        for mirRNA, ls_expresion in dict_miRNAListExpressions.iteritems():
            
            #print mirRNA
            
            if mixedCondition not in dict_groupSizes:
                dict_groupSizes[mixedCondition]=len(ls_expresion)
            
            ls_expresion_new=[]
            for exprstr in ls_expresion:
                expr=float(exprstr)
                if expr>lowCutoff:
                    ls_expresion_new.append(expr)
                else:
                    ls_expresion_new.append(1.001)
            
            #if na valw log2 kai meta mean
            if ls_expresion_new:
                
                #print str(ls_expresion_new)
                ls_expresion_log2=[np.log2(x) for x in ls_expresion_new]
                #print str(ls_expresion_log2)
                avgLog2Expr=float(sum(ls_expresion_log2))/float(len(ls_expresion_log2))
                #print str(avgLog2Expr)
                dict_miRNAAVGLog2Expression[mirRNA]=avgLog2Expr
                '''Debug: avg'''
                #print mixedCondition+" "+mirRNA+" "+str(avgLog2Expr)
                i=i+1
        #print mixedCondition+" numExpressed: "+str(i)
        dict_mixedConditionmiRNALog2AVG[mixedCondition]=dict_miRNAAVGLog2Expression
        
    dataTuple=dict_mixedConditionmiRNALog2AVG
    dataTuple=(dict_mixedConditionmiRNALog2AVG, dict_groupSizes)
    
    return dataTuple


 
            