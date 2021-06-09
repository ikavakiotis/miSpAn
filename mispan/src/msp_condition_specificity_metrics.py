
import sys
import numpy as np


def shannon_entropy_hg(dict_mirCondEpr,ofile_tsi):

    
    
    dict_mirnashg={}
    
    for mirna, dict_condExpr in dict_mirCondEpr.iteritems():
        
        
        '''For each microRNA collect all expression values. One for each condition'''
        ls_expression=[]
        for cond, expression in dict_condExpr.iteritems():
            ls_expression.append(expression)
            #print cond+str(ls_expression)
        maxVal=max(ls_expression)
        summary_cond=sum(ls_expression)
        #sys.exit()
        if maxVal==0.0:
            continue
        
        #if mirna=='hsa-let-7a-2-3p':
        #    print "This: "+str(ls_expression)+" "+ str(maxVal)+" "+str(summary_cond)
        #print  mirna+ " "+str(maxVal) 
        
        n_cond=0
        summary=0
        
        for cond, expression in dict_condExpr.iteritems():
            
            '''How many condition exist'''
            n_cond=n_cond+1
            
            pi=float(expression)/float(summary_cond)
            log2pi=np.log2(pi)
            
            formula=pi+log2pi
            
            summary=summary+ (-1)*formula
        
        hg=summary
        
        
        transHg=1-hg/np.log2(n_cond)
        
        dict_mirnashg[mirna]=transHg
        
    fw = open (ofile_tsi, "w")
    
    
    d_view = [ (tsi,miRNA) for miRNA,tsi in dict_mirnashg.iteritems() ]
    d_view.sort(reverse=True) # natively sort tuples by first element
    for tsi,miRNA in d_view:
        fw.write(str(miRNA)+"\t"+str(tsi)+"\n")
        
    fw.close()
    
    return dict_mirnashg


def tsi(dict_mirCondEpr,ofile_tsi):
    
    
    
    dict_mirnasTsi={}
    
    for mirna, dict_condExpr in dict_mirCondEpr.iteritems():
        
        
        '''For each microRNA collect all expression values. One for each condition'''
        ls_expression=[]
        for cond, expression in dict_condExpr.iteritems():
            ls_expression.append(expression)
            #print cond+str(ls_expression)
        maxVal=max(ls_expression)
        summary=sum(ls_expression)
        
        #if mirna=='hsa-let-7a-2-3p':
        #    print str(ls_expression)+" "+ str(maxVal)+" "+str(summary)
        #sys.exit()
        if maxVal==0.0:
            continue
        #print  mirna+ " "+str(maxVal) 
        tsi=float(maxVal)/float(summary)
              
        dict_mirnasTsi[mirna]=tsi
        
    fw = open (ofile_tsi, "w")
    
    
    d_view = [ (tsi,miRNA) for miRNA,tsi in dict_mirnasTsi.iteritems() ]
    d_view.sort(reverse=True) # natively sort tuples by first element
    for tsi,miRNA in d_view:
        fw.write(str(miRNA)+"\t"+str(tsi)+"\n")
        
    fw.close()
    
    return dict_mirnasTsi


def tau_tsi(dict_mirCondEpr,ofile_tsi):

    
    
    dict_mirnasTsi={}
    
    for mirna, dict_condExpr in dict_mirCondEpr.iteritems():
        
        
        '''For each microRNA collect all expression values. One for each condition'''
        ls_expression=[]
        for cond, expression in dict_condExpr.iteritems():
            ls_expression.append(expression)
            #print cond+str(ls_expression)
        maxVal=max(ls_expression)
        #sys.exit()
        if maxVal==0.0:
            continue
        #print  mirna+ " "+str(maxVal) 
        
        n_cond=0
        summary=0
        
        for cond, expression in dict_condExpr.iteritems():
            
            '''How many condition exist'''
            n_cond=n_cond+1
            
            summary=summary+(1-float(expression)/float(maxVal))
        
        tsi=float(summary)/float(n_cond-1)
        
        
        dict_mirnasTsi[mirna]=tsi
        
    fw = open (ofile_tsi, "w")
    
    
    d_view = [ (tsi,miRNA) for miRNA,tsi in dict_mirnasTsi.iteritems() ]
    d_view.sort(reverse=True) # natively sort tuples by first element
    for tsi,miRNA in d_view:
        fw.write(str(miRNA)+"\t"+str(tsi)+"\n")
        
    fw.close()
    
    return dict_mirnasTsi