'''Comments or suggestions please contact Ioannis Kavakiotis: ikavakiotis@gmail.com'''


import sys
sys.path.append('./src')
import msp_pipeline  as pl 
import os

def mispan():
    
    workingDirName = os.path.dirname(os.path.realpath(__file__)) + os.sep
    
    pl.pipeline(workingDirName)
    


if __name__ == "__main__":
    
    mispan()