# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:30:20 2023

@author: timtic
"""

import  numpy     as np
import  time
from    PolArStat_CV_Script_Write_Outputs        import PolArStat_CV_Write_Read as PCVWR


class ACQUISITION_BACKEND():
    ###############################################################################
    #   Initialize some variables in the class "ACQUISITION_BACKEND()"
    ###############################################################################
    running        = False
    Read_Finished  = False
    StopMeasure    = False
    Looper_ON      = False
    DEVICE         = None
    CV_Inputs      = None
    Storage_Array  = np.zeros((1,5))
    
    ###############################################################################
    #   Reset variables, if they have been modified
    ###############################################################################
    def CLEAR_STORAGE():
        ACQUISITION_BACKEND.Storage_Array  = np.zeros((1,5))   # clear memory
        ACQUISITION_BACKEND.Looper_ON      = False             # reset variable
        ACQUISITION_BACKEND.StopMeasure    = False             # reset variable
        ACQUISITION_BACKEND.running        = False             # reset variable
        ACQUISITION_BACKEND.Read_Finished  = False             # reset variable
    
    ###############################################################################
    #    The following will read the line of a serial output, provided by
    #    the device in action (given that there is a device connected).  
    ###############################################################################      
    def LINE_READER_FOR_THREAD():
        while True :
            if ACQUISITION_BACKEND.running == False:
                time.sleep(0.05)
            if ACQUISITION_BACKEND.running == True and ACQUISITION_BACKEND.Read_Finished == False:
                try:
                    data = ACQUISITION_BACKEND.DEVICE.readline()[:-2]
                    if data:
                        if data != b'999999':   # As soon as the Arduino sends 999999, the measurement is done
                            DECODED       = np.array(data.decode("utf-8").split('\t'))
                            DECODED_FLOAT = DECODED.astype(np.float)
                            ACQUISITION_BACKEND.Storage_Array = np.vstack([ACQUISITION_BACKEND.Storage_Array, DECODED_FLOAT])
                        if data == b'999999':
                            ACQUISITION_BACKEND.running           = False
                            ACQUISITION_BACKEND.Looper_ON         = False
                            ACQUISITION_BACKEND.Read_Finished     = True
                            PCVWR.SAVE_DATA_AFTER_STOP(STOP_STATE = "Success", EXP_PARAMS = ACQUISITION_BACKEND.CV_Inputs, DATA = ACQUISITION_BACKEND.Storage_Array)
                            ACQUISITION_BACKEND.DEVICE.close()
                except:
                    ACQUISITION_BACKEND.running     = False
                    ACQUISITION_BACKEND.Looper_ON   = False
                    ACQUISITION_BACKEND.StopMeasure = True   
                    PCVWR.SAVE_DATA_AFTER_STOP(STOP_STATE = "Interrupt_or_fail", EXP_PARAMS = ACQUISITION_BACKEND.CV_Inputs, DATA = ACQUISITION_BACKEND.Storage_Array)
                    ACQUISITION_BACKEND.DEVICE.close()
                    
    
    ###############################################################################
    #    The following function takes an array and reduces it by its first comuln 
    #    indices by averaging equal ones. This is what BioLogic does when it says
    #    average over n-percent of the step. Here, we just have (around 4 points).
    ###############################################################################    
    def ARRAYCONDENSER(ARRAY):
        NLB   = 0
        Out   = np.zeros(len(ARRAY[0,1::]))
        for i in range(len(ARRAY[::,0])-1):
            if ARRAY[i+1,0] != ARRAY[i,0]:
                B   = np.mean(ARRAY[NLB:i+1,1::], axis = 0)
                B[0] = ARRAY[i,1]
                Out = np.vstack([Out,B])        
                NLB = i+1
        return Out[1::,::]
                   
            
            
                
