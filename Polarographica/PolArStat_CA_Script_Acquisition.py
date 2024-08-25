# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:56:57 2023

@author: gisel
"""

import  numpy     as np
import  time
from    PolArStat_CA_Script_Write_Outputs        import PolArStat_CA_Write_Read as PCAWR


class CA_ACQUISITION_BACKEND():
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
        CA_ACQUISITION_BACKEND.Storage_Array  = np.zeros((1,5))   # clear memory
        CA_ACQUISITION_BACKEND.Looper_ON      = False             # reset variable
        CA_ACQUISITION_BACKEND.StopMeasure    = False             # reset variable
        CA_ACQUISITION_BACKEND.running        = False             # reset variable
        CA_ACQUISITION_BACKEND.Read_Finished  = False             # reset variable
    
    ###############################################################################
    #    The following will read the line of a serial output, provided by
    #    the device in action (given that there is a device connected).  
    ###############################################################################      
    def LINE_READER_FOR_THREAD():
        while True :
            if CA_ACQUISITION_BACKEND.running == False:
                time.sleep(0.05)
            if CA_ACQUISITION_BACKEND.running == True and CA_ACQUISITION_BACKEND.Read_Finished == False:
                try:
                    data = CA_ACQUISITION_BACKEND.DEVICE.readline()[:-2]
                    if data:
                        if data != b'999999':   # As soon as the Arduino sends 999999, the measurement is done
                            DECODED       = np.array(data.decode("utf-8").split('\t'))
                            DECODED_FLOAT = DECODED.astype(np.float)
                            CA_ACQUISITION_BACKEND.Storage_Array = np.vstack([CA_ACQUISITION_BACKEND.Storage_Array, DECODED_FLOAT])
                        if data == b'999999':
                            CA_ACQUISITION_BACKEND.running        = False
                            CA_ACQUISITION_BACKEND.Looper_ON      = False
                            CA_ACQUISITION_BACKEND.Read_Finished  = True
                            PCAWR.SAVE_DATA_AFTER_STOP(STOP_STATE = "Success", EXP_PARAMS = CA_ACQUISITION_BACKEND.CV_Inputs, DATA = CA_ACQUISITION_BACKEND.Storage_Array)
                            CA_ACQUISITION_BACKEND.DEVICE.close()
                except:
                    CA_ACQUISITION_BACKEND.running        = False
                    CA_ACQUISITION_BACKEND.Looper_ON      = False
                    CA_ACQUISITION_BACKEND.StopMeasure    = True   
                    PCAWR.SAVE_DATA_AFTER_STOP(STOP_STATE = "Interrupt_or_fail", EXP_PARAMS = CA_ACQUISITION_BACKEND.CV_Inputs, DATA = CA_ACQUISITION_BACKEND.Storage_Array)
                    CA_ACQUISITION_BACKEND.DEVICE.close()
                    