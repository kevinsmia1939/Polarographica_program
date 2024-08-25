# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:30:20 2023

@author: timtic
"""

import serial 

class Serial_Communication():
    
    portsl = None
    ###############################################################################
    # The following function  will find available Arduino-ports which are able for 
    # serial communication. This function wil be embedded into a GUI later on.
    ###############################################################################
    def FindPorts():
        ports = ['COM%s' % (i+1) for i in range(256)]
        result = []
        for port in ports:
            try:
                s = serial.Serial(port)
                s.close()
                result.append(port)
            except (OSError, serial.SerialException):
                pass
        Serial_Communication.portsl = result

