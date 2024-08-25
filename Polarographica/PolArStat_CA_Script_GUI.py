# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:30:20 2023

@author: timtic
"""

########################################################################################################
########################################################################################################
#                    Importing required modules
########################################################################################################
########################################################################################################
import numpy as np
import struct
"""
########################################################################################################
# tkinter for building GUIs 
########################################################################################################"""
import tkinter as tk 
from   tkinter import ttk                           # Python 3
from   tkinter import *                             # Python 3
from   tkinter import messagebox                    # Python 3
from   tkinter import Menu                          # Python 3
from   tkinter import filedialog                    # Python 3
from   tkinter.scrolledtext  import ScrolledText    # Python 3
########################################################################################################
# matplotlib related modules for graphics
########################################################################################################
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
########################################################################################################
# Implement the default Matplotlib key bindings.
########################################################################################################
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
"""
########################################################################################################
# time for timed events, serial for serial
# communication of Arduino with the PC. The 
# serial requires the installation of pyserial
# on the operating PC.
########################################################################################################"""   
import time
import serial
########################################################################################################
from datetime import datetime        # import datetome for checking recent cal file     
########################################################################################################
# Everything for threading the Data aqcuisition out of the GUI
########################################################################################################
from    threading                                import Thread
from    PolArStat_CA_Script_Acquisition          import CA_ACQUISITION_BACKEND  as CAAQB
from    PolArStat_CA_Script_Write_Outputs        import PolArStat_CA_Write_Read as PCAWR
from    PolArStat_Serial_Communication_Functions import Serial_Communication    as SERCO
from    PolArStat_Serial_Communication_Functions import *    



def PotentiostatScript_CA():
    ####################################################################################################
    #    Define plotting options to make the plot look cool :)
    ####################################################################################################
    font = {'family': 'Times New Roman', 'color':  'black','weight': 'normal','size': 15,}
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    plt.rcParams['font.sans-serif'] = ['Times new Roman']
    ####################################################################################################
    # Initialize some global variables
    ####################################################################################################
    global arduino       ;   arduino       = None      # initialize the device as None until it is defined
    global UpCounter     ;   UpCounter     = 0
    global CalData       ;   CalData       = 0
    global Zero_IDX_E    ;   Zero_IDX_E    = 13250
    global Zero_IDX_I    ;   Zero_IDX_I    = 13250
    global Initiate_Plot ;   Initiate_Plot = True      # For first setup of a plot 
    global thread_2
    ####################################################################################################
    # Check, if there is a recent calibration file. If so, load it and use it. If not, 
    # show a respective waring and start the interface
    ####################################################################################################
    try: 
        CalFile    = datetime.today().strftime('CALIBRATION\%Y_%m_%d_Cal.txt')
        CalData    = np.genfromtxt(CalFile , skip_header = 0, skip_footer = 1)
        Neg_IDX_E  = np.average(CalData[CalData[::,1]==1,3])
        Pos_IDX_E  = np.average(CalData[CalData[::,1]==2,3])
        Neg_IDX_I  = np.average(CalData[CalData[::,1]==1,4])
        Pos_IDX_I  = np.average(CalData[CalData[::,1]==2,4])
        Zero_IDX_E = int(0.5*(Neg_IDX_E + Pos_IDX_E))
        Zero_IDX_I = int(0.5*(Neg_IDX_I + Pos_IDX_I))
    except:
            messagebox.showwarning(title="No Recent Calibration found!", message="You did not provide a recent calibration file!\n Your data might be incorrect.")
            txt_filename = ""
         
    ####################################################################################################
    # Change the variables in the class of the PolArStat_CA_Script_Write_Outputs script
    # according to the classical values or the values from the calibration file
    ####################################################################################################
    PCAWR.Zero_IDX_E = Zero_IDX_E
    PCAWR.Zero_IDX_I = Zero_IDX_I         
         
    ####################################################################################################
    # The following function asks for a filename, where the data will be stored
    ####################################################################################################
    def AskSaveFile_AndStart_CA():
        try:
            txt_filename = filedialog.asksaveasfile(defaultextension='.txt').name #öffnet fenster um nach speicherort zu fragen
        except:
            messagebox.showwarning(title="No File selected!", message="You did not select a file to save to!")
            txt_filename = ""
        else:
            outfile_CA_Data   = open(txt_filename, 'w')
            PCAWR.Output_File = outfile_CA_Data
            START_CA()
         
    ######################################################################################################
    # The following functionwill find available Arduino-ports which are able for 
    # serial communication. This function wil be embedded into a GUI later on.
    ###################################################################################################### 
    def StartSerial():
        global arduino                                                                           # TT comment: Initiate the global variable arduino
        try:
            arduino = serial.Serial(portList.get(), baudrate = 115200, timeout=.1) # TT comment: define the communication with the Arduino
            time.sleep(2)
            PCAWR.WRITE_TEXT_OUTPUT(Output_text, "Serial communication at Port %s established successfully" %portList.get())
        except:
            messagebox.showerror(title="Serial Communication failed!", 
                                 message="Unable to communicate with selected Port or no Port selected")
        else:
            randfloat  = 11.01                                      # construct a random float to send
            SEND_BYTES = b'\x44\x66' + struct.pack('f', randfloat)  # third entry of stuff to send is the packed version of the generated float
            arduino.write(SEND_BYTES)                               # Here, the bytes are send to the arduino
            time.sleep(2)                                           # sleep for 2 seconds
            SER_OUT_1  = arduino.readline()[:-2]                    # since Arduino returns the suff with println, at which end is always \r\n, skip the last to, to get what matters
            SER_OUT_2  = arduino.readline()[:-2]                    # Read the float which was send to Arduino, analogue, skip the last two entries
            if (SER_OUT_1 == SEND_BYTES) & (float(SER_OUT_2) == np.round(randfloat, decimals = 2)):
                pass                                                # pass the loop, if connection is established
            else:
                messagebox.showwarning(title="Serial Connection Corrupted!", 
                                       message="It seems like the data recieved over Serial is corrupted!")
    
    ###################################################################################################### 
    # The following function will set the inputs for a CA-measurement
    ###################################################################################################### 
    def Set_CA_Params():
        global ReadResist
        try:
            E_1          = float(E_1_Eingabe.get())
            E_2          = float(E_2_Eingabe.get())
            E_3          = float(E_3_Eingabe.get())
            E_4          = float(E_4_Eingabe.get())
            E_5          = float(E_5_Eingabe.get())
            t_1          = (1e3)*float(t_1_Eingabe.get())  # Transmission in ms
            t_2          = (1e3)*float(t_2_Eingabe.get())  # Transmission in ms
            t_3          = (1e3)*float(t_3_Eingabe.get())  # Transmission in ms
            t_4          = (1e3)*float(t_4_Eingabe.get())  # Transmission in ms
            t_5          = (1e3)*float(t_5_Eingabe.get())  # Transmission in ms
            Repetitions  = float(Repetitions_Eingabe.get())
            ReadResist   = float(R_read_Eingabe.get())
            #==============================================================================================
            # See, if there are experimental notes, specified by the user
            #==============================================================================================
            try:
                Exper_Notes = str(Exper_Notes_Eingabe.get())
                if Exper_Notes != "":
                    PCAWR.Exper_Notes = Exper_Notes
                if Exper_Notes == "":
                    PCAWR.Exper_Notes = "No notes specified"
            except:
                pass
            #==============================================================================================
            # Prepare sending the Info to Arduino
            #==============================================================================================
            SEND_E_1   = b'\x12' + b'\x17' + struct.pack('f', E_1)          # b'\x12' means modify CA params. Rest, look in Arduino script
            SEND_E_2   = b'\x12' + b'\x18' + struct.pack('f', E_2)          # b'\x12' means modify CA params. Rest, look in Arduino script
            SEND_E_3   = b'\x12' + b'\x19' + struct.pack('f', E_3)          # b'\x12' means modify CA params. Rest, look in Arduino script
            SEND_E_4   = b'\x12' + b'\x20' + struct.pack('f', E_4)          # b'\x12' means modify CA params. Rest, look in Arduino script
            SEND_E_5   = b'\x12' + b'\x21' + struct.pack('f', E_5)          # b'\x12' means modify CA params. Rest, look in Arduino script
            SEND_t_1   = b'\x12' + b'\x22' + struct.pack('f', t_1)          # b'\x12' means modify CA params. Rest, look in Arduino script
            SEND_t_2   = b'\x12' + b'\x23' + struct.pack('f', t_2)          # b'\x12' means modify CA params. Rest, look in Arduino script
            SEND_t_3   = b'\x12' + b'\x24' + struct.pack('f', t_3)          # b'\x12' means modify CA params. Rest, look in Arduino script
            SEND_t_4   = b'\x12' + b'\x25' + struct.pack('f', t_4)          # b'\x12' means modify CA params. Rest, look in Arduino script
            SEND_t_5   = b'\x12' + b'\x26' + struct.pack('f', t_5)          # b'\x12' means modify CA params. Rest, look in Arduino script
            SEND_REP   = b'\x12' + b'\x27' + struct.pack('f', Repetitions)  # b'\x12' means modify CA params. Rest, look in Arduino script
            #==============================================================================================
            # Change the variable ReadResist in the PolArStat_CA_Script_Write_Outputs script according to
            # its actual value which is used in the experiments
            #==============================================================================================
            PCAWR.ReadResist = ReadResist
            #==============================================================================================
            # Start sending the Info to Arduino
            #==============================================================================================
            arduino.write(SEND_E_1)
            time.sleep(2) 
            arduino.write(SEND_E_2)
            time.sleep(2) 
            arduino.write(SEND_E_3)
            time.sleep(2) 
            arduino.write(SEND_E_4)
            time.sleep(2) 
            arduino.write(SEND_E_5)
            time.sleep(2) 
            arduino.write(SEND_t_1)
            time.sleep(2) 
            arduino.write(SEND_t_2)
            time.sleep(2) 
            arduino.write(SEND_t_3)
            time.sleep(2) 
            arduino.write(SEND_t_4)
            time.sleep(2) 
            arduino.write(SEND_t_5)
            time.sleep(2) 
            arduino.write(SEND_REP)
            time.sleep(2) 
            #==============================================================================================
            # In case of successful transmission, print statement   
            #==============================================================================================
            PCAWR.WRITE_TEXT_OUTPUT(Output_text, "\n Data transmission complete :) \n---------------------------------------------")
        except:
            messagebox.showerror(title="Incomplete Parameters!", 
                                 message="All inputs have to be filled out in order to run a CV!")
        
    ###################################################################################################### 
    # The following function will retreive the transmitted parameters for a CA-measurement
    ###################################################################################################### 
    def Get_CA_Params():
        '''
        # the readline reads the already converted input,         
        # the[:-2] skips the last \r\n of println in Arduino. 
        # The float converts the output into a float.
        '''
        #==============================================================================================
        # Initialize an Array for daata storage
        #==============================================================================================
        Readback_Array     = np.zeros(12)     # E1, E2, E3, E4, E5, t1, t2, t3, t4, t5, nRep, Rread
        #==============================================================================================
        # Send command to Arduino for retreiving the inputs which specify the CA
        #==============================================================================================
        GETTING_COMMAND_BYTES  = b'\x13\x00\x00\x00\x00\x00'     
        arduino.write(GETTING_COMMAND_BYTES)
        #==============================================================================================
        # Wait two seconds to get back the data
        #==============================================================================================
        time.sleep(2) 
        #==============================================================================================
        # Fill the "Readback_Array" with the data obtained from the device
        #==============================================================================================
        Readback_Array[0]   = float(arduino.readline()[:-2])
        Readback_Array[1]   = float(arduino.readline()[:-2])
        Readback_Array[2]   = float(arduino.readline()[:-2])
        Readback_Array[3]   = float(arduino.readline()[:-2])
        Readback_Array[4]   = float(arduino.readline()[:-2])
        Readback_Array[5]   = (1e-3)*float(arduino.readline()[:-2])
        Readback_Array[6]   = (1e-3)*float(arduino.readline()[:-2])
        Readback_Array[7]   = (1e-3)*float(arduino.readline()[:-2])
        Readback_Array[8]   = (1e-3)*float(arduino.readline()[:-2])
        Readback_Array[9]   = (1e-3)*float(arduino.readline()[:-2])
        Readback_Array[10]  = float(arduino.readline()[:-2])
        Readback_Array[11]  = ReadResist
        #==============================================================================================
        # Write retreived data in the variable CA_Inputs of the PolArStat_CA_Script_Acquisition
        #==============================================================================================
        CAAQB.CV_Inputs     = Readback_Array 
        #==============================================================================================
        # Write retreived data to the output window
        #==============================================================================================
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"The following Parameters were set! \n")   
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_1 \t= \t %s V vs. Re"         %Readback_Array[0] )  
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_2 \t= \t %s V vs. Re"         %Readback_Array[1] )  
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_3 \t= \t %s V vs. Re"         %Readback_Array[2] )  
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_4 \t= \t %s V vs. Re"         %Readback_Array[3] )  
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_5 \t= \t %s V vs. Re"         %Readback_Array[4] )  
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t1 \t= \t %s s"                 %Readback_Array[5] )
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t1 \t= \t %s s"                 %Readback_Array[6] )
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t1 \t= \t %s s"                 %Readback_Array[7] )
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t1 \t= \t %s s"                 %Readback_Array[8] )
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t1 \t= \t %s s"                 %Readback_Array[9] )
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"n-Reps \t = \t %s "             %Readback_Array[10]) 
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"R_read \t = \t %s Ohm"          %Readback_Array[11]) 
       

    ####################################################################################################
    # The following function will start a CA measurment
    #################################################################################################### 
    def START_CA():
        #===============================================================================================
        # reset the global variable up counter, which is required for writing output-data
        #===============================================================================================
        global UpCounter    ;   UpCounter= 0
        #===============================================================================================
        # retreive CA-inputs from the device with the Get_CA_Params() function and
        # write them to the output-text field called "Output_text"
        #===============================================================================================
        Get_CA_Params()
        #===============================================================================================
        # write a header in the field "Output_text" for monitoring data
        #===============================================================================================
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"\n Preparing CA completed.\n\n") 
        PCAWR.WRITE_TEXT_OUTPUT(Output_text,"Loop.No.\t Step.No.\t t/s\t E/V\t I/mA\n")
        #===============================================================================================
        # If there is a device - named arduino - talking to the PC, clear memory and run the 
        # looper function for data acquisition
        #===============================================================================================
        if arduino is not None:   
            SEND_BYTES = b'\x14\x00\x00\x00\x00\x00'   #bytes for running a CA
            arduino.write(SEND_BYTES)
            time.sleep(3)
            Starter = arduino.readline()[:-2]
            while Starter != b'10101010':          # proceed only, if Arduino tells that conditioning-loop is done by sending b'101010'
                time.sleep(0.001)                  # Wait for one millisecond before trying again  
                Starter = arduino.readline()[:-2]  # read again and enter while loop again (if not fulfilled)
            CAAQB.CLEAR_STORAGE()
            CAAQB.DEVICE  = arduino
            CAAQB.running = True
            time.sleep(1)
            #-------------------------------------------------------------------------------------------
            CAAQB.Looper_ON = True
            Serial_Looper(plot1, plot2, canvas)
        else:
            messagebox.showwarning(title="No device in list!", 
                                   message="No device selected. Initialize serial communication before\running an experiment!")
        
    ####################################################################################################
    # Function, to stop a CA measurement. The Portlist_Refresh() function is defined in the 
    # PolArStat_Serial_Communication_Functions script and re-imported above
    ####################################################################################################   
    def Stop_CA():
        if CAAQB.running == False:
            messagebox.showwarning(title="No active measurement!", 
                                   message="There is no measurement which can be stopped!")
        else:
            if CAAQB.running == True:
                PCAWR.WRITE_TEXT_OUTPUT(Output_text,"\n Status: Measurment was manually stopped by User.\n") 
            CAAQB.StopMeasure = True
            CAAQB.Looper_ON   = False
            if arduino is not None:
                arduino.close()
            SERCO.FindPorts()    
     
    ####################################################################################################
    # The following function will kill thread_1 and exit GUI cleanly.
    ####################################################################################################
    def _quit():
        thread_2.join()
        PotentiostatWindow.quit()
        PotentiostatWindow.destroy()
        
    ####################################################################################################
    # Function, calling itself
    ####################################################################################################   
    def Serial_Looper(plot1, plot2, canvas):          # This function will be called in a loop back and forth with the main_loop of the GUI to update the output
        if CAAQB.Looper_ON == True:
            #===========================================================================================
            # Everything important for updating the plot
            #===========================================================================================
            plot1.clear()
            plot2.clear()
            InputArray = CAAQB.Storage_Array
            if len(InputArray[::,0]) <= 10000:
                plot1.plot((1e-3)*InputArray[2::,2], -0.12452*(InputArray[2::,4]-Zero_IDX_I)/ReadResist, color = 'blue' )                       # 3300/26500 = 0.12452   
                plot2.plot((1e-3)*InputArray[2::,2], -0.000249*(InputArray[2::,3]-Zero_IDX_E), color = 'red' )                                  # 6.6/26500 = 0.000249
            if len(InputArray[::,0]) > 10000:
                plot1.plot((1e-3)*InputArray[2:10:,2], -0.12452*(InputArray[2:10:,4]-Zero_IDX_I)/ReadResist, color = 'blue' )                       # 3300/26500 = 0.12452   
                plot2.plot((1e-3)*InputArray[2:10:,2], -0.000249*(InputArray[2:10:,3]-Zero_IDX_E), color = 'red' )                                  # 6.6/26500 = 0.000249
            plot1.set_xlabel("$t$ in s", fontsize = 13)
            plot1.set_ylabel("$I$ in mA", fontsize = 13)
            plot1.grid(which='both', linestyle = '--')
            plot1.tick_params(direction = 'in', length=4, width=0.5, colors='k', labelsize = 13)
            plot2.set_xlabel("$t$ in s", fontsize = 13)
            plot2.set_ylabel("$E$ vs. RE in V", fontsize = 13) 
            plot2.tick_params(direction = 'in', length=4, width=0.5, colors='k', labelsize = 13)
            canvas.draw()
            #===========================================================================================
            # Everything important for updating the serial Output window
            #===========================================================================================
            global UpCounter
            PCAWR.WRITE_DATA_OUTPUT(WHERE_TO_WRITE = Output_text, DATA_TO_WRITE = InputArray[UpCounter+1::,::])
            UpCounter = len(InputArray[::,0])
            Output_text.see(tk.END)
            #===========================================================================================
            # Callback to looper function - only if the respective bool is True
            #===========================================================================================
            PotentiostatWindow.after(100, Serial_Looper, plot1, plot2, canvas)
        
    ####################################################################################################
    # Start the conditional data-acquisition on its own thread, that the GUI does not freeze
    # if data collection is busy
    ####################################################################################################
    thread_2 = Thread(target=CAAQB.LINE_READER_FOR_THREAD)
    thread_2.start()
    ####################################################################################################
    # Run the PolArStat_CA_Script in a mainloop
    ####################################################################################################
    if __name__ == 'PolArStat_CA_Script_GUI':
        ################################################################################################
        # The following part is the main loop which calls the GUI for interfacing with 
        #the Arduino. The Arduino has to contain the firmware before this script works.
        ################################################################################################
        #===============================================================================================
        # The following part is the main loop which calls the GUI for interfacing with 
        # the Arduino. The Arduino has to contain the firmware before this script works.
        # First, the size of the GUI-window will be defined. Subsequently, all input field are set
        #===============================================================================================
        PotentiostatWindow = tk.Tk()
        PotentiostatWindow.title("Main Interface for CV-measurement with PolArStat Potentiostat")
        PotentiostatWindow.geometry('800x460')
        menu = Menu(PotentiostatWindow)                                               
        PotentiostatWindow.config(menu=menu)
        #==============================================================================================
        # Set all buttons and the dropdown menu for initializing serial communication
        #==============================================================================================
        labelTop = tk.Label(PotentiostatWindow, text = "Choose COM Port")
        labelTop.place(x = 25, y = 10)
        #-----------------------------------------------------------------------------------------------
        portsl = SERCO.portsl
        #-----------------------------------------------------------------------------------------------
        portList = ttk.Combobox(PotentiostatWindow, postcommand=lambda: portList.configure(values=SERCO.portsl))
        portList.place(x = 140, y = 10, width = 135, height = 19)
         #----------------------------------------------------------------------------------------------       
        Refresh_btn = ttk.Button(PotentiostatWindow, text="Refresh Portlist", command=SERCO.FindPorts)   
        Refresh_btn.place(x = 25, y = 35, width = 120, height = 25)
        #-----------------------------------------------------------------------------------------------
        SSer_btn = ttk.Button(PotentiostatWindow, text="Start Serial", command=StartSerial)   
        SSer_btn.place(x = 155, y = 35, width = 120, height = 25)
        #==============================================================================================
        # Set all the input buttons for CV measurement and place them in the GUI
        #==============================================================================================
        E_1_Label = Label(PotentiostatWindow,text="E1 vs. Ref [V]*")
        E_1_Label.place(x = 25, y = 70)
        E_1_Eingabe = Entry(PotentiostatWindow)
        E_1_Eingabe.insert(END, -0.5) 
        E_1_Eingabe.place(x = 110, y = 70, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        E_2_Label = Label(PotentiostatWindow,text="E2 vs. Ref [V]*")
        E_2_Label.place(x = 25, y = 95)
        E_2_Eingabe = Entry(PotentiostatWindow)
        E_2_Eingabe.insert(END, 0.5) 
        E_2_Eingabe.place(x = 110, y = 95, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        E_3_Label = Label(PotentiostatWindow,text="E3 vs. Ref [V]*")
        E_3_Label.place(x = 25, y = 120)
        E_3_Eingabe = Entry(PotentiostatWindow)
        E_3_Eingabe.insert(END, 0.0) 
        E_3_Eingabe.place(x = 110, y = 120, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        E_4_Label = Label(PotentiostatWindow,text="E4 vs. Ref [V]*")
        E_4_Label.place(x = 25, y = 145)
        E_4_Eingabe = Entry(PotentiostatWindow)
        E_4_Eingabe.insert(END, 0.0) 
        E_4_Eingabe.place(x = 110, y = 145, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        E_5_Label = Label(PotentiostatWindow,text="E5 vs. Ref [V]*")
        E_5_Label.place(x = 25, y = 170)
        E_5_Eingabe = Entry(PotentiostatWindow)
        E_5_Eingabe.insert(END, 0.0) 
        E_5_Eingabe.place(x = 110, y = 170, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        t_1_Label = Label(PotentiostatWindow,text="t1 [s]*")
        t_1_Label.place(x = 175, y = 70)
        t_1_Eingabe = Entry(PotentiostatWindow)
        t_1_Eingabe.insert(END, 10) 
        t_1_Eingabe.place(x = 222, y = 70, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        t_2_Label = Label(PotentiostatWindow,text="t2 [s]*")
        t_2_Label.place(x = 175, y = 95)
        t_2_Eingabe = Entry(PotentiostatWindow)
        t_2_Eingabe.insert(END, 10) 
        t_2_Eingabe.place(x = 222, y = 95, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        t_3_Label = Label(PotentiostatWindow,text="t3 [s]*")
        t_3_Label.place(x = 175, y = 120)
        t_3_Eingabe = Entry(PotentiostatWindow)
        t_3_Eingabe.insert(END, 0) 
        t_3_Eingabe.place(x = 222, y = 120, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        t_4_Label = Label(PotentiostatWindow,text="t4 [s]*")
        t_4_Label.place(x = 175, y = 145)
        t_4_Eingabe = Entry(PotentiostatWindow)
        t_4_Eingabe.insert(END, 0) 
        t_4_Eingabe.place(x = 222, y = 145, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        t_5_Label = Label(PotentiostatWindow,text="t5 [s]*")
        t_5_Label.place(x = 175, y = 170)
        t_5_Eingabe = Entry(PotentiostatWindow)
        t_5_Eingabe.insert(END, 0) 
        t_5_Eingabe.place(x = 222, y = 170, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        Repetitions_Label = Label(PotentiostatWindow,text="Repetitions*")
        Repetitions_Label.place(x = 25, y = 215)
        Repetitions_Eingabe = Entry(PotentiostatWindow)
        Repetitions_Eingabe.insert(END, 1)
        Repetitions_Eingabe.place(x = 150, y = 215)
        #-----------------------------------------------------------------------------------------------
        R_read_Label = Label(PotentiostatWindow,text="R read [Ohm]*")
        R_read_Label.place(x = 25, y = 240)
        R_read_Eingabe = Entry(PotentiostatWindow)
        R_read_Eingabe.insert(END, 120)
        R_read_Eingabe.place(x = 150, y = 240)
        #==============================================================================================
        # Define all buttons in the PolArStat CA-GUI
        #==============================================================================================
        CASet_btn = ttk.Button(PotentiostatWindow, text="Send CA inputs", command=Set_CA_Params)   
        CASet_btn.place(x = 25, y = 370, width = 120, height = 25)
        #-----------------------------------------------------------------------------------------------
        CAGet_btn = ttk.Button(PotentiostatWindow, text="Check CA inputs", command=Get_CA_Params)   
        CAGet_btn.place(x = 155, y = 370, width = 120, height = 25)
        #-----------------------------------------------------------------------------------------------
        CARun_btn = ttk.Button(PotentiostatWindow, text="Run CA", command=AskSaveFile_AndStart_CA)   
        CARun_btn.place(x = 25, y = 405, width = 120, height = 50)
        #-----------------------------------------------------------------------------------------------
        CAStop_btn = ttk.Button(PotentiostatWindow, text="Stop CA", command=Stop_CA)   
        CAStop_btn.place(x = 155, y = 405, width = 120, height = 50)
        #==============================================================================================
        # Define all text-containing windows in the GUI
        #==============================================================================================
        Exper_Notes_Label = Label(PotentiostatWindow,text="Experimental Notes")
        Exper_Notes_Label.place(x = 25, y = 270)
        Exper_Notes_Eingabe = Entry(PotentiostatWindow)
        Exper_Notes_Eingabe.place(x = 25, y = 295, width = 250, height = 60)
        #-----------------------------------------------------------------------------------------------
        Output_text_Label = Label(PotentiostatWindow,text="Monitor serial output data")
        Output_text_Label.place(x = 325, y = 350)
        Output_text = tk.scrolledtext.ScrolledText(PotentiostatWindow,  wrap = tk.WORD,  font = ("Times New Roman", 9)) 
        Output_text.place(x = 325, y = 370,  width = 450,  height = 85) 
        Output_text.focus() 
        #==============================================================================================
        # Initiate the Plot in the main GUI of the PolArStat CA-Script
        #==============================================================================================
        if Initiate_Plot == True:
            global plot1
            global plot2
            global canvas
            fig = Figure(figsize = (6.5,4.5), dpi = 68)
            plot1 = fig.add_subplot(111)
            x = np.array([0])
            y = np.array([0])
            plot1.plot((x,y), linewidth = 0)                     # 3300/26500 = 0.12452   
            plot1.set_xlabel("$t$ in s", fontsize = 13)
            plot1.set_ylabel("$I$ in mA", fontsize = 13, color = 'blue')
            plot1.grid(which='both', linestyle = '--')
            plot1.tick_params(direction = 'in', length=4, width=0.5, colors='k', labelsize = 13)
            plot2 = plot1.twinx()
            plot2.set_ylabel("$E$ vs. RE in V", fontsize = 13, color = 'red')
            plot2.tick_params(direction = 'in', length=4, width=0.5, colors='k', labelsize = 13)
            canvas = FigureCanvasTkAgg(fig, PotentiostatWindow)  
            canvas.draw()
            canvas.get_tk_widget().place(x = 325, y = 10)
            Initiate_Plot = False
        #==============================================================================================
        # Run the Potentiostat window in a mainloop
        #==============================================================================================
        PotentiostatWindow.mainloop()
        #==============================================================================================
        # If mainloop is killed, kill port (here named arduino)
        #==============================================================================================
        arduino.close()   
    

        
        























    