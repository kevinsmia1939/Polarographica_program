# coding: utf-8

#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================
#Author:       Tim Tichter
#Date:         2019.07.16
#Function:     POLAROGRAPHICAS Cyclic Voltammetry Functions

#=======================================================================================================================
#=======================================================================================================================
#importing all required modules from Python
#=======================================================================================================================
#=======================================================================================================================
from tkinter                           import *
from tkinter.filedialog                import askopenfilename
from tkinter.filedialog                import asksaveasfilename
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases          import key_press_handler
from matplotlib.figure                 import Figure
from scipy.interpolate                 import interp1d
from scipy.interpolate                 import InterpolatedUnivariateSpline
from scipy.optimize                    import curve_fit
from scipy.optimize                    import nnls
from scipy.linalg                      import toeplitz
from scipy.special                     import kv, iv, gamma
from scipy.integrate                   import quad
from scipy                             import fft
from scipy.signal                      import hilbert
from cmath                             import *
from math                              import ceil,floor
import mpmath as mp
mp.dps = 25; mp.pretty = True

#=======================================================================================================================
#=======================================================================================================================
#define some global variables, functions and warning functions
#=======================================================================================================================
#=======================================================================================================================

global F
global R
global File_Was_Loaded
F = 96485.0
R = 8.314
File_Was_Loaded = 0
def cot(phi):
    return 1.0/tan(phi)
def csc(phi):
    return 1.0/sin(phi)
def coth(x):
    return 1/tanh(x)

#=======================================================================================================================    
#=======================================================================================================================
#=======================================================================================================================



#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================
def BaseGetter():
    global BaseData
    Delim = Delimiter
    
    if Delim == 1:
       
        BasePath = askopenfilename()
        BaseData = np.genfromtxt(BasePath, delimiter='\t')

        
    if Delim == 0:
       
        BasePath = askopenfilename()
        BaseData = np.genfromtxt(BasePath, delimiter='')
        
        
    global BasePotenzial
    global BaseStrom
    
    
    if FirstComesxx  == 1:
        BasePotenzial     =  BaseData[FromRowxx:ToRowxx:Readeveryxx,0::2]  * UmrPot              
        BaseStrom         =  BaseData[FromRowxx:ToRowxx:Readeveryxx,1::2]  * UmRStrom

    if FirstComesxx  == 0:
        BasePotenzial     =  BaseData[FromRowxx:ToRowxx:Readeveryxx,1::2]  * UmrPot
        BaseStrom         =  BaseData[FromRowxx:ToRowxx:Readeveryxx,0::2]  * UmRStrom
        
    if BasePotenzial[0,0] > BasePotenzial[1,0]:
        BasePotenzial     = BasePotenzial[::-1]
        BaseStrom         = BaseStrom[::-1]


# In[ ]:

def Open_KL_File():
    root = Toplevel()
    root.title("Your Data")
    global File_Was_Loaded
    File_Was_Loaded = 1
    global ExpData   #needed later to differentiate experimental and simulated data of ACCV
    ExpData      = 1
    Delim        = Delimiter
    
    global data
    
    
    if Delim == 1:
       
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='\t')

        
    if Delim == 0:
       
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='')
    
    
    #============================================================================
    #Plotting of loaded file
    #============================================================================
    
    #============================================================================
    #define entry window
    #============================================================================
    
    f = Figure(figsize=(5, 5), dpi=100)
    b = f.add_subplot(111)
        
  
        
    #============================================================================
    
    global Potenzial
    global Strom
    
    
    if FirstComesxx  == 1:
        Potenzial     =  data[FromRowxx:ToRowxx:Readeveryxx,0::2] * UmrPot            
        Strom         =  data[FromRowxx:ToRowxx:Readeveryxx,1::2] * UmRStrom
    if FirstComesxx  == 0:
        Potenzial     =  data[FromRowxx:ToRowxx:Readeveryxx,1::2] * UmrPot
        Strom         =  data[FromRowxx:ToRowxx:Readeveryxx,0::2] * UmRStrom
            
    
    
    global Weite
    Weite = Potenzial.shape[0]
    global NumMeas
    NumMeas = Potenzial.shape[1]
    
    for i in range(Weite):
        Stromarrays             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))                          
        Potenzialarrays         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:])) 
        LenPotArrays            = len(Potenzialarrays)
        
        
        b.plot(Potenzial, 0.001*Strom, linestyle='-',marker='',color='k')
        b.set_xlabel('E vs. Ref. / V', fontsize=12)
        b.set_ylabel('I / mA', fontsize=12)
        for axis in ['top','bottom','left','right']:
            b.spines[axis].set_linewidth(2)
            b.spines[axis].set_color('k')
        b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
         
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    KouteckyLevichWindow()
    
    


# In[ ]:

def Get_KL_Data():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Get-Koutecky-Levich Data File")                         
    Fenster.geometry("400x270")                                            

    FromRowxxx_label      = Label(Fenster,text="Read from row")
    FromRowxxx_label.grid(row=0, column=0)
    FromRowxxx_Eingabe    = Entry(Fenster)                                               
    FromRowxxx_Eingabe.grid(row=0, column=1)

    ToRowxxx_label       = Label(Fenster,text="Read to row")
    ToRowxxx_label.grid(row=1, column=0)
    ToRowxxx_Eingabe     = Entry(Fenster)
    ToRowxxx_Eingabe.grid(row=1, column=1)
    
    Readeveryxxx_label   = Label(Fenster,text="Read every")
    Readeveryxxx_label.grid(row=2, column=0)
    Readeveryxxx_Eingabe = Entry(Fenster)
    Readeveryxxx_Eingabe.grid(row=2, column=1)
    
    
    UmRStrom_Label = Label(Fenster,text="Current Factor to be Microampere")
    UmRStrom_Label.grid(row=3, column=0)
    UmRStrom_Eingabe = Entry(Fenster)
    UmRStrom_Eingabe.grid(row=3, column=1)
        
    UmrPot_Label = Label(Fenster,text="Potential Factor to be Volt")
    UmrPot_Label.grid(row=4, column=0)
    UmrPot_Eingabe = Entry(Fenster)
    UmrPot_Eingabe.grid(row=4, column=1)
    

    
        
    DesiReactxxx_Label  = Label(Fenster,text="Desired Reaction")
    DesiReactxxx_Label.grid(row=5, column=0)
    var1 = IntVar()
    Checkbutton(Fenster, text="Oxidation", variable=var1).grid(row=5,column=1, sticky=W)
    var2 = IntVar()
    Checkbutton(Fenster, text="Reduction", variable=var2).grid(row=5,column=2, sticky=W)
    
    
    
    #First Comes
    FirstComesxxx_Label  = Label(Fenster,text="First comes")
    FirstComesxxx_Label.grid(row=6, column=0)
    var3 = IntVar()
    Checkbutton(Fenster, text="Potential", variable=var3).grid(row=6,column=1, sticky=W)
    var4 = IntVar()
    Checkbutton(Fenster, text="Current", variable=var4).grid(row=6,column=2, sticky=W)
    
        
    Delimiterxxx_Label  = Label(Fenster,text="Delimiter")
    Delimiterxxx_Label.grid(row=7, column=0)
    
    var5 = IntVar()
    Checkbutton(Fenster, text="Tab", variable=var5).grid(row=7,column=1, sticky=W)
    var6 = IntVar()
    Checkbutton(Fenster, text="Space", variable=var6).grid(row=7,column=2, sticky=W)



    def AcceptKL():
              
        global FromRowxx
        FromRowxx     = (int(FromRowxxx_Eingabe.get()))
        global ToRowxx
        ToRowxx       = (int(ToRowxxx_Eingabe.get()))
        global Readeveryxx
        Readeveryxx   = (int(Readeveryxxx_Eingabe.get()))
        global UmRStrom#
        UmRStrom      = (float(UmRStrom_Eingabe.get()))
        global UmrPot#
        UmrPot        = (float(UmrPot_Eingabe.get()))
        
        global DesiReactxx
        global FirstComesxx

        
        DesiReactxx   = var1.get()
        FirstComesxx  = var3.get()
        
        global Delimiter
        Delimiter         = var5.get()
    
    
    def NextKL():
        Open_KL_File()
        
        def quit():
            Fenster.destroy()
        quit()
    
    
    
    
    Accept = Button(Fenster, text="Accept",command=AcceptKL)
    Accept.grid(row=8, column=0) 
    
    Next = Button(Fenster, text="Next",command=NextKL)
    Next.grid(row=9, column=0)  
    
#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================



def KouteckyLevichWindow():
    
    
    Fenster = Toplevel()                                                         
    Fenster.title("Koutecky-Levich Curve-Getter")                         
    Fenster.geometry("400x400")
    
    OCP_Label = Label(Fenster,text="OCPs vs E \nRef in V")
    OCP_Label.grid(row=0, column=0)
    
    RotRates_Label = Label(Fenster,text="Rotation rates \n in rpm")
    RotRates_Label.grid(row=0, column=1)
    
    Concent_Label = Label(Fenster,text="Concentrations as fract.\n of the highest one")
    Concent_Label.grid(row=0, column=2)
    
    #IF ONLY ONE OCP
    #______________________________________
    var1 = IntVar()
    Checkbutton(Fenster, text="One OCP", variable=var1).grid(row=80, column=0, sticky=W)
    OneOCP_Eingabe = Entry(Fenster)                                               
    OneOCP_Eingabe.grid(row=81, column=0)
    
    
    #IF ONLY ONE Rotation rate
    #______________________________________
    var2 = IntVar()
    Checkbutton(Fenster, text="One Rotation", variable=var2).grid(row=80, column=2, sticky=W)
    OneRot_Eingabe = Entry(Fenster)                                               
    OneRot_Eingabe.grid(row=81, column=2)
    

    OCPs     = []
    RotRates = []
    Concents = []

    for i in range(NumMeas):
        
        en = Entry(Fenster)
        en.grid(row=i+1, column=0)
        OCPs.append(en)

        en1 = Entry(Fenster)
        en1.grid(row=i+1, column=1)
        RotRates.append(en1)
        
        en2 = Entry(Fenster)
        en2.grid(row=i+1, column=2)
        Concents.append(en2)
    
    
    
    def GetOCPs():    
        OnlyOneOCP = var1.get()

        OCPots = []
        
        if OnlyOneOCP ==0:
            for entry in OCPs:
                OCPots.append(float(entry.get()))
        
        if OnlyOneOCP ==1:
            for j in range(NumMeas):
                OCPots.append(float(OneOCP_Eingabe.get()))

        global OCPArray
        OCPArray = np.array(OCPots) 

    
    
    def GetRots():    
        Rots= []
        for entry in RotRates:
            Rots.append(float(entry.get()))
        global RotRatesArray
        RotRatesArray = np.array(Rots) 
        global VarParRot
        global VarParCon
        VarParRot = 1
        VarParCon = 0
        
        
    def GetConcents():    
        OnlyOneRot = var2.get()
        
        if OnlyOneRot ==1:
            global RotRate
            RotRate = (float(OneRot_Eingabe.get()))
        
        Concs= []
        for entry in Concents:
            Concs.append(float(entry.get()))
        global ConcentrationsArray
        ConcentrationsArray = np.array(Concs) 
        global VarParRot
        global VarParCon
        VarParRot = 0
        VarParCon = 1
        
   
        
    def Next():
        KL_NextLevel2()
        
        def quit():
            Fenster.destroy()
        quit()
    
    
    
    
    button=Button(Fenster,text="Accept OCPs",command=GetOCPs).grid(row=82,column=0)
    button=Button(Fenster,text="Accept RotRates",command=GetRots).grid(row=82,column=1)
    button=Button(Fenster,text="Accept Concentrations",command=GetConcents).grid(row=82,column=2)
    button=Button(Fenster,text="Next",command=Next).grid(row=83,column=1)


# In[7]:

def KL_NextLevel2():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Koutecky-Levich Analysis Param-Getter")                         
    Fenster.geometry("400x400")

    
    Getter_Label = Label(Fenster, text="Get From KL")
    Getter_Label.grid(row=0, column=0)
    var1 = IntVar()
    Checkbutton(Fenster, text="get n", variable=var1).grid(row=0, column=1, sticky=W)
    var2 = IntVar()
    Checkbutton(Fenster, text="get c", variable=var2).grid(row=0,column=2, sticky=W)
    var3 = IntVar()
    Checkbutton(Fenster, text="get D", variable=var3).grid(row=0,column=3, sticky=W)
    
        
    n_Label = Label(Fenster, text="n")
    n_Label.grid(row=1, column=0)
    n_Eingabe = Entry(Fenster)                                               
    n_Eingabe.grid(row=1, column=1)

    A_Label = Label(Fenster,text="A in cm^2")
    A_Label.grid(row=2, column=0)
    A_Eingabe = Entry(Fenster)
    A_Eingabe.grid(row=2, column=1)
    
    c_Label = Label(Fenster,text="c in mol/cm^3")
    c_Label.grid(row=3, column=0)
    c_Eingabe = Entry(Fenster)
    c_Eingabe.grid(row=3, column=1)
    
    D_Label = Label(Fenster,text="D in cm^2/s")
    D_Label.grid(row=4, column=0)
    D_Eingabe = Entry(Fenster)
    D_Eingabe.grid(row=4, column=1)
    
    viscosity_Label = Label(Fenster,text="viscosity in cm^2/s")
    viscosity_Label.grid(row=5, column=0)
    viscosity_Eingabe = Entry(Fenster)
    viscosity_Eingabe.grid(row=5, column=1)
    
    T_Label = Label(Fenster,text="T in °C")
    T_Label.grid(row=6, column=0)
    T_Eingabe = Entry(Fenster)
    T_Eingabe.grid(row=6, column=1)
    

    
    def AcceptParams():
        global n
        global A
        global T 
        global Visc
        global D
        global c
        global get_n
        global get_c
        global get_D
        
    
        get_n = var1.get()
        get_c = var2.get()
        get_D = var3.get()
        
        if get_n ==1:
            if get_c ==0:
                if get_D ==0:
              
                    #n         = (float(n_Eingabe.get()))
                    A         = (float(A_Eingabe.get()))
                    T         = (float(T_Eingabe.get())) + 273.15
                    Visc      = (float(viscosity_Eingabe.get()))
                    D         = (float(D_Eingabe.get()))
                    c         = (float(c_Eingabe.get()))
        
        if get_n ==0:
            if get_c ==1:
                if get_D ==0:
              
                    n         = (float(n_Eingabe.get()))
                    A         = (float(A_Eingabe.get()))
                    T         = (float(T_Eingabe.get())) + 273.15
                    Visc      = (float(viscosity_Eingabe.get()))
                    D         = (float(D_Eingabe.get()))
                    #c         = (float(c_Eingabe.get()))
        
        if get_n ==0:
            if get_c ==0:
                if get_D ==1:
              
                    n         = (float(n_Eingabe.get()))
                    A         = (float(A_Eingabe.get()))
                    T         = (float(T_Eingabe.get())) + 273.15
                    Visc      = (float(viscosity_Eingabe.get()))
                    #D         = (float(D_Eingabe.get()))
                    c         = (float(c_Eingabe.get()))
    
    
    
    def Next():
        KL_NextLevel3()
        
        def quit():
            Fenster.destroy()
        quit()   
        
    def Back():
        KouteckyLevichWindow()
        
        def quit():
            Fenster.destroy()
        quit()   
    
    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=8,column=0)
    button=Button(Fenster,text="Next",command=Next).grid(row=9,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=9,column=0)
    
    


# In[8]:

def KL_NextLevel3():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Koutecky-Levich Analysis Level 3")                         
    Fenster.geometry("400x500")
    
    
    
    
    KLStart_Label = Label(Fenster,text="KL-StartPot")
    KLStart_Label.grid(row=0, column=0)
    KLStart_Eingabe = Entry(Fenster)
    KLStart_Eingabe.grid(row=0, column=1)
    
    KLEnd_Label = Label(Fenster,text="KL-EndPot")
    KLEnd_Label.grid(row=1, column=0)
    KLEnd_Eingabe = Entry(Fenster)
    KLEnd_Eingabe.grid(row=1, column=1)
    
    REFF_Label = Label(Fenster,text="Refining factor")
    REFF_Label.grid(row=2, column=0)
    REFF_Eingabe = Entry(Fenster)
    REFF_Eingabe.grid(row=2, column=1)
    
    Ru_Label = Label(Fenster,text="Ru in Ohm")
    Ru_Label.grid(row=3, column=0)
    Ru_Eingabe = Entry(Fenster)
    Ru_Eingabe.grid(row=3, column=1)
    
    
    #RefToCotti
    #_________________________
    
    
    if get_n == 1:
    
    
        var1 = IntVar()
        Checkbutton(Fenster, text="Ref to Cottrell", variable=var1).grid(row=4, column=0, sticky=W)

    
        CottiSlope_Label = Label(Fenster,text="Cott Slope microamp/s^0.5")
        CottiSlope_Label.grid(row=5, column=0)
        CottiSlope_Eingabe = Entry(Fenster)
        CottiSlope_Eingabe.grid(row=5, column=1)
    
      
    #Smoothing
    #________________
    var4 = IntVar()
    Checkbutton(Fenster, text="Smooth data", variable=var4).grid(row=6, column=0, sticky=W)
    
    SmoothFact_Label = Label(Fenster,text="Smoothing factor")
    SmoothFact_Label.grid(row=7, column=0)
    SmoothFact_Eingabe = Entry(Fenster)
    SmoothFact_Eingabe.grid(row=7, column=1)
    
    SmoothFrom_Label = Label(Fenster,text="Smooth from")
    SmoothFrom_Label.grid(row=8, column=0)
    SmoothFrom_Eingabe = Entry(Fenster)
    SmoothFrom_Eingabe.grid(row=8, column=1)
    
    SmoothTo_Label = Label(Fenster,text="Smooth to")
    SmoothTo_Label.grid(row=9, column=0)
    SmoothTo_Eingabe = Entry(Fenster)
    SmoothTo_Eingabe.grid(row=9, column=1)
    
    ManuCorrCurr_Label = Label(Fenster,text="Cap. Curr. microamp")
    ManuCorrCurr_Label.grid(row=12, column=0)
    ManuCorrCurr_Eingabe = Entry(Fenster)
    ManuCorrCurr_Eingabe.grid(row=12, column=1)
    
    #Corrections and sep Plots
    #__________________________
    
    var5 = IntVar()
    Checkbutton(Fenster, text="Auto Correction", variable=var5).grid(row=10, column=0, sticky=W)
    var6 = IntVar()
    Checkbutton(Fenster, text="Separate KL-Plot", variable=var6).grid(row=14, column=0, sticky=W)
    var7 = IntVar()
    Checkbutton(Fenster, text="Manual Correction", variable=var7).grid(row=11, column=0, sticky=W)
    var8 = IntVar()
    Checkbutton(Fenster, text="Baseline Correction", variable=var8).grid(row=13, column=0, sticky=W)
    var9 = IntVar()
    Checkbutton(Fenster, text="Separate LSV-Plot", variable=var9).grid(row=14, column=1, sticky=W)
    var10 = IntVar()
    Checkbutton(Fenster, text="Show uncorr LSV", variable=var10).grid(row=15, column=0, sticky=W)
    var11 = IntVar()
    Checkbutton(Fenster, text="Save as txt", variable=var11).grid(row=15, column=1, sticky=W)
    
    
    
    

    def AcceptParams():
        
        global KLStart 
        global KLEnd        
        global REFF
        global Ru 
        global RefToCotti    
        global CottiSlope    
        global Cotti_D     
        global Cotti_c 
        global Smooth_Data  
        global SmoothFact     
        global SmoothFrom    
        global SmoothTo     
        global AutoCorr 
        global OnlyKLPlot
        global Instead_n_D
        global Instead_n_c
        global ManuCorr
        global ManuCorrCurr
        global OnlyLSVPlot
        global BaseCorr
        global ShowUncorrLSV
        global AsTxtSaver
        
        
        
        KLStart    = (float(KLStart_Eingabe.get()))
        KLEnd      = (float(KLEnd_Eingabe.get()))
        REFF       = (float(REFF_Eingabe.get()))
        Ru         = (float(Ru_Eingabe.get()))

        
        
        Smooth_Data   = var4.get()
        AutoCorr      = var5.get()
        OnlyKLPlot    = var6.get()
        ManuCorr      = var7.get()
        BaseCorr      = var8.get()
        OnlyLSVPlot   = var9.get()
        ShowUncorrLSV = var10.get()
        AsTxtSaver    = var11.get()
        
        
        if ManuCorr == 1:
            ManuCorrCurr  = (float(ManuCorrCurr_Eingabe.get()))
        
        
        
        if  Smooth_Data == 1:
            SmoothFact  = (float(SmoothFact_Eingabe.get()))
            SmoothFrom  = (float(SmoothFrom_Eingabe.get()))
            SmoothTo    = (float(SmoothTo_Eingabe.get()))
        
        
        if get_n == 1:
            
            RefToCotti = var1.get()
            
            if RefToCotti  == 1:
                CottiSlope  = (float(CottiSlope_Eingabe.get()))

                
                
    def Next():
        
        if BaseCorr == 0:
            KL_NextLevel4()
        if BaseCorr == 1:
            BaseGetter()
            KL_NextLevel4()
        
  
        
    def Back():
        KL_NextLevel2()
        
        def quit():
            Fenster.destroy()
        quit()   
    
    
    
    
    
    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=17,column=1)
    button=Button(Fenster,text="Next",command=Next).grid(row=18 ,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=19,column=1)
    


# In[ ]:

def KL_NextLevel4():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Koutecky-Levich Analysis Level 4")                         
    Fenster.geometry("1200x600")
    
    #DIESE ZWEI FUNKTIONEN MÜSSEN HIER DRIN DEFINIERT WERDEN!!!
        
    def GeradenFit(xWert,ySchnitt,Steig):
        return  ySchnitt + Steig*xWert

    def EchtesPotFinden(WoSieSuchenSoll,WasSieSuchenSoll):                         
        idxStart = (np.abs(WoSieSuchenSoll-WasSieSuchenSoll)).argmin()               
        return WoSieSuchenSoll[idxStart]           
    
    
    f = Figure(figsize=(12, 6), dpi=100)
    #f.subplots_adjust(hspace=0.4)
    
    if OnlyLSVPlot == 1:
        LSVFenster = Toplevel()                                                         
        LSVFenster.title("LSV-Plot")                         
        LSVFenster.geometry("700x700")
        
        LSVPlot = Figure(figsize=(6, 6), dpi=100)
        
    
    global Length
    global ShortLength
    
    global ILimMittel 
    global InvILimMittel
    global InvILimFITTED
    
    #um später den Offset rausholen zu können
    global KL_Offset
    
    Length = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:])))*REFF)
    ShortLength = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:]))))
    
    global KORRPOTENZIALSUPPERARRAY
    global KORRSTROMSUPPERARRAY
    global UNKORRPOTENZIALSUPPERARRAY
    global UNKORRSTROMSUPPERARRAY
    
    KORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,Length])
    KORRSTROMSUPPERARRAY     = np.empty([NumMeas,Length])
    
    UNKORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,ShortLength])
    UNKORRSTROMSUPPERARRAY     = np.empty([NumMeas,ShortLength])
    
    
    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #FOR SCHLEIFE WENN ROT RATES GEÄNDERT WERDEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    
    if VarParRot == 1:
    
        OCPs = OCPArray
        
        #________________________________________________________________      
        
              
        ILimMittel    = np.empty(int(NumMeas))
        

        for i in range(int(NumMeas)):
            StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
            UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            
            
            if DesiReactxx == 0:
                Stromarrays                = StromarraysROH[::] 
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::] 
                Potenzialarrays            = PotenzialarraysROH[::]   -  Ru*Stromarrays/1000000
            
            
            if DesiReactxx == 1:
                Stromarrays                = StromarraysROH
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
                Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
            
            
            
            if BaseCorr == 1:
                BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
                if DesiReactxx == 0:
                    BaseStromArrays        = BaseStromArraysROH[::]
                if DesiReactxx == 1:
                    BaseStromArrays        = BaseStromArraysROH[::]
                Stromarrays            = Stromarrays - BaseStromArrays
        
            LenPotArrays            = len(Potenzialarrays)
        
            AufgelPotArrays        = np.linspace(Potenzialarrays[0], Potenzialarrays[-1], num=int(LenPotArrays*REFF))#, endpoint=True)
            Interpolation          = interp1d(Potenzialarrays, Stromarrays, kind='cubic')                  
            AuflgelStrarrays       = Interpolation(AufgelPotArrays)                                        
        
            if Smooth_Data == 1:
            
                from scipy.interpolate import UnivariateSpline
            
                SmoothFromPotential = EchtesPotFinden(Potenzialarrays,SmoothFrom)
                SmoothToPotential   = EchtesPotFinden(Potenzialarrays,SmoothTo)
            
                IdxSmoothFromPotential = np.asscalar(np.where(Potenzialarrays == SmoothFromPotential) [0])
                IdxSmoothToPotential   = np.asscalar(np.where(Potenzialarrays == SmoothToPotential) [0])
            
            
                spl = UnivariateSpline(Potenzialarrays, Stromarrays)
                spl.set_smoothing_factor(SmoothFact)
                AufgelPotArrays  = np.linspace(Potenzialarrays[IdxSmoothFromPotential],Potenzialarrays[IdxSmoothToPotential] , int(REFF*len(Potenzialarrays)))
                AuflgelStrarrays = spl(AufgelPotArrays)
            
            
            #___________________________________________________________
            #AUFGELÖSTES POTENZIALARRAY UM OHMSCHEN ANTEIL KORRIGIEREN
            #___________________________________________________________
            KorrAufgPotArr         = AufgelPotArrays #Widerstandskorrektur erfolgt schon weiter oben!!!
        
            KorrAufgStrArrays      = AuflgelStrarrays  
            if AutoCorr == 1:
                
                OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
                OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
                KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
           
            if ManuCorr == 1:
                KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr 
        
   
            KLStartPotenzial       = EchtesPotFinden(KorrAufgPotArr,KLStart)
            KLEndPotenzial         = EchtesPotFinden(KorrAufgPotArr,KLEnd)
            
            KLStartIndex           = np.asscalar(np.where(KorrAufgPotArr == KLStartPotenzial) [0])
            KLEndIndex             = np.asscalar(np.where(KorrAufgPotArr == KLEndPotenzial) [0])
    

            
            if KLStartIndex < KLEndIndex :
                Strommittelwerte       = np.mean(KorrAufgStrArrays[KLStartIndex:KLEndIndex])
            if KLStartIndex > KLEndIndex :
                Strommittelwerte       = np.mean(KorrAufgStrArrays[KLEndIndex:KLStartIndex])
    
            ILimMittel[i] = Strommittelwerte
            
            
            for j in range(Length):
                
                KORRPOTENZIALSUPPERARRAY[i,j] = KorrAufgPotArr[j] 
                KORRSTROMSUPPERARRAY[i,j]     = KorrAufgStrArrays[j]
                      
            for j in range(ShortLength):
                
                UNKORRPOTENZIALSUPPERARRAY[i,j] = UnkorrPotenzialarrays[j] 
                UNKORRSTROMSUPPERARRAY[i,j]     = StromarraysROH[j]
            
            
        
            
            
        
            bild1 = f.add_subplot(121)
            if ShowUncorrLSV == 1:
                bild1.plot(UnkorrPotenzialarrays,StromarraysROH, linestyle='-',marker='',color='k')
            bild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
            bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            bild1.set_ylabel('I [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild1.spines[axis].set_linewidth(2)
            bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            bild1.axvline(OCPs[i],color='k')
            bild1.axhline(0, color='k')
            bild1.axvline(KLStart, color='r', linestyle='--')
            bild1.axvline(KLEnd, color='r', linestyle='--')
            
            
            if OnlyLSVPlot == 1:
                LSVbild1 = LSVPlot.add_subplot(111)
                if ShowUncorrLSV == 1:
                    LSVbild1.plot(UnkorrPotenzialarrays,StromarraysROH, linestyle='-',marker='',color='k')
                LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
                LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
                LSVbild1.set_ylabel('I [$\mu$A]', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    LSVbild1.spines[axis].set_linewidth(2)
                LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                LSVbild1.axvline(OCPs[i],color='k')
                LSVbild1.axhline(0, color='k')
                LSVbild1.axvline(KLStart, color='r', linestyle='--')
                LSVbild1.axvline(KLEnd, color='r', linestyle='--')


        global InvWurzelRot
        
        
        InvWurzelRot    = 1/(RotRatesArray**0.5)
        InvILimMittel   = 1/ILimMittel
        
     
        KLFitting, pcov = curve_fit(GeradenFit, InvWurzelRot,InvILimMittel)   
        
        KL_Offset = KLFitting[0]
        
        
        InvILimFITTED = np.empty(int(NumMeas))
        for i in range(int(NumMeas)):
            InvILimFITTED[i] = GeradenFit(InvWurzelRot[i], *KLFitting)
        
        
        bild2 = f.add_subplot(122)
        
        def KLROTPLOTTER():
        
            bild2.plot(InvWurzelRot,GeradenFit(InvWurzelRot, *KLFitting),color='r',linestyle='-',marker='')
            bild2.plot(InvWurzelRot,InvILimMittel,linestyle='',marker='.')
            bild2.set_xlabel(r'$\sqrt{\omega^{-1}}$' '[min$^{-0.5}$]', fontsize=12)
            bild2.set_ylabel('I$^{-1}$ [$\mu$A$^{-1}$]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild2.spines[axis].set_linewidth(2)
            bild2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
            #________________________________________________________
            #________________________________________________________
            #ANNOTOATIONS FOR PARAMETERS AND CALCULATION OF n D or c
            #________________________________________________________
            #________________________________________________________
        
    
        
            #IF SCHLEIFE WENN n BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
            
            global n_KL
            global D_KL
            global c_KL
            
            if get_n ==1:
      
                n_KL = np.absolute(1/(-0.201*1000000*F*A*(D**0.66666666)*Visc**(-0.16666666)*c*KLFitting[1]))

                if RefToCotti == 1:
                    
                    BeschreibefensterCotti = Toplevel()                                                         
                    BeschreibefensterCotti.title("Information")                         
                    BeschreibefensterCotti.geometry("200x200")    
                    T = Text(BeschreibefensterCotti, height=6, width=30)
                    T.pack()
                    T.insert(END, "D got referred to Cottrell")
                    
                    Di = (np.absolute(0.201*CottiSlope*KLFitting[1]*np.pi**0.5 * Visc**(-0.16666666)))**(-6)
                
                    n_KL = np.absolute(1/(-0.201*1000000*F*A*(Di**0.66666666)*Visc**(-0.16666666)*c*KLFitting[1]))
                      
    
                
                bild2.annotate('n     =%8.3f' % n_KL, xy=(0.5, 0.93), xycoords='axes fraction')           
        
        
            #IF SCHLEIFE WENN D BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_D ==1:
            
                D_KL = 1000000*(np.absolute(1/(-0.201*1000000*F*A*n*Visc**(-0.16666666)*c*KLFitting[1])))**1.5
            
                
                bild2.annotate('D =%5.3f*10$^{-6} cm^{2}/s$' % D_KL, xy=(0.5, 0.93), xycoords='axes fraction') 
        
            #IF SCHLEIFE WENN c BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_c ==1:
            
                c_KL = 1000000*np.absolute(1/(-0.201*1000000*F*A*(D**0.66666666)*Visc**(-0.16666666)*n*KLFitting[1]))
            
                
                bild2.annotate('c =%5.3f*10$^{-6} mol/cm^{3}$' % c_KL, xy=(0.5, 0.93), xycoords='axes fraction') 
        
        KLROTPLOTTER()

    
        canvas = FigureCanvasTkAgg(f, master=Fenster)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, Fenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        if OnlyLSVPlot ==1:
            canvas = FigureCanvasTkAgg(LSVPlot, master=LSVFenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, LSVFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        

 
        
        

        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #ONLY KL PLOT OPTION
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        
        
        if OnlyKLPlot == 1:
            
            KLFenster = Toplevel()                                                         
            KLFenster.title("Koutecky-Levich Plot")                         
            KLFenster.geometry("700x700")
        
            KLPlot = Figure(figsize=(6, 6), dpi=100)
            
            bild2 = KLPlot.add_subplot(111)
            
            KLROTPLOTTER()
                
            canvas = FigureCanvasTkAgg(KLPlot, master=KLFenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, KLFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        if AsTxtSaver == 1:
            KLROTASTXTsaver()

    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #FOR SCHLEIFE WENN CONCENTRATIONS GEÄNDERT WERDEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________ 
    if VarParCon == 1:
    
        OCPs = OCPArray
        
        
        
        #________________________________________________________________      
        
               
        ILimMittel    = np.empty(int(NumMeas))
        

        for i in range(int(NumMeas)):
            StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
            UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            
            
            if DesiReactxx == 0:
                Stromarrays                = StromarraysROH[::]
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::]
                Potenzialarrays            = PotenzialarraysROH[::]   -  Ru*Stromarrays/1000000
            
            
            if DesiReactxx == 1:
                Stromarrays                = StromarraysROH
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
                Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
            
            
            
            if BaseCorr == 1:
                BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
                if DesiReactxx == 0:
                    BaseStromArrays        = BaseStromArraysROH[::]
                if DesiReactxx == 1:
                    BaseStromArrays        = BaseStromArraysROH[::]
                Stromarrays            = Stromarrays - BaseStromArrays
            
            LenPotArrays            = len(Potenzialarrays)
        
            AufgelPotArrays        = np.linspace(Potenzialarrays[0], Potenzialarrays[-1], num=int(LenPotArrays*REFF))#, endpoint=True)
            Interpolation          = interp1d(Potenzialarrays, Stromarrays, kind='cubic')                  
            AuflgelStrarrays       = Interpolation(AufgelPotArrays)                                        
        
            if Smooth_Data == 1:
            
                from scipy.interpolate import UnivariateSpline
            
                SmoothFromPotential = EchtesPotFinden(Potenzialarrays,SmoothFrom)
                SmoothToPotential   = EchtesPotFinden(Potenzialarrays,SmoothTo)
            
                IdxSmoothFromPotential = np.asscalar(np.where(Potenzialarrays == SmoothFromPotential) [0])
                IdxSmoothToPotential   = np.asscalar(np.where(Potenzialarrays == SmoothToPotential) [0])
            
            
                spl = UnivariateSpline(Potenzialarrays, Stromarrays)
                spl.set_smoothing_factor(SmoothFact)
                AufgelPotArrays  = np.linspace(Potenzialarrays[IdxSmoothFromPotential],Potenzialarrays[IdxSmoothToPotential] , int(REFF*len(Potenzialarrays)))
                AuflgelStrarrays = spl(AufgelPotArrays)
            
            
            #___________________________________________________________
            #AUFGELÖSTES POTENZIALARRAY UM OHMSCHEN ANTEIL KORRIGIEREN
            #___________________________________________________________
            KorrAufgPotArr         = AufgelPotArrays # Korrektur erfolgte schon weiter oben
        
            KLStartPotenzial       = EchtesPotFinden(KorrAufgPotArr,KLStart)
            KLEndPotenzial         = EchtesPotFinden(KorrAufgPotArr,KLEnd)
            
            KLStartIndex           = np.asscalar(np.where(KorrAufgPotArr == KLStartPotenzial) [0])
            KLEndIndex             = np.asscalar(np.where(KorrAufgPotArr == KLEndPotenzial) [0])
    
            KorrAufgStrArrays      = AuflgelStrarrays  
            if AutoCorr == 1:
                
                OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
                OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
                KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
                #print OCPotenziale
            
            if ManuCorr == 1:
                KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr    
                
                

            
            if KLStartIndex < KLEndIndex :
                Strommittelwerte       = np.mean(KorrAufgStrArrays[KLStartIndex:KLEndIndex])
            if KLStartIndex > KLEndIndex :
                Strommittelwerte       = np.mean(KorrAufgStrArrays[KLEndIndex:KLStartIndex])
    
            ILimMittel[i] = Strommittelwerte
        
        
        
            for j in range(Length):
                
                KORRPOTENZIALSUPPERARRAY[i,j] = KorrAufgPotArr[j] 
                KORRSTROMSUPPERARRAY[i,j]     = KorrAufgStrArrays[j]
                      
            for j in range(ShortLength):
                
                UNKORRPOTENZIALSUPPERARRAY[i,j] = UnkorrPotenzialarrays[j] 
                UNKORRSTROMSUPPERARRAY[i,j]     = StromarraysROH[j]
            
        
            bild1 = f.add_subplot(121)
            if ShowUncorrLSV == 1:
                bild1.plot(UnkorrPotenzialarrays,StromarraysROH, linestyle='-',marker='',color='k')
            bild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
            bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            bild1.set_ylabel('I [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild1.spines[axis].set_linewidth(2)
            bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            bild1.axvline(OCPs[i],color='k')
            bild1.axhline(0, color='k')
            bild1.axvline(KLStart, color='r', linestyle='--')
            bild1.axvline(KLEnd, color='r', linestyle='--')
        
        
            if OnlyLSVPlot == 1:
                LSVbild1 = LSVPlot.add_subplot(111)
                if ShowUncorrLSV == 1:
                    LSVbild1.plot(UnkorrPotenzialarrays,StromarraysROH, linestyle='-',marker='',color='k')
                LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
                LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
                LSVbild1.set_ylabel('I [$\mu$A]', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    LSVbild1.spines[axis].set_linewidth(2)
                LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                LSVbild1.axvline(OCPs[i],color='k')
                LSVbild1.axhline(0, color='k')
                LSVbild1.axvline(KLStart, color='r', linestyle='--')
                LSVbild1.axvline(KLEnd, color='r', linestyle='--')
            
        
        
        global InvConcArray
            
        InvConcArray    = 1/ConcentrationsArray
        InvILimMittel   = 1/ILimMittel
        
     
        KLFitting, pcov = curve_fit(GeradenFit, InvConcArray,InvILimMittel)   

        KL_Offset = KLFitting[0]
        

        InvILimFITTED = np.empty(int(NumMeas))
        for i in range(int(NumMeas)):
            InvILimFITTED[i] = GeradenFit(InvConcArray[i], *KLFitting)
        
        bild2 = f.add_subplot(122)
        
        def KLCONCPLOTTER():
            bild2.plot(InvConcArray,GeradenFit(InvConcArray, *KLFitting),color='r',linestyle='-',marker='')
            bild2.plot(InvConcArray,InvILimMittel,linestyle='',marker='.')
            bild2.set_xlabel('$c^{-1}}$' '[cm$^{3}$ mol$^{-1}$]', fontsize=12)
            bild2.set_ylabel('I$^{-1}$ [$\mu$A$^{-1}$]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild2.spines[axis].set_linewidth(2)
            bild2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
            #________________________________________________________
            #________________________________________________________
            #ANNOTATIONS FOR PARAMETERS AND CALCULATION OF n D or c
            #________________________________________________________
            #________________________________________________________
        
        
            #IF SCHLEIFE WENN n BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
            global n_KL
            global D_KL
            global c_KL
            
            if get_n ==1:
         
                n_KL = np.absolute(1/(-0.201*1000000*F*A*(D**0.66666666)*Visc**(-0.16666666)*c *RotRate**0.5 *KLFitting[1]))

                if RefToCotti == 1:
               
                    BeschreibefensterCotti = Toplevel()                                                         
                    BeschreibefensterCotti.title("Information")                         
                    BeschreibefensterCotti.geometry("200x200")    
                    T = Text(BeschreibefensterCotti, height=6, width=30)
                    T.pack()
                    T.insert(END, "D got referred to Cottrell")
                    
                    Di = (np.absolute(0.201*CottiSlope*KLFitting[1]*np.pi**0.5 * Visc**(-0.16666666)))**(-6)
                
                    n_KL = np.absolute(1/(-0.201*1000000*F*A*(Di**0.66666666)*Visc**(-0.16666666) *c *RotRate**0.5 *KLFitting[1]))
      
                bild2.annotate('n     =%8.3f' % n_KL, xy=(0.5, 0.93), xycoords='axes fraction')           
        
        
            #IF SCHLEIFE WENN D BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_D ==1:
            
                D_KL = 1000000*(np.absolute(1/(-0.201*1000000*F*A*n*Visc**(-0.16666666)*c *RotRate**0.5 *KLFitting[1])))**1.5
            
                bild2.annotate('D =%5.3f*10$^{-6} cm^{2}/s$' % D_KL, xy=(0.5, 0.93), xycoords='axes fraction') 
        
            #IF SCHLEIFE WENN c BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_c ==1:
            
                bild2.annotate('GET c NOT POSSIBLE HERE', xy=(0.5, 0.93), xycoords='axes fraction') 
        

        KLCONCPLOTTER()
        
        
        canvas = FigureCanvasTkAgg(f, master=Fenster)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, Fenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        if OnlyLSVPlot ==1:
            canvas = FigureCanvasTkAgg(LSVPlot, master=LSVFenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, LSVFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        

        
                
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #ONLY KL PLOT OPTION
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        
        
        if OnlyKLPlot == 1:
            
            KLFenster = Toplevel()                                                         
            KLFenster.title("Koutecky-Levich Plot")                         
            KLFenster.geometry("700x700")
        
            KLPlot = Figure(figsize=(6, 6), dpi=100)
            
            bild2 = KLPlot.add_subplot(111)
            
            KLCONCPLOTTER()
                
            canvas = FigureCanvasTkAgg(KLPlot, master=KLFenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, KLFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        if AsTxtSaver == 1:
            KLCONASTXTsaver()
        
            




# In[ ]:

def KLROTASTXTsaver():

 
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")
        if get_n == 1:
            f.write("n-KL")
            f.write("\t")
            f.write(str(np.asscalar(n_KL)))
            f.write("\n")
            f.write("KL_Offset")
            f.write("\t")
            f.write(str(np.asscalar(KL_Offset)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("D in cm^2/s")
            f.write("\t")
            f.write(str(D))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Visc in cm^2/s")
            f.write("\t")
            f.write(str(Visc))
            f.write("\n")
            f.write("c in mol/cm^3")
            f.write("\t")
            f.write(str(c))
            f.write("\n")
            f.write("KL_StartPot in V")
            f.write("\t")
            f.write(str(KLStart))
            f.write("\n")
            f.write("KL_EndPot in V")
            f.write("\t")
            f.write(str(KLEnd))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base")
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
            if RefToCotti  ==1:
                f.write("D ref. to Cottrell")
                f.write("\n")
                
            
  

            
        if get_c == 1:
            f.write("c-KL")
            f.write("\t")
            f.write(str(0.000001*np.asscalar(c_KL)))
            f.write("\n")
            f.write("KL_Offset")
            f.write("\t")
            f.write(str(np.asscalar(KL_Offset)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
            f.write("\n")
            f.write("D in cm^2/s")
            f.write("\t")
            f.write(str(D))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Visc in cm^2/s")
            f.write("\t")
            f.write(str(Visc))
            f.write("\n")
            f.write("KL_StartPot in V")
            f.write("\t")
            f.write(str(KLStart))
            f.write("\n")
            f.write("KL_EndPot in V")
            f.write("\t")
            f.write(str(KLEnd))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base") 
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
                
                
                
        if get_D == 1:
            f.write("D-KL")  
            f.write("\t")
            f.write(str(0.000001*np.asscalar(D_KL)))
            f.write("\n")
            f.write("KL_Offset")
            f.write("\t")
            f.write(str(np.asscalar(KL_Offset)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Visc in cm^2/s")
            f.write("\t")
            f.write(str(Visc))
            f.write("\n")
            f.write("c in mol/cm^3")
            f.write("\t")
            f.write(str(c))
            f.write("\n")
            f.write("KL_StartPot in V")
            f.write("\t")
            f.write(str(KLStart))
            f.write("\n")
            f.write("KL_EndPot in V")
            f.write("\t")
            f.write(str(KLEnd))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base") 
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
        
        
        
            
        
        
        for i in range(3):
            f.write("\n")
        
        f.write("KL-Plot Data")
        f.write("\n")
        f.write("\n")
        f.write("OCPs in V")   
        f.write("\t")
        f.write("RotRates in rpm")   
        f.write("\t")
        f.write("I-Lim in MicAmp")   
        f.write("\t")
        f.write("Inv root RotRates")   
        f.write("\t")
        f.write("Inv I Lim ")   
        f.write("\t")
        f.write("Fit Inv I Lim")   
        f.write("\t")
        f.write("\n")
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(RotRatesArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(ILimMittel[i])))
            f.write("\t")
            f.write(str(np.asscalar(InvWurzelRot[i])))
            f.write("\t")
            f.write(str(np.asscalar(InvILimMittel[i])))
            f.write("\t")
            f.write(str(np.asscalar(InvILimFITTED[i])))
            f.write("\t")
            f.write("\n")
            

        
        for i in range(3):
            f.write("\n")
#_________________________________________________________________________________            
        f.write("Korrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Rot-Speed")
        f.write("\t")
    
        for jj in range(len(RotRatesArray)):
        
            f.write(str(np.asscalar(RotRatesArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(len(RotRatesArray)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(Length)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(len(RotRatesArray)):
                b = str(np.asscalar(KORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(KORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")
#__________________________________________________________________________________________________-        
        for i in range(3):
            f.write("\n")
            
        f.write("Uncorrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Rot-Speed")
        f.write("\t")
    
        for jj in range(len(RotRatesArray)):
        
            f.write(str(np.asscalar(RotRatesArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(len(RotRatesArray)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(ShortLength)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(len(RotRatesArray)):
                b = str(np.asscalar(UNKORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(UNKORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")

    root.destroy()




# In[ ]:

def KLCONASTXTsaver():

 
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")
        if get_n == 1:
            f.write("n-KL")
            f.write("\t")
            f.write(str(np.asscalar(n_KL)))
            f.write("\n")
            f.write("KL_Offset")
            f.write("\t")
            f.write(str(np.asscalar(KL_Offset)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("D in cm^2/s")
            f.write("\t")
            f.write(str(D))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Visc in cm^2/s")
            f.write("\t")
            f.write(str(Visc))
            f.write("\n")
            f.write("Rotation speed in rpm")
            f.write("\t")
            f.write(str(RotRate))
            f.write("\n")
            f.write("KL_StartPot in V")
            f.write("\t")
            f.write(str(KLStart))
            f.write("\n")
            f.write("KL_EndPot in V")
            f.write("\t")
            f.write(str(KLEnd))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base")
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
            if RefToCotti  ==1:
                f.write("D ref. to Cottrell")
                f.write("\n")
                
            
  

            
        if get_c == 1:
            f.write("c-KL")
            f.write("\t")
            f.write("Getting c not possible here")
            f.write("\n")
            f.write("KL_Offset")
            f.write("\t")
            f.write(str(np.asscalar(KL_Offset)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
            f.write("\n")
            f.write("Rotation speed in rpm")
            f.write("\t")
            f.write(str(RotRate))
            f.write("\n")
            f.write("D in cm^2/s")
            f.write("\t")
            f.write(str(D))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Visc in cm^2/s")
            f.write("\t")
            f.write(str(Visc))
            f.write("\n")
            f.write("KL_StartPot in V")
            f.write("\t")
            f.write(str(KLStart))
            f.write("\n")
            f.write("KL_EndPot in V")
            f.write("\t")
            f.write(str(KLEnd))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base") 
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
                
                
                
        if get_D == 1:
            f.write("D-KL")  
            f.write("\t")
            f.write(str(0.000001*np.asscalar(D_KL)))
            f.write("\n")
            f.write("KL_Offset")
            f.write("\t")
            f.write(str(np.asscalar(KL_Offset)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
            f.write("\n")
            f.write("Rotation speed in rpm")
            f.write("\t")
            f.write(str(RotRate))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Visc in cm^2/s")
            f.write("\t")
            f.write(str(Visc))
            f.write("\n")
            f.write("KL_StartPot in V")
            f.write("\t")
            f.write(str(KLStart))
            f.write("\n")
            f.write("KL_EndPot in V")
            f.write("\t")
            f.write(str(KLEnd))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base") 
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")

            
        
        
        for i in range(3):
            f.write("\n")
        
        f.write("KL-Plot Data")
        f.write("\n")
        f.write("\n")
        f.write("OCPs in V")   
        f.write("\t")
        f.write("Concentrations in mol/cm^3")   
        f.write("\t")
        f.write("I-Lim in MicAmp")   
        f.write("\t")
        f.write("Inv Concentrations")   
        f.write("\t")
        f.write("Inv I Lim ")   
        f.write("\t")
        f.write("Fit Inv I Lim")   
        f.write("\t")
        f.write("\n")
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(c*ConcentrationsArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(ILimMittel[i])))
            f.write("\t")
            f.write(str(np.asscalar(1/c*ConcentrationsArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(InvILimMittel[i])))
            f.write("\t")
            f.write(str(np.asscalar(InvILimFITTED[i])))
            f.write("\t")
            f.write("\n")
            

        
        for i in range(3):
            f.write("\n")
#_________________________________________________________________________________            
        f.write("Korrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Concentrations as Fract of highest")
        f.write("\t")
    
        for jj in range(len(ConcentrationsArray)):
        
            f.write(str(np.asscalar(ConcentrationsArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(len(ConcentrationsArray)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(Length)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(len(ConcentrationsArray)):
                b = str(np.asscalar(KORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(KORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")
#__________________________________________________________________________________________________-        
        for i in range(3):
            f.write("\n")
            
        f.write("Uncorrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Concentrations as Fract of highest")
        f.write("\t")
    
        for jj in range(len(ConcentrationsArray)):
        
            f.write(str(np.asscalar(ConcentrationsArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(len(ConcentrationsArray)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(ShortLength)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(len(ConcentrationsArray)):
                b = str(np.asscalar(UNKORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(UNKORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")

    root.destroy()



