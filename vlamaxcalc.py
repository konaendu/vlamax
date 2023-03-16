import numpy as np
import matplotlib.pyplot as plt
from tkinter import *
from tkinter.filedialog import askopenfilename
import fitparse
from scipy.integrate import solve_ivp

# Define Functions
def plotLactate(VO2max, vLamax, Eco, Ks1, Ks2, VolRel, bw):
    # generate oxygen steady-states
    VO2ss = np.arange(1, VO2max, 0.01)

    # potentiate constants
    Ks1 = Ks1 ** 2
    Ks2 = Ks2 ** 3

    # calculate ADP corresponding to VO2ss (eq. 4b)
    ADP = np.sqrt((Ks1 * VO2ss) / (VO2max - VO2ss))

    # calculate steady state gross lactic acid (pyruvate) formation rate (eq. 3)
    vLass = 60 * vLamax / (1 + (Ks2 / ADP ** 3))

    # calculate lactate combustion (eq. 6)
    LaComb = (0.01576 / VolRel) * VO2ss

    # calculate net lactate formation
    vLanet = abs(vLass - LaComb)

    # Calculate overall demand in O2-equivalents and Intensity
    overall_demand = ((vLass * (VolRel * bw) * (5.5) / bw) + VO2ss) # (1 / 4.3) * 22.4) #4.3=4.0727272 in inscyd
    Intensity = overall_demand / Eco

    # calculate the crossing point (AT) between gross lactate production and combustion
    arg_sAT = np.argmin(vLanet)
    sAT = Intensity[arg_sAT]

    # speed at FatMax
    s_Fmx = Intensity[np.argmax(vLanet[0:arg_sAT])]

    # percentage of VO2max
    pcVO2maxAT = sAT * Eco / VO2max

    if Eco > 5:
        plt.plot(Intensity, vLanet, label = 'Lack of Pyruvate & Lactate Accumulation')
        plt.xlabel('Speed [m\s]')
        plt.ylabel('mmol/L/min')
        plt.title(f'AT is at {round(sAT, 2)} m/s, {round(pcVO2maxAT*100, 2)} % of VO2max, Fatmax is at {round(s_Fmx, 2)} m/s')
        plt.legend()
        plt.ylim(0,8)
    else:
        plt.plot(Intensity, vLanet, label='Lack of Pyruvate & Lactate Accumulation')
        plt.xlabel('Power [W]')
        plt.ylabel('mmol/L/min')
        plt.title(f'AT is at {round(sAT, 2)} W, {round(pcVO2maxAT * 100, 2)} % of VO2max, Fatmax is at {round(s_Fmx, 2)} W')
        plt.legend()
        plt.ylim(0, 8)
    plt.show()


def plotClass(VO2max, VLamax, Eco, Ks, Kss, VolRel, Kel, bw):
    Kel = Kel
    Kss = Kss ** 3
    Ks = Ks ** 2
    VO2ss = np.arange(1, VO2max, 0.01)
    ADP = np.sqrt((Ks * VO2ss) / (VO2max - VO2ss))
    vLass = 60 * VLamax / (1 + (Kss / ADP ** 3))
    overall_demand = ((vLass * (VolRel * bw) * (5.5) / bw) + VO2ss)
    Intensity = overall_demand / Eco
    Class = np.sqrt((VLamax * Kel * 60) / (((0.01576/ VolRel) * VO2ss) * (1 + (Kss / ((Ks * VO2ss) / (VO2max - VO2ss)) ** (3 / 2))) - (VLamax * 60)))

    if Eco > 5:
        plt.plot(Intensity, Class, label = 'Steady-State Lactate')
        plt.xlabel('Speed [m\s]')
        plt.ylabel('mmol/L/min')
        plt.legend()
        plt.ylim(0,8)
    else:
        plt.plot(Intensity, Class, label='Steady-State Lactate')
        plt.xlabel('Power [W]')
        plt.ylabel('mmol/L/min')
        plt.legend()
        plt.ylim(0, 8)
    plt.show()

def Macronutrients(VO2max, vLamax, Eco, Ks1, Ks2, VolRel, bw):
    # generate oxygen steady-states
    VO2ss = np.arange(1, VO2max, 0.01)

    # potentiate constants
    Ks1 = Ks1 ** 2
    Ks2 = Ks2 ** 3

    # calculate ADP corresponding to VO2ss (eq. 4b)
    ADP = np.sqrt((Ks1 * VO2ss) / (VO2max - VO2ss))

    # calculate steady state gross lactic acid (pyruvate) formation rate (eq. 3)
    vLass = 60 * vLamax / (1 + (Ks2 / ADP ** 3))

    # calculate lactate combustion (eq. 6)
    LaComb = (0.01576 / VolRel) * VO2ss

    # calculate net lactate formation
    vLanet = abs(vLass - LaComb)

    # Calculate overall demand in O2-equivalents and Intensity
    overall_demand = ((vLass * (VolRel * bw) * (5.5) / bw) + VO2ss)
    Intensity = overall_demand / Eco

    # calculate the crossing point (AT) between gross lactate production and combustion
    arg_sAT = np.argmin(vLanet)
    sAT = Intensity[arg_sAT]

    CHO_util = vLass * (bw * VolRel) * 60 / 1000 / 2 * 162.14                       # 162.14 molare masse glykogen
    Fat_util = (vLanet[:arg_sAT] * VolRel) / 0.01576 * bw * 60 * 4.65 / 9.5 / 1000  # 4.65 kcal/L VO2 stearinsäure, 9.5 g/kcal

    if Eco > 5:
        fig, ax = plt.subplots()
        ax.plot(Intensity, CHO_util, label = 'Carb.', color='darkgoldenrod')
        ax.set_xlabel('Speed [m\s]')
        ax.set_ylabel('Carb. g/h')
        ax.legend(loc='upper left')
        ax.set_ylim(0,400)
        ax1 = ax.twinx()
        ax1.plot(Intensity[:arg_sAT], Fat_util, label = 'Fat', color='green')
        ax1.set_ylabel('Fat g/h')
        ax1.set_ylim(0, 100)
        ax1.legend()

    else:
        fig, ax = plt.subplots()
        ax.plot(Intensity, CHO_util, label='Carb.', color='darkgoldenrod')
        ax.set_xlabel('Power [W]')
        ax.set_ylabel('Carb. g/h')
        ax.legend(loc='upper left')
        ax.set_ylim(0, 400)
        ax1 = ax.twinx()
        ax1.plot(Intensity[:arg_sAT], Fat_util, label='Fat', color='green')
        ax1.set_ylabel('Fat g/h')
        ax1.set_ylim(0, 100)
        ax1.legend()
    plt.show()

def get_values():
    VO2max = eval(entry_VO2max.get())
    VLamax = eval(entry_VLamax.get())
    Volrel = eval(entry_Volrel.get())
    Eco = eval(entry_Eco.get())
    bw = eval(entry_Bodyweight.get())
    if entry_KS1.get() == 'Ks1':
        Ks1 = 0.25
    else:
        Ks1 = eval(entry_KS1.get())

    if entry_KS2.get() == 'Ks2':
        Ks2 = 1.1
    else:
        Ks2 = eval(entry_KS2.get())

    if entry_KEL.get() == 'Kel':
        Kel = 2.5
    else:
        Kel = eval(entry_KEL.get())

    return VO2max, VLamax, Eco, Volrel, bw, Ks1, Ks2, Kel

def Adenosine_Phosphates(CHEP_System, ATP, VO2a, VO2max, Pi, Lam, Sa=7, Sc=25, M3=0.96):
    # calculating hydrogen ion concentration
    pCO2 = 40 + 55 * (VO2a / VO2max)
    dbuff = 1 / 54
    pHm = 7.85 + dbuff * (0.8 * Pi - Lam) - 0.55 * np.log10(pCO2)
    H = 10 ** (-pHm)

    # Calculating Phosphates
    PCr = CHEP_System - ATP
    M1 = H * 1.66 * 10 ** 9
    Q = M1 * (PCr / (Sc - PCr))
    ADP = Sa * Q / (M3 + Q + Q ** 2)
    ATP = ADP * Q
    #AMP = Sa - ADP - ATP
    AMP = ADP ** 2
    Pi = Sc - PCr

    return PCr, Pi, ATP, ADP, AMP, M1, H, pHm

def Model(t, y, ATP_Belastung, VO2max, vlamax, Body_weight, Musclemass, Muscle_Aktiv, Ks1, Ks2, Ks3, Kelox, KpyrO2, T,
          VolRel_m, VolRel_b, Sc, Sa, M3, ATP_Ruhe, Vrel, PCr, Pi, ATP, ADP, AMP, M1, H, kdiff):

    GP, VO2a, Lam, Lab = y

    #Additional parameters of funktion
    #vlaoxm = ((KpyrO2 * VO2a) / (1 + ((VolRel_m ** 2 * Kelox) / Lam ** 2)))
    #vlaoxb = (KpyrO2 * VO2a) / (1 + ((Vrel ** 2 * Kelox) / (Lab ** 2)))
    #vlaoxm = ((0.021 * VO2a) / (1 + (Kelox/Lam**2))) * (2/3)
    #vlaoxb = ((0.021 * VO2a) / (1 + (Kelox/Lam**2))) * (1/3)
    # ADP_Gly = np.sqrt(abs((Ks1*VO2a)/(VO2max-VO2a)))

    #VO2ss = VO2max / (1 + (Ks1 / ADP ** 2))
    # gl = (1 / (1 + (H ** 3 / Ks3))) * (vlamax / (1 + (Ks2 / (ADP_Gly**3))))
    # gl = (1 / (1 + (H ** 3 / Ks3))) * (vlamax / (1 + (Ks2 / (ADP * AMP))))
    # gl = (1 / (1 + (H ** 3/ Ks3))) * (vlamax / (1 + Ks2 * (((VO2max - VO2a)/(Ks1 * VO2a)) ** 1.5)))

    #### System of ODEs

    # 1. ODE
    dGP = VO2a * 0.233 + (
                (1 / (1 + (H ** 3 / Ks3))) * (vlamax / (1 + (Ks2 / (ADP * AMP))))) * 1.4 - ATP_Belastung - ATP_Ruhe

    # 2. ODE
    dVO2a = ((VO2max / (1 + (Ks1 / ADP ** 2))) - VO2a) / T

    # 3. ODE
    dLam = ((VolRel_m ** -1) * (((1 / (1 + (H ** 3 / Ks3))) * (vlamax / (1 + (Ks2 / (ADP * AMP))))) - (
                (0.021 * VO2a) / (1 + (Kelox / Lam ** 2))) * (2 / 3))) - kdiff * (Lam - Lab)

    # 4. ODE
    dLab = Vrel * (kdiff * (Lam - Lab) - ((0.021 * VO2a) / (1 + (Kelox / Lam ** 2))) * (1 / 3))
    return np.array([dGP, dVO2a, dLam, dLab])


def call_Mader(ATP_Belastung=np.repeat(250, 30*60), y0s=np.array([25, 350 / 60 / 18.2, 1.1, 1.1]),
               VO2max=150, vlamax= 60, Body_weight=70, Musclemass=100, Muscle_Aktiv=0.26,
               Ks1=0.035 ** 2, Ks2=0.15 ** 3, Ks3=10 ** -20.2, Kelox=2**2, KpyrO2=0.01475, T=10, VolRel_m=0.75,
               VolRel_b=0.45, Sc=25, Sa=7, M3=0.96, ATP_temp=5, Pi_temp=5, H_temp=10 ** -7, BMR=350, ATP_Ruhe=350/18.2/60*0.233):

    ATP = ATP_temp
    H = H_temp
    Pi = Pi_temp
    t_span = np.array([0, 1])
    CHEP_System = np.array([y0s[0]])
    OxP = np.array([y0s[1]])
    Lactate_Muscle = np.array([y0s[2]])
    Lactate_Blood = np.array([y0s[3]])
    ATP_store = np.array([ATP_temp])
    pHm_store = np.array([np.log10(H_temp) * -1])
    AMP_store = np.array([0.0001])
    ADP_store = np.array([0.0001])

    # converting units
    vlamax = vlamax / 60
    VO2max = VO2max / 60
    Vrel = Muscle_Aktiv / (VolRel_b - Muscle_Aktiv)


    for n in range(0, len(ATP_Belastung)):

        # Phosphate
        if n == 0:
            PCr, Pi, ATP, ADP, AMP, M1, H, pHm = Adenosine_Phosphates(CHEP_System, ATP, y0s[1], VO2max, Pi, y0s[2], Sa=7, Sc=25, M3=0.96)
            kdiff = 0.065 * (y0s[3] ** -1.4)

        else:
            PCr, Pi, ATP, ADP, AMP, M1, H, pHm = Adenosine_Phosphates(x.y[0][1], ATP, x.y[1][1], VO2max, Pi, x.y[2][1], Sa=7, Sc=25, M3=0.96)
            kdiff = 0.065 * (x.y[3][1] ** -1.4)

        # RK4
        x = solve_ivp(Model, t_span, y0s, args=(
            ATP_Belastung[n], VO2max, vlamax, Body_weight, Musclemass, Muscle_Aktiv, Ks1, Ks2, Ks3, Kelox, KpyrO2, T,
            VolRel_m, VolRel_b, Sc, Sa, M3, ATP_Ruhe, Vrel, PCr, Pi, ATP, ADP, AMP, M1, H, kdiff), first_step=1, method="RK45", dense_output=True, vectorized=True)

        # storing data
        CHEP_System = np.append(CHEP_System, x.y[0][1])
        OxP = np.append(OxP, x.y[1][1])
        Lactate_Muscle = np.append(Lactate_Muscle, x.y[2][1])
        Lactate_Blood = np.append(Lactate_Blood, x.y[3][1])
        ATP_store = np.append(ATP_store, ATP)
        pHm_store = np.append(pHm_store, pHm)
        AMP_store = np.append(AMP_store, AMP)
        ADP_store = np.append(ADP_store, ADP)

        # Setting new start values
        y0s = np.array([x.y[0][1], x.y[1][1], x.y[2][1], x.y[3][1]])

    # reporting ending status
    if x.message == 'The solver successfully reached the end of the integration interval.':
        print('Calculation successfull')
    else:
        print('Beware. An Error occured during the calculation!')

    return CHEP_System, OxP, Lactate_Muscle, Lactate_Blood, ATP_store, pHm_store, AMP_store, ADP_store

def fitfile(VO2max, vLamax, Eco, Ks1, Ks2, VolRel, bw):
    file = askopenfilename()
    fitfile = fitparse.FitFile(str(file))

    datas = np.array([])
    # Iterate over all messages of type "record"
    # (other types include "device_info", "file_creator", "event", etc)
    for record in fitfile.get_messages("record"):
        # Records can contain multiple pieces of data (ex: timestamp, latitude, longitude, etc)
        for data in record:
            #    print(data.value)
            # Print the name and value of the data (and the units if it has any)
            if data.name == 'power':
                a = data.value
                if a == None:
                    datas = np.append(datas, datas[len(datas) - 1])
                else:
                    datas = np.append(datas, a)
    # Plotting Power data
    time = np.arange(0, len(datas), 1) / 60
    plt.plot(time, datas)
    plt.xlabel('Time [min]')
    plt.ylabel('Power [W]')
    plt.show()


    fat = ((100 * VolRel - 55.13) * - 2)
    active_muscle = (-0.42 * fat + 46.2) * 0.65 / 100
    VO2max_mm = VO2max * bw / (active_muscle * bw)
    VLamax_mm = (vLamax * (VolRel * bw) * 60) / (active_muscle * bw)
    datas = datas / (active_muscle * bw) * (0.048541/12.5 * (Eco * bw))
    CHEP_System, OxP, Lactate_Muscle, Lactate_Blood, ATP_store, pHm_store, AMP_store, ADP_store = call_Mader(datas,
                                                                                                             VO2max=VO2max_mm,
                                                                                                             vlamax=VLamax_mm,
                                                                                                             VolRel_b=VolRel,
                                                                                                             Body_weight=bw,
                                                                                                             Muscle_Aktiv=active_muscle,
                                                                                                             ATP_Ruhe=3.5 * bw / 60 / (
                                                                                                            active_muscle * bw) * 0.233)

    plt.plot(np.arange(0, len(CHEP_System), 1) / 60, CHEP_System-ATP_store)
    plt.axhline(3, color='grey', linestyle='--')
    plt.show()
    H = 10 ** - pHm_store
    Ks3 = 10 ** -20.2
    Ks2 = 0.15 ** 3

    vlamax_pc = (((1 / (1 + (H ** 3 / Ks3))) * (VLamax_mm / (1 + (Ks2 / (ADP_store * AMP_store)))))) / VLamax_mm * 100
    VO2_pc = OxP / (VO2max_mm/60) * 100

    time = np.arange(0,len(vlamax_pc), 1) / 60
    plt.plot(time, vlamax_pc, label='VLamax')
    plt.plot(time, VO2_pc, label='VO2max')
    plt.legend()
    plt.xlabel('Time [min]')
    plt.ylabel('%')
    plt.show()

    plt.plot(time, Lactate_Muscle, label='Lactate Muscle')
    plt.plot(time, Lactate_Blood, label='Lactate Blood')
    plt.xlabel('Time [min]')
    plt.ylabel('mmol/L')
    plt.legend()
    plt.show()


    overall_demand = datas + (3.5 * bw / 60 / (active_muscle * bw) * 0.233)

    vLamax_con = (((1 / (1 + (np.delete(H, 0)  ** 3 / Ks3))) * (VLamax_mm / (1 + (Ks2 / (np.delete(ADP_store,0)  * np.delete(AMP_store,0))))))) * 1.4 / 60
    OxP_con = np.delete(OxP, 0)  * 0.233
    Phosphates = overall_demand - vLamax_con  - OxP_con
    time = np.arange(0, len(Phosphates), 1) / 60
    plt.stackplot(time, vLamax_con/overall_demand*100, OxP_con/overall_demand*100, Phosphates/overall_demand*100, labels=['Gly.','OxP.','CHEP'])
    plt.legend()
    plt.ylim(0,100)
    plt.show()

    Overall_O2 =np.round(np.sum(OxP) * (active_muscle * bw) / 1000, 1) # in Litre
    Overall_mean_O2pc = np.round(np.mean(VO2_pc), 1)
    Overall_mean_vLa = np.round(np.mean(vlamax_pc), 1)

    CHO = np.sum((((1 / (1 + (H ** 3 / Ks3))) * ((VLamax_mm/60) / (1 + (Ks2 / (ADP_store * AMP_store)))))) / 1000 / 2 * 162.14 * (bw*active_muscle))
    #fat = np.sum(OxP * (bw*active_muscle)) / 1000 * 4.65 / 9.5



    root_temp = Tk()
    root_temp.title("Summary")
    sum = Label(root_temp, text=f'Mean % of VO2max: {Overall_mean_O2pc} %\n'
                                f'Mean % of VLamax: {Overall_mean_vLa} %\n'
                                f'Overall consumed oxygen: {Overall_O2} L \n '
                                f'Overall CHO combustion: {round(CHO,1)} g')
    sum.pack()

    root_temp.mainloop()

# GUI

root = Tk()
root.title("Metabolic Simulator v2")
root.geometry("200x300")

## creating entry boxes
entry_VO2max = Entry(root)
entry_VO2max.insert(0, "VO2max")
entry_VLamax = Entry(root)
entry_VLamax.insert(0, "VLamax")
entry_Eco = Entry(root)
entry_Eco.insert(0, "O2 per Unit of Loading [ml/min/kg]")
entry_Volrel = Entry(root)
entry_Volrel.insert(0, "Volrel")
entry_Bodyweight = Entry(root)
entry_Bodyweight.insert(0, "Body Weight [KG]")
entry_KS1 = Entry(root)
entry_KS1.insert(0, 'Ks1')
entry_KS2 = Entry(root)
entry_KS2.insert(0, 'Ks2')
entry_KEL = Entry(root)
entry_KEL.insert(0, 'Kel')

# Visualising Entry boxes
## 1. Column - Personal Data
label = Label(root, text="Create your Athlete")
label.pack()
entry_VO2max.pack()
entry_VLamax.pack()
entry_Eco.pack()
entry_Bodyweight.pack()
entry_Volrel.pack()

## 2. Column - Constants
label_optional = Label(root, text='Constants (optional)')
label_optional.pack()
entry_KS1.pack()
entry_KS2.pack()
entry_KEL.pack()

# Creating Buttons
CLass_Button = Button(root, text="CLass", width=15, height=1, command=lambda: plotClass(float(get_values()[0]), float(get_values()[1]), float(get_values()[2]), float(get_values()[5]), float(get_values()[6]), float(get_values()[3]), float(get_values()[7]), float(get_values()[4])))
LOP_Button = Button(root, text="  AT  ", width=15, height=1, command=lambda: plotLactate(float(get_values()[0]), float(get_values()[1]), float(get_values()[2]), float(get_values()[5]), float(get_values()[6]), float(get_values()[3]), float(get_values()[4])))
Macro_Button = Button(root, text= "Macronutrients", width=15, height=1, command=lambda: Macronutrients(float(get_values()[0]), float(get_values()[1]), float(get_values()[2]), float(get_values()[5]), float(get_values()[6]), float(get_values()[3]), float(get_values()[4])))
ImportFIT_Button = Button(root, text= "Import Fitfile...", width=15, height=1, command=lambda: fitfile(float(get_values()[0]), float(get_values()[1]), float(get_values()[2]), float(get_values()[5]), float(get_values()[6]), float(get_values()[3]), float(get_values()[4])))
#Diagnostik = Button(root, text= "Diagnostic", width=15, height=1, command=lambda: fitfile(float(get_values()[0]), float(get_values()[1]), float(get_values()[2]), float(get_values()[5]), float(get_values()[6]), float(get_values()[3]), float(get_values()[4])))

LOP_Button.pack()
CLass_Button.pack()
Macro_Button.pack()
ImportFIT_Button.pack()
#Diagnostik.pack()

root.mainloop()

