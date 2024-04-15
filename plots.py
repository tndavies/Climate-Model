import matplotlib.pyplot as plt
from climate import *
import scienceplots
import numpy as np
from data import Wildfire_Damages

def Get_ClimateRecordLength():
    return list(Historic_Temperatures)[-1] - list(Historic_Temperatures)[0]

def Sample_Years(arr: list):
    return np.array(arr)[::365]

def Save_Figure(fig, filename: str):
    path = ".\\" + filename + ".pdf"
    fig.savefig(path, dpi=96, bbox_inches="tight")

plt.style.use('science')
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.size'] = 12

# ========================================

def Fig_GlobalTemperatureCalibration():    
    Sim_Duration = Get_ClimateRecordLength()

    Uniform_Sim = Simulate_Climate(Sim_Specification(Duration=Sim_Duration, InitialTempDist=None))
    Stablised_Sim = Simulate_Climate(Sim_Specification(Duration=Sim_Duration, InitialTempDist=Equilibrium_Config))
    Obs_Times, Obs_GATs = Serialise(Historic_Temperatures)
    
    fig, (ax_uniform, ax_stable) = plt.subplots(nrows=2,ncols=1,sharex=True)

    fig.supxlabel("Year")
    fig.supylabel("Temperature (K)")

    Uniform_GATs = Average(Uniform_Sim, lambda x: x, (-90, 90))
    ax_uniform.plot(Sample_Years(Uniform_Sim.times), Sample_Years(Uniform_GATs), "-", color="red") 
    ax_uniform.plot(Obs_Times, Obs_GATs, "k--")

    Stable_GATs = Average(Stablised_Sim, lambda x: x, (-90, 90))
    ax_stable.plot(Sample_Years(Stablised_Sim.times), Sample_Years(Stable_GATs), "-", color="seagreen") 
    ax_stable.plot(Obs_Times, Obs_GATs, "k--")

    plt.show()

    # Print out the comparison metrics for the uniform and stable calibration cases.
    print("\n --- Comparison Metrics ---")
    print("Initial Uniform Temperature Distribution: {X}".format(X=Calc_ComparisonMetric(Uniform_Sim)))
    print("Initial Stable Temperature Distribution: {X}".format(X=Calc_ComparisonMetric(Stablised_Sim)))
    
def Fig_RCPs():
    fig, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False)

    axis.set_ylabel("Concentration (ppm)")
    axis.set_xlabel("Year")

    # Plots the CO2 data
    Observation_Times, Observed_Concentrations = Serialise(Historic_Co2)
    axis.plot(Observation_Times, Observed_Concentrations, "--", color="black")

    # Plots RCP 8.5
    Times, RCP85_Concentrations = Serialise(RCP85_Data)
    axis.plot(Times, RCP85_Concentrations, color="firebrick")

    # Plots RCP 6.0
    Times, Concentrations = Serialise(RCP6_Data)
    axis.plot(Times, Concentrations, color="darkorange")

    # Plots RCP 4.5
    Times, Concentrations = Serialise(RCP45_Data)
    axis.plot(Times, Concentrations, color="royalblue")

    # Plots RCP 2.6
    Times, RCP26_Concentrations = Serialise(RCP26_Data)
    axis.plot(Times, RCP26_Concentrations, color="seagreen")

    axis.fill_between(Times, RCP26_Concentrations, RCP85_Concentrations, color="slategrey", alpha=0.25)

    plt.show()
    
def Fig_ClimateRecords():
    fig, (co2, gat) = plt.subplots(nrows=2,ncols=1,sharex=True)

    fig.supxlabel("Year")

    Times, Gats = Serialise(Historic_Temperatures)
    gat.plot(Times, Gats, "k-")
    gat.set_ylabel("Temperature (K)")

    Times, Concentrations = Serialise(Historic_Co2)
    co2.plot(Times, Concentrations, "k-")
    co2.set_ylabel("Concentration (ppm)")

    plt.show()
    
def Fig_AlbedoModel():
    fig, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False)

    axis.set_xlabel("Temperature (K)")
    axis.set_ylabel("Albedo")

    temps = np.linspace(180,300,1000)
    albedos = [calc_Albedo(T) for T in temps]
    
    axis.plot(temps, albedos, "-", color="black")
    axis.axvline(273.0, linestyle="--", color="blue", alpha=0.45)

    plt.show()
    
def Fig_Irradiance():
    fig, (ax_Q, ax_irr, ax_dec) = plt.subplots(nrows=3,ncols=1,sharex=True)
    fig.supxlabel("Day")

    times = np.linspace(0, 86400*365, 1000)
    Plot_Times = np.divide(times, 86400)
    
    qstar = [calc_RadiationIntensity(t) for t in times]
    ax_Q.plot(Plot_Times, qstar, "-", color="black")
    ax_Q.axhline(1361.0, linestyle="--", color="black", alpha=0.4)
    ax_Q.set_ylabel(r"Q$_\star$ (Wm$^{-2}$)")

    decl = [calc_Declination(t) for t in times]
    ax_dec.plot(Plot_Times, np.degrees(decl), "-", color="black")
    ax_dec.set_ylabel(r"$\delta$")

    ax_irr.set_ylabel(r"S (Wm$^{-2}$)")
    lats = [-np.pi/2, 0, np.pi/2]
    cols = ["darkblue", "lightcoral", "cornflowerblue"]
    styles = ["--", "-", "--"]
    for idx,L in enumerate(lats):
        diurnal_fluxes = [calc_AverageRadiation(L, t) for t in times]
        ax_irr.plot(Plot_Times, diurnal_fluxes, linestyle=styles[idx], color=cols[idx])

    plt.show()

def Fig_AbsorptionModel():
    fig, (ax_vapour,ax_potency) = plt.subplots(nrows=1,ncols=2,sharex=False)

    # fixed Co2-level, varying temp
    temps = np.linspace(220, 300, 1000)	
    absoprtions = [calc_AtmosphericAbsorption(T, Co2_Norm, 0.0) for T in temps]
    ax_vapour.plot(temps, absoprtions, "-", color="black")
    ax_vapour.set_ylabel(r"$\Gamma$")
    ax_vapour.set_xlabel("Temperature (K)")
 
    # Fix temp and vary CO2 level, using optimal potency
    pCO2 = np.linspace(Co2_Norm, 1000, 1000)
    absoprtions = [calc_AtmosphericAbsorption(Temp_Norm, CO2_Level, Best_Beta) for CO2_Level in pCO2]
    ax_potency.plot(pCO2, absoprtions, "-", color="red", linewidth=1.5, label=str(Best_Beta))
  
    # Fix temp and vary CO2 level, using differing potencies
    potencies = np.linspace(0,0.4,4)
    for idx, B in enumerate(potencies):
        absoprtions = [calc_AtmosphericAbsorption(Temp_Norm, CO2_Level, B) for CO2_Level in pCO2]
        ax_potency.plot(pCO2, absoprtions, "--", linewidth=1.8, color="black", alpha=1/(1+idx**2.2))
 
    ax_potency.set_ylabel(r"$\Gamma$")
    ax_potency.set_xlabel("Concentration (ppm)")

    plt.show() 

def Fig_HeatCapacityTest():
    fig, (ax_hcaps, ax_ice) = plt.subplots(nrows=2,ncols=1,sharex=False)

    fig.supxlabel("Temperature (K)")

    temps = np.linspace(180, 300, 1000)
    Ocean_Lat, Land_Lat = np.radians(-55.0), np.radians(-90.0)
    ocean_hc = [calc_AverageHeatCapacity(Ocean_Lat, T) for T in temps]
    land_hc = [calc_AverageHeatCapacity(Land_Lat, T) for T in temps]
    # hc_ices = [C_Transitioning_Ice if(T >= 263 and T <= 273) else C_Land for T in temps]

    ax_hcaps.plot(temps, np.divide(ocean_hc,1e6), "-", color="steelblue", linewidth=2.0)
    ax_hcaps.plot(temps, np.divide(land_hc, 1e6), "-", color="peru", linewidth=2.0)
    # ax_hcaps.plot(temps, np.divide(hc_ices, 1e6), "--", color="black", alpha=0.5) 
    ax_hcaps.set_ylabel(r"$\tilde{C}$ (MJ K$^{-1}$)", fontsize=17)

    ice_fracs = [Calc_IceFraction(T) for T in temps]
    ax_ice.plot(temps, ice_fracs, "-", color="black")
    ax_ice.set_ylabel(r"f$_i$", fontsize=18)

    plt.show()
    
def Fig_Geography():
    fig, (ax) = plt.subplots(nrows=1,ncols=1,sharex=False)

    ax.set_xlabel(r"$\lambda$", fontsize=16)
    ax.set_ylabel(r"f$_L$", fontsize=16)
    
    lats = np.linspace(-90, 90, 1000)
    fl = [1.0 - Get_OceanFraction(lat, to_degrees=False) for lat in lats]
    ax.plot(lats, fl, "-", color="black")

    plt.show()

def Fig_AntarcticaCorrection():
    fig, (ax_temp, ax_alb) = plt.subplots(nrows=2,ncols=1,sharex=True)

    ax_temp.set_ylabel("Temperature (K)")
    ax_alb.set_ylabel("Albedo")
    fig.supxlabel("Year")

    Sim_Duration = Get_ClimateRecordLength()
    SamplesPerYear = 2
    N = int(365 / SamplesPerYear)

    Corrected_Sim = Simulate_Climate(Sim_Specification(Sim_Duration, Altitude_Correction=True, InitialTempDist=Equilibrium_Config))
    temps = Average(Corrected_Sim, lambda x: x, Antarctic_Bounds)
    albedos = Average(Corrected_Sim, calc_Albedo, Antarctic_Bounds)
    Corrected_AvgAlbedo = np.mean(albedos)
    ax_temp.plot(np.array(Corrected_Sim.times)[::N], np.array(temps)[::N], "-", color="black")
    ax_alb.plot(np.array(Corrected_Sim.times)[::N], np.array(albedos)[::N], ":", color="black", alpha=0.48)

    Uncorrected_Sim = Simulate_Climate(Sim_Specification(Sim_Duration, Altitude_Correction=False, InitialTempDist=Equilibrium_Config))
    temps = Average(Uncorrected_Sim, lambda x: x, Antarctic_Bounds)
    albedos = Average(Uncorrected_Sim, calc_Albedo, Antarctic_Bounds)
    Uncorrected_AvgAlbedo = np.mean(albedos)
    ax_temp.plot(np.array(Uncorrected_Sim.times)[::N], np.array(temps)[::N], "--", color="tomato")
    ax_alb.plot(np.array(Uncorrected_Sim.times)[::N], np.array(albedos)[::N], ":", color="tomato", alpha=0.48)

    ax_temp.axhline(y=273, linestyle="-", color='b', linewidth=1.2, alpha=0.4, label="Freezing point")
    ax_alb.axhline(y=Corrected_AvgAlbedo, linestyle="-", color='k', linewidth=2.5)
    ax_alb.axhline(y=Uncorrected_AvgAlbedo, linestyle="-", color='r', linewidth=2.5)

    plt.show()
    
def Fig_DipDependance():
    fig, (ax) = plt.subplots(nrows=1,ncols=1,sharex=False)

    # We simulate a few decades before the record starts, as to sinmulate the
    # climate in a constant atmospheric CO2 concentration, with the uniform
    # temperature IC, varying the CO2 to show how the dip responds.

    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Temperature (K)")

    Sim_Duration = 35
    for idx,Const_CO2 in enumerate(np.linspace(0.5*Co2_Norm, 2*Co2_Norm, 5)): 
        Sim = Simulate_Climate(Sim_Specification(Sim_Duration, Initial_Year=0, Prerecord_Co2Level=Const_CO2))
        Gats = Average(Sim, lambda x: x, (-90, 90))
        ax.plot(Sample_Years(Sim.times), Sample_Years(Gats), "k-", alpha=1/(1+idx**1.2))
    
    plt.show()

def Fig_Forecasts():
    Target_Year = 2100
    Duration = 1 + (Target_Year - list(Historic_Temperatures)[0])
    Sim_RCP85 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP85, InitialTempDist=Equilibrium_Config))
    Sim_RCP45 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP45, InitialTempDist=Equilibrium_Config))
    Sim_RCP6 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP6, InitialTempDist=Equilibrium_Config))
    Sim_RCP26 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP26, InitialTempDist=Equilibrium_Config))
    Sim_Context = Sim_RCP85
    
	# GAT Forecast plot
  	# ----------------------------------------------------------------
    fig, (ax) = plt.subplots(nrows=1,ncols=1,sharex=False)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Temperature (K)")

    LastYearOnRecord = (list(Historic_Temperatures)[-1] - Sim_Context.spec.Initial_Year)*31536000
    N = round(LastYearOnRecord / Sim_Context.spec.Time_Step)

    Times = np.array(Sim_RCP85.times)[N::365]
    GATs_RCP85 = np.array(Average(Sim_RCP85, lambda x: x, (-90,90)))[N::365]
    ax.plot(Times, GATs_RCP85, "-", color="firebrick")

    Times = np.array(Sim_RCP6.times)[N::365]
    GATs = np.array(Average(Sim_RCP6, lambda x: x, (-90,90)))[N::365]
    ax.plot(Times, GATs, "-", color="darkorange")

    Times = np.array(Sim_RCP45.times)[N::365]
    GATs = np.array(Average(Sim_RCP45, lambda x: x, (-90,90)))[N::365]
    ax.plot(Times, GATs, "-", color="royalblue")
    
    Times = np.array(Sim_RCP26.times)[N::365]
    GATs_RCP26 = np.array(Average(Sim_RCP26, lambda x: x, (-90,90)))[N::365]
    ax.plot(Times, GATs_RCP26, "-", color="seagreen")

    ax.fill_between(Times, GATs_RCP26, GATs_RCP85, color="slategrey", alpha=0.25)

    plt.show()
    
    # Temperature Distributions in the year 2100
  	# ----------------------------------------------------------------
    # fig, (ax) = plt.subplots(nrows=1,ncols=1,sharex=False)
    # ax.set_xlabel(r"$\lambda$", fontsize=16)
    # ax.set_ylabel("Temperature (K)")

    # def get_avgTDist(sim: Sim_Result):
    #     avg_td = []
    #     tdists_2100 = np.array(sim.tps)[-365:]
    #     for idx in range(len(Sim_Context.lats)):
    #         td_zonal_temps = [td[idx] for td in tdists_2100]
    #         avg_td.append(np.mean(td_zonal_temps))
    #     return avg_td
        
    # ax.plot(np.degrees(Sim_Context.lats), get_avgTDist(Sim_RCP85), "-", color="firebrick", linewidth=1.8) 
    # ax.plot(np.degrees(Sim_Context.lats), get_avgTDist(Sim_RCP6), "-", color="darkorange", linewidth=1.8) 
    # ax.plot(np.degrees(Sim_Context.lats), get_avgTDist(Sim_RCP45), "-", color="royalblue", linewidth=1.8) 
    # ax.plot(np.degrees(Sim_Context.lats), get_avgTDist(Sim_RCP26), "-", color="seagreen", linewidth=1.8) 

    # avg_td = []
    # tdists_LastYrOnRecord = np.array(Sim_Context.tps)[N-365:N]
    # for idx in range(len(Sim_Context.lats)):
    #     td_zonal_temps = [td[idx] for td in tdists_LastYrOnRecord]
    #     avg_td.append(np.mean(td_zonal_temps))
    # ax.plot(np.degrees(Sim_Context.lats), avg_td, "--", color="black", linewidth=2.2, alpha=0.5) 

    # plt.show()
     
def Fig_AntarcticaForecast():
    fig, (ax_temp, ax_alb) = plt.subplots(nrows=2,ncols=1,sharex=False)
    
    ax_temp.set_ylabel(r"$\sigma$ (\%)", fontsize=14)
    ax_alb.set_ylabel("Albedo")
    fig.supxlabel("Year")

    Target_Year = 2100
    Duration = Target_Year - list(Historic_Temperatures)[0]
    Sim_RCP85 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP85, InitialTempDist=Equilibrium_Config))
    Sim_RCP6 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP6, InitialTempDist=Equilibrium_Config))
    Sim_RCP45 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP45, InitialTempDist=Equilibrium_Config))
    Sim_RCP26 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP26, InitialTempDist=Equilibrium_Config))

    def plot(Sim: Sim_Result, col: str):	
        Antarctic_Temps = Average(Sim, lambda x: x, Antarctic_Bounds)
        Antarctic_Albedo = Average(Sim, calc_Albedo, Antarctic_Bounds)

        Plot_Offset = Get_ClimateRecordLength()
        Annual_Temps = np.array_split(Antarctic_Temps, Duration)
        Annual_Albedos = np.array_split(Antarctic_Albedo, Duration)

        # Finds the 365-day albedo average, for each year in the Annual_Albedos 
        Times, AvgAlbedos = [], []
        Year_Number = Sim.spec.Initial_Year
        for year in Annual_Albedos:
            AvgAlbedos.append(np.mean(year))
            Times.append(Year_Number)
            Year_Number += 1
        ax_alb.plot(np.array(Times)[Plot_Offset:], np.array(AvgAlbedos)[Plot_Offset:], "-", linewidth=1.5, color=col)

        # Counts how many days are >0C for each set of 365 dats in Annual_Temps
        Times, Frequencies = [], []
        Year_Number = Sim.spec.Initial_Year
        for year in Annual_Temps:
            DaysAboveFreezing = 0
            for T in year: 
                if T > 273.15: DaysAboveFreezing += 1
            Frequencies.append(100*(DaysAboveFreezing / 365))
            Times.append(Year_Number)
            Year_Number += 1

        Sample_Times = np.array(Times)[Plot_Offset:]
        Sample_Freqs = np.array(Frequencies)[Plot_Offset:]
        SmoothFit = Polynomial.fit(Sample_Times, Sample_Freqs, 2)
        Plot_Times = np.linspace(Sample_Times[0], Sample_Times[-1], 1000)
        ax_temp.plot(Plot_Times, SmoothFit(Plot_Times), linewidth=1.5, color=col)

    plot(Sim_RCP85, "firebrick")
    plot(Sim_RCP6, "darkorange")
    plot(Sim_RCP45, "royalblue")
    plot(Sim_RCP26, "seagreen")
     
    plt.show()
 
Fig_AntarcticaForecast()