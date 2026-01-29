import matplotlib.pyplot as plt
import numpy as np

def printMEACztpFromCLAST():
    f = open('cloudy_clima/CLIMA/IO/clima_last.tab','r')
    data = f.readlines()
    f.close()

    alt,t,p = [],[],[]
    for line in data[1:]:
        line = line.split()
        alt.append(float(line[0]))
        p.append(np.log10(float(line[1])*101325))
        t.append(float(line[2]))
    alt = alt[::-1]
    p = p[::-1]
    t = t[::-1]

    zgap = alt[-1]-alt[-2]
    pgap = p[-1]-p[-2]

    z = alt[-1]
    while z < 100:
        z += zgap
        alt.append(z)
        p.append(p[-1]+pgap)
        t.append(t[-1])

    for i in range(len(p)):
        aa = np.round(alt[i],decimals=6)
        pp = np.round(p[i],decimals=6)
        tt = np.round(t[i],decimals=6)
        print(f"{aa}    {pp}      {tt}")

def compMEACandPhotochem():
    fig,ax = plt.subplots(ncols=1)

    ### T-P Profiles
    
    #   Coupled Model
    ztp_meac = open('hu-code-sr/scenario_library/guzman-marmolejo/fco2_1e-1/TP.dat','r')
    z,t,p = [],[],[]
    for line in ztp_meac:
        line = line.split()
        z.append(float(line[0]))
        t.append(float(line[2]))
        p.append(10**float(line[1]))
    ax.plot(t,p,c='black',label='CLIMA/MEAC')

    #   Photochem
    ztp_pc = open('../../pc/ex_evoatmosphere/atmosphere.txt','r')
    z,t,p = [],[],[]
    for line in ztp_pc.readlines()[1:]:
        line = line.split()
        t.append(float(line[3]))
        p.append(float(line[1])*1e5)
    ax.plot(t,p,c='red',label="Photochem")

    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("Pressure [Pa]")
    ax.set_ylim(1.5e5,1)
    ax.set_yscale('log')
    ax.legend()
    
    ### Kzz profiles

    # #   Coupled model
    # kzz_meac = open('hu-code-sr/scenario_library/guzman-marmolejo/fco2_1e-1/kzz.dat','r')
    # z,kzz = [],[]
    # for line in kzz_meac:
    #     line = line.split()
    #     z.append(float(line[0]))
    #     kzz.append(float(line[1]))
    # ax[1].plot(kzz,z,c='black',label="CLIMA/MEAC")

    # #   Photochem
    # z,kzz = [],[]
    # for line in ztp_pc.readlines()[1:]:
    #     z.append(float(line[0]))
    #     kzz.append(float(line[4]))
    # ax[1].plot(kzz,z,c='red',label="Photochem")

    ### SED's

    # #   MEAC
    # #   The solar spectrum employed by MEAC uses 1nm-wide bins up til 630, then goes to 2nm-wide bins
    # spec_meac = open('hu-code-sr/Data/solar00.txt','r')
    # wvln,flux = [],[]
    # for line in spec_meac:
    #     line = line.split()
    #     wvln.append(float(line[0]))
    #     flux.append(float(line[1]))
    # ax.plot(wvln,flux,c='black',label="CLIMA/MEAC")

    # #   Photochem
    # #   Spectral bins are 1nm-wide up to 120nm, then 0.05nm-wide thru 400nm and beyond
    # spec_pc = open('../../pc/ex_evoatmosphere/Sun_now.txt','r')
    # wvln,flux = [],[]
    # for line in spec_pc.readlines()[1:]:
    #     line = line.split()
    #     wvln.append(float(line[0]))
    #     flux.append(float(line[1])/1e3)
    # ax.plot(wvln,flux,c='red',label="Photochem",zorder=-100)

    # ax.set_xlabel("Wavelength [nm]")
    # ax.set_ylabel("Intensity [W/m^2/nm]")
    # # ax.set_yscale("log")
    # ax.set_xlim(100,400)
    # ax.legend()

    plt.tight_layout()
    plt.show()

def writePhotochemTPasMEACTP():
    #   Photochem TP profile
    ztp_pc = open('../../pc/ex_evoatmosphere/atmosphere.txt','r')
    z,t,p = [],[],[]
    for line in ztp_pc.readlines()[1:]:
        line = line.split()
        z.append(float(line[0]))
        t.append(float(line[3]))
        p.append(np.log10(float(line[1])*1e5))

    #   Fill missing altitudes up to 100km by drawing isotherm
    zz = z[-1]
    pp = p[-1]
    dz = z[-1]-z[-2]
    dp = p[-1]-p[-2]
    while zz<100:
        zz += dz
        pp += dp
        z.append(zz)
        t.append(t[-1])
        p.append(pp)

    #   MEAC format TP profile
    f = open('TP.dat','w')
    for i in range(len(t)):
        f.write(np.format_float_positional(z[i],precision=6,min_digits=6)+'\t')
        f.write(np.format_float_positional(p[i],precision=6,min_digits=6)+'\t')
        f.write(np.format_float_positional(t[i],precision=6,min_digits=6)+'\n')


    f.close()

def binPhotochemSpectrum():
    # Bin the solar spectrum from Photochem into bins with the same centers and
    #   widths as the bins in the MEAC solar spectrum. Only 100nm-400nm.

    fig,ax = plt.subplots(nrows=2)

    #   Photochem
    #   Spectral bins are 1nm-wide up to 120nm, then 0.05nm-wide thru 400nm and beyond
    spec_pc = open('../../pc/ex_evoatmosphere/Sun_now.txt','r')
    wvln,flux = [],[]
    for line in spec_pc.readlines()[1:]:
        line = line.split()
        wvln.append(float(line[0]))
        flux.append(float(line[1])/1e3)
    ax[0].plot(wvln,flux,c='grey',label="Photochem")

    new_wvln,new_flux = [],[]
    for i in range(100,401):
        bin_total = 0
        bin_count = 0
        for j in range(len(wvln)):
            if np.floor(wvln[j])==i:
                bin_total += flux[j]
                bin_count += 1
        new_wvln.append(i+0.5)
        new_flux.append(bin_total/bin_count)
    # ax[0].plot(new_wvln,new_flux,c='red',label="Photochem, binned")

    #   MEAC
    #   The solar spectrum employed by MEAC uses 1nm-wide bins up til 630, then goes to 2nm-wide bins
    spec_meac = open('hu-code-sr/Data/solar0.txt','r')
    wvln,flux = [],[]
    for line in spec_meac.readlines()[:301]:
        line = line.split()
        wvln.append(float(line[0]))
        flux.append(float(line[1]))
    ax[0].plot(wvln,flux,c='black',label="CLIMA/MEAC")

    #   Relative difference of Photochem from CLIMA/MEAC
    ax[1].hlines(1,xmin=100,xmax=400,colors='black')
    new_flux = np.asarray(new_flux)
    flux = np.asarray(flux)
    ax[1].plot(wvln,new_flux/flux,c='red',label="Photochem, binned")
    ax[1].set_xlim(100,160)
    
    ave = np.mean(new_flux/flux)
    print(ave)

    ax[0].set_ylabel("Flux [W/m^2/nm]")
    ax[0].set_yscale('log')
    ax[0].legend()
    ax[0].set_xlim(100,160)
    ax[1].set_xlabel("Wavelength [nm]")
    ax[1].set_ylabel("Relative flux")
    ax[1].legend()
    plt.tight_layout()
    plt.savefig('fig',dpi=300)
    plt.show()

# printMEACztpFromCLAST()
writePhotochemTPasMEACTP()
# compMEACandPhotochem()
# binPhotochemSpectrum()

