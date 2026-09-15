import matplotlib.pyplot as plt
import numpy as np

meac_start      = 'hu-code-sr/scenario_library/co2_1e-6_ch4_1e8/ConcentrationSTD_base.dat'
meac_final      = 'hu-code-sr/scenario_library/co2_1e-6_ch4_1e8/ConcentrationSTD.dat'
clima_ztp_start = 'outputs/co2_1e-6_ch4_1e8/clima-out/ztp_0.dat'
clima_ztp_final = 'outputs/co2_1e-6_ch4_1e8/clima-out/ztp_14.dat'
clima_sp_final_root = 'cloudy_clima/CLIMA/IO/Profiles_out/'

fig,ax = plt.subplots(nrows=2,ncols=2)

# Axis [0,0]: Initial TP profiles
init_clima = np.genfromtxt(clima_ztp_start,skip_header=1)
init_meac   = np.genfromtxt(meac_start,skip_header=2)
ax[0,0].plot(init_clima[:,1],init_clima[:,2],c='black',label='Cloudy-CLIMA')
ax[0,0].plot(init_meac[:,3],init_meac[:,4],c='black',linestyle='dashed',label='MEAC')
ax[0,0].set_yscale('log')
ax[0,0].invert_yaxis()

# Axis [1,0]: Initial mixing ratios (none for CLIMA
# CLIMA
pres = init_clima[:,2]
ax[1,0].plot(init_clima[:,3],pres,c='gold',label='C2H6')
ax[1,0].plot(init_clima[:,4],pres,c='orange',label='CH4')
ax[1,0].plot(init_clima[:,5],pres,c='red',label='CO2')
ax[1,0].plot(init_clima[:,6],pres,c='purple',label='H2')
ax[1,0].plot(init_clima[:,7],pres,c='cyan',label='H2O')
ax[1,0].plot(init_clima[:,8],pres,c='green',label='N2')
ax[1,0].plot(init_clima[:,9],pres,c='blue',label='O2')
ax[1,0].plot(init_clima[:,10],pres,c='navy',label='O3')

# MEAC
nsum = np.sum(init_meac[:,4:],axis=1)
ax[1,0].plot(init_meac[:,35]/nsum,init_meac[:,4],c='gold',linestyle='dashed')
ax[1,0].plot(init_meac[:,25]/nsum,init_meac[:,4],c='orange',linestyle='dashed')
ax[1,0].plot(init_meac[:,56]/nsum,init_meac[:,4],c='red',linestyle='dashed')
ax[1,0].plot(init_meac[:,57]/nsum,init_meac[:,4],c='purple',linestyle='dashed')
ax[1,0].plot(init_meac[:,11]/nsum,init_meac[:,4],c='cyan',linestyle='dashed')
ax[1,0].plot(init_meac[:,59]/nsum,init_meac[:,4],c='green',linestyle='dashed')
ax[1,0].plot(init_meac[:,58]/nsum,init_meac[:,4],c='blue',linestyle='dashed')
ax[1,0].plot(init_meac[:,6]/nsum,init_meac[:,4],c='navy',linestyle='dashed')

# init_meac   = np.genfromtxt('outputs/co2_1e-6_ch4_1e8/meac-out/mr_0.dat',skip_header=1)
# ax[1,0].plot(init_meac[:,4],init_meac[:,2],c='gold',linestyle='dashed')
# ax[1,0].plot(init_meac[:,6],init_meac[:,2],c='orange',linestyle='dashed')
# ax[1,0].plot(init_meac[:,8],init_meac[:,2],c='red',linestyle='dashed')
# ax[1,0].plot(init_meac[:,10],init_meac[:,2],c='purple',linestyle='dashed')
# ax[1,0].plot(init_meac[:,12],init_meac[:,2],c='cyan',linestyle='dashed')
# ax[1,0].plot(init_meac[:,14],init_meac[:,2],c='green',linestyle='dashed')
# ax[1,0].plot(init_meac[:,16],init_meac[:,2],c='blue',linestyle='dashed')
# ax[1,0].plot(init_meac[:,18],init_meac[:,2],c='navy',linestyle='dashed')

ax[1,0].set_xscale('log')
ax[1,0].set_yscale('log')
ax[1,0].set_xlim(left=1e-18)
ax[1,0].invert_yaxis()

# Axis [0,1]: Final TP profiles
init_clima = np.genfromtxt(clima_ztp_final,skip_header=1)
init_meac   = np.genfromtxt(meac_final,skip_header=2)
ax[0,1].plot(init_clima[:,1],init_clima[:,2],c='black',label='Cloudy-CLIMA')
ax[0,1].plot(init_meac[:,3],init_meac[:,4],c='black',linestyle='dashed',label='MEAC')
ax[0,1].set_yscale('log')
ax[0,1].invert_yaxis()

# Axis [1,1]: Final mixing ratios
init_clima = np.genfromtxt(clima_ztp_final,skip_header=1)
init_meac   = np.genfromtxt('outputs/co2_1e-6_ch4_1e8/meac-out/mr_14.dat',skip_header=1)

# CLIMA
clima_sp = np.genfromtxt(clima_sp_final_root+'C2H6.dat')
ax[1,1].plot(clima_sp,init_clima[:,2],c='gold',label='C2H6')
clima_sp = np.genfromtxt(clima_sp_final_root+'CH4.dat')
ax[1,1].plot(clima_sp,init_clima[:,2],c='orange',label='CH4')
clima_sp = np.genfromtxt(clima_sp_final_root+'CO2.dat')
ax[1,1].plot(clima_sp,init_clima[:,2],c='red',label='CO2')
clima_sp = np.genfromtxt(clima_sp_final_root+'H2.dat')
ax[1,1].plot(clima_sp,init_clima[:,2],c='purple',label='H2')
clima_sp = np.genfromtxt(clima_sp_final_root+'H2O.dat')
ax[1,1].plot(clima_sp,init_clima[:,2],c='cyan',label='H2O')
clima_sp = np.genfromtxt(clima_sp_final_root+'N2.dat')
ax[1,1].plot(clima_sp,init_clima[:,2],c='green',label='N2')
clima_sp = np.genfromtxt(clima_sp_final_root+'O2.dat')
ax[1,1].plot(clima_sp,init_clima[:,2],c='blue',label='O2')
clima_sp = np.genfromtxt(clima_sp_final_root+'O3.dat')
ax[1,1].plot(clima_sp,init_clima[:,2],c='navy',label='O3')

# MEAC
ax[1,1].plot(init_meac[:,4],init_meac[:,2],c='gold',linestyle='dashed')
ax[1,1].plot(init_meac[:,6],init_meac[:,2],c='orange',linestyle='dashed')
ax[1,1].plot(init_meac[:,8],init_meac[:,2],c='red',linestyle='dashed')
ax[1,1].plot(init_meac[:,10],init_meac[:,2],c='purple',linestyle='dashed')
ax[1,1].plot(init_meac[:,12],init_meac[:,2],c='cyan',linestyle='dashed')
ax[1,1].plot(init_meac[:,14],init_meac[:,2],c='green',linestyle='dashed')
ax[1,1].plot(init_meac[:,16],init_meac[:,2],c='blue',linestyle='dashed')
ax[1,1].plot(init_meac[:,18],init_meac[:,2],c='navy',linestyle='dashed')

ax[1,1].set_xscale('log')
ax[1,1].set_xlim(left=1e-18)
ax[1,1].set_yscale('log')
ax[1,1].invert_yaxis()

ax[0,0].set_title("Initial",size='x-large')
ax[0,1].set_title("Final",size='x-large')

for a in ax[0].flatten():
    a.set_xlabel("Temperature [K]",size='large')

for a in ax[1].flatten():
    a.set_xlabel("Mixing ratio",size='large')

for a in ax.flatten():
    a.legend()
    a.set_ylim(1e5,1e-2)
    a.set_ylabel("Pressure [Pa]",size='large')

fig.set_figwidth(12)
fig.set_figheight(8)
plt.tight_layout()
plt.savefig('convergence_tp_mr.png',dpi=200)
plt.show()