# SW and LW range for plot
l_sw_min = 0.3e-4 # cm
l_sw_max = 3.0e-4 # cm
l_lw_min = 3.0e-4 # cm
l_lw_max = 100.0e-4 # cm

def add_SW_LW_plot_areas(ax):
    ax.axvspan(l_sw_min, l_sw_max, facecolor='lightskyblue', alpha=0.2, label='Shortwave region')
    ax.axvspan(l_lw_min, l_lw_max , facecolor='lightcoral', alpha=0.2, label='Longwave region')

def wavelength_x(ax):
    ax.set_xlabel(f'Wavelength $\lambda$ / cm')