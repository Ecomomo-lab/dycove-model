import time
from dycove import Reporter

r = Reporter()  # get singleton reporter instance

def print_model_time_info(simstate):
    r.report(f"Hydrodynamic model duration: {int(simstate.hydro_sim_days)} days")
    r.report(f"Eco-morphodynamic model duration: {simstate.veg_sim_years} years ")

def print_runtime_updates(simstate, i):
    elapsed_units = 1., "seconds", 0      
    simstate.times_elapsed.append(time.time() - simstate.time_0)
    if simstate.times_elapsed[-1] > 60:
        elapsed_units = 60., "minutes", 1
        if simstate.times_elapsed[-1] > 3600:
            elapsed_units = 3600., "hours", 2
    t_str = round(simstate.times_elapsed[-1]/elapsed_units[0], elapsed_units[2])

    r.report(f"Current hydrodynamic model time = {round(simstate.hydrotime_seconds/86400., 1)} days")
    r.report(f"Current eco-morphodynamic model time = {round(simstate.hydrotime_seconds/86400./simstate.days_per_year*simstate.vegfac, 1)} years")
    r.report(f"Current eco-morphodynamic model date = {simstate.vegtime_date}")
    r.report(f"Total time elapsed: {t_str} {elapsed_units[1]}")        
    
    proj_units_avg = 1., "seconds", 0
    if i == 0:
        time_per_ets_avg = simstate.times_elapsed[-1]        # (time so far)
    else:
        time_per_ets_avg = simstate.times_elapsed[-1] / i    # (time so far) / (n ETS so far)
    proj_time_avg = time_per_ets_avg * simstate.n_veg_steps  # (avg time per ETS) * (n ETS total)
    if proj_time_avg > 60:
        proj_units_avg = 60., "minutes", 1
        if proj_time_avg > 3600:
            proj_units_avg = 3600., "hours", 2

    t_str = round(proj_time_avg/proj_units_avg[0], proj_units_avg[2])
    r.report(f"Projected total run time based on average yield step: {t_str} {proj_units_avg[1]}"
                        "\n\n##### --------------------------------------------- #####\n")

