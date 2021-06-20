import random
import matplotlib.pyplot as plt
import numpy as np
from math import inf as infinity.  # CAREFUL! WHAT IS THIS USED FOR?

# 1. VARIABLE DEFINITION
N = "test"
s_part = 2
e_part = 0
catalyst = "Without"
if e_part == 2:
    catalyst = "With"
n_steps = 10000000         # NUMBER STEPS
p_part = 0                 # PRODUCT PARTICLES
n_part = s_part + e_part   # TOTAL NUMBER PARTICLES
dt = 0.006		  # TIME STEP -- TOO LARGE?
lambda_constant = 0.65     # FOR MODIFIED VELOCITY-VERLET
T_init = 0.7               # TEMPERATURE
n_part_inv = 1.0 / n_part
dt_over_2 = dt/2.0

box_width = 24.0                   # DIMENSIONS BOX
box_width_over_2 = box_width / 2.0

k_a = 2.0                  # STIFFNESS BREAKABLE
l_a = 0.0		  # REST_LENGTH BREAKABLE
z_a = 3.0                  # INTERACTION CUTOFF

k_r = 1.0                  # STIFFNESS UNBREAKABLE
l_r = 6.0                  # REST_LENGTH ONCE BROKEN


k_ar = k_a + k_r                     # SPRINGS ARE IN SERIES
l_ar = (k_r * l_r) / (k_a + l_r)     # SPRINGS ARE IN SERIES

k_e = 9999                 # STIFFNESS CATALYST -- NOT RELEVANT
l_e = z_a                  # REST_LENGTH CATALYST -- INTERACTION RANGE


k_i = 13                   # FLEXIBILITY CATALYST-SUBSTRATE INTERACTION
l_i = 0                    # REST LENGTH ""
z_i = (z_a - l_ar)/2.0     # INTERACTION RANGE

mass = 1.0                  # MASS - NOT RELEVANT
radius = 0.1               # PARTICLE RADIUS - OVITO PLOTTING

gamma = 10.0
gamma_inv = 1.0/gamma
sigma_brown = np.sqrt(2*T_init*gamma)

print "SIMULATION WITH T    :: ", T_init
print "FRICTION COEFFICIENT :: ", gamma
print "SIGMA NOISE          :: ", sigma_brown

# variables defined to save time
dt_2_over_2_mass = dt ** 2 / (2.0 * mass)
sigma_brown_divided_by_sqrt_dt = sigma_brown/np.sqrt(dt)
mass_n_part_inv = 1.0 / (mass * n_part)
inv_dt = 1.0/dt
s_part_list = [s_part]
k_ar_div_2 = k_ar * 0.5
k_i_div_2 = k_i * 0.5
k_r_div_2 = k_r * 0.5
e_zi = 0.5 * k_i * (z_i - l_i)**2
e_za = 0.5 * k_a * (z_a - l_a)**2
##print(e_za)

# 2. CREATE DATAFILE
nom_fichier = "datafile"+str(catalyst)+"Catalyst_dt_"+str(dt)+"_sigmaBrown_"+str(sigma_brown)+"_gamma_"+str(gamma)+"_Tinit_"+str(T_init)+"_"+N

with open ("Logfile_"+str(nom_fichier)+".txt","w") as logfile:
    logfile.write("##LANGEVIN BROWNIAN MOTION 1D SPRING REACTION WITHOUT CATALYST##"+"\n")
    logfile.write("Spart_"+str(s_part)+"\n")
    logfile.write("Epart_"+str(e_part)+"\n")
    logfile.write("Npart_"+str(n_part)+"\n")
    logfile.write("Nsteps_"+str(n_steps)+"\n")
    logfile.write("boxwidht_"+str(box_width)+"\n")
    logfile.write("Tinit_"+str(T_init)+"\n")
    logfile.write("dt_"+str(dt)+"\n")
    logfile.write("lambdaconstant_"+str(lambda_constant)+"\n")
    logfile.write("sigmaBrown_"+str(sigma_brown)+"\n")
    logfile.write("gamma_"+str(gamma)+"\n")
    logfile.write("mass_"+str(mass)+"\n")
    logfile.write("##Springs parameters##"+"\n")
    logfile.write("##Substrate sissile and non sissile bond##"+"\n")
    logfile.write("kar_"+str(k_ar)+"\n")
    logfile.write("lar_"+str(l_ar)+"\n")
    logfile.write("za_"+str(z_a)+"\n")
    logfile.write("##Substrate non-sissile bond##"+"\n")
    logfile.write("kr_"+str(k_r)+"\n")
    logfile.write("lr_"+str(l_r)+"\n")
    logfile.write("zr_infinite"+"\n")
    logfile.write("catalyst_0"+"\n")
    logfile.write("##No_catalyst_bonds##")


#ratio_file = open(str(nom_fichier)+"_ratio.txt","w")
#ratio_file.flush()
#ratio_file.write("Simulation_number;ForwardReactionTime"+"\n")

# ---------------------------
# INITIALIZATION OF ARRAYS
# - Two-dimensional motion --- RELEVANT ONLY FOR OVITO
y = np.zeros(n_part)
z = np.zeros(n_part)
vel_y = np.zeros(n_part)
vel_z = np.zeros(n_part)
# ---------------------------

for A in range(1):

    well_counts = 0
    First_well = True
    r_t_file = open(str(nom_fichier)+"_R(t).txt","w")
    r_t_file.flush()
    r_t_file.write("r;time;ep;temp"+"\n")

    time = 0
    time_list = []
    time_list.append(time)
    T = T_init
    Temp_average = T
    Temp_average_list=[T]
    Temp_list = [T]

    forward_reaction_time = None
    forward_reaction = False
    backward_reaction_time = None
    backward_reaction = False
    over_l_r = False


    # ---------------------------
    # INITIALIZATION OF ARRAYS
    # ---------------------------
    prev_x = np.zeros(n_part)
    x = np.zeros(n_part)
    next_x = np.zeros(n_part)
    vel_x = np.zeros(n_part)
    f_x = np.zeros(n_part)
    prev_f_x = np.zeros(n_part)

    # ---------------------------

    # ---------------------------
    # INITIALIZATION OF TYPES
    # ---------------------------
    type= s_part*['S'] + e_part*['E']
    type[0] = 'S_0'
    type[1] = 'S_1'
    if e_part == 2:
        type[2] = 'C_0'
        type[3] = 'C_1'

    # ---------------------------


    # ---------------------------
    # INITIALIZATION OF POSITIONS
    # ---------------------------
    sumv = 0.0
    sumv2 = 0.0
    e_tot = 0.0
    en = 0
    r_list=[]

    # MAKE TEST ---> SINGLE PARTICLE IN A WELL
    x[0] = np.sqrt(5/3)
    x[1] = 0
    r = x[0] - x[1]
    #r = x[0]

    if e_part == 2:
        x[2] = 0.0
        x[3] = x[2] + l_e
        type[2], type[3] = 'C_0','C_1'


    #type[0],type[1] = 'S','S'

    # ---------------------------
    # INITIALIZATION OF VELOCITY
    # ---------------------------
    # velocity between -0,5 and 0,5
    x_sum, y_sum, z_sum = 0, 0, 0
    for i in range(s_part+p_part):
        vel_x[i]= random.uniform(-0.5, 0.5)
        sumv2 += vel_x[i] ** 2
        sumv += vel_x[i]
        x_sum += x[i]


    sumv = sumv * n_part_inv
    sumv2 = sumv2 * n_part_inv
    fs = np.sqrt(3 * T / sumv2)
    #ec_list.append(sumv2/2)

    for i in range(n_part):
        vel_x[i]= (vel_x[i] - sumv) * fs
        x[i]= x[i] - x_sum * n_part_inv


## PUT THIS SUBROUTINES AT THE BEGINNING OF THE CODE AND CREATE A CLASS TO AFFORD ERRORS ##
    def subroutine_force(en):

        global r,s_part, e_part,forward_reaction_time,backward_reaction_time,forward_reaction,backward_reaction,r,over_l_r,time, well_counts, First_well,first,second
        first = False
        second = False
        en = 0

        for k in range(n_part):
            prev_f_x[k], f_x[k] = f_x[k], 0.0.  # save old forces, set to zero new ones

        #Interaction between S and E particles#
        if e_part == 2:
            x_r1 = x[0] - x[2] # S0 interacts with C2
            x_r2 = x[1] - x[3] # S1 interacts with C3
            x_r1 = x_r1 - box_width * round(x_r1 / box_width) # BC
            x_r2 = x_r2 - box_width * round(x_r2 / box_width) # BC

            if abs(x_r1) < z_i :
                ff = -k_i * (x_r1 - l_i)
                f_x[0] += ff
                f_x[2] -= ff
                en += k_i_div_2 * (abs(x_r1) - l_i)**2 - e_zi
                first = True

            if abs(x_r2) < z_i :
                ff = -k_i * (x_r2 - l_i)
                f_x[1] += ff
                f_x[3] -= ff
                en += k_i_div_2 * (abs(x_r2) - l_i)**2 - e_zi
                second = True

            if second and first:
                print("2 interactions  ",time)

        ## Force between S pairs (1 pair = 1 S molecule) ##
        #
        r = x[0] - x[1]
        #r = x[0]
        #r = r - box_width * round(r / box_width)
        #print("loop " + str(r))


        if abs(r) <= z_a:
            ff = - k_ar * (r - l_ar)  # Calculation of the force
            f_x[0] -= ff  # Incrementation of the force for S1
            f_x[1] += ff  # Incrementation of the force for S2
            en += k_ar_div_2 * (abs(r) - l_ar)**2 - e_za  # Calculation of the energy
            type[0], type[1] = 'S_0', 'S_1' # Change the type if needed


        ## Force between P pairs (1 pair = 1 S molecule) ##
        else:
            ff = - k_r * (r - l_r)  # Calculation of the force
            f_x[0] -= ff  # Incrementation of the force for P1
            f_x[1] += ff  # Incrementation of the force for S2
            en += k_r_div_2 * (abs(r) - l_r) ** 2  # Calculation of the energy
            type[0], type[1] = 'P_0', 'P_1' # Change the type if needed




        if abs(r) > 3.1 and not forward_reaction:
            forward_reaction = True
            forward_reaction_time = time


        ## Force between E pairs (1 pair = 1 S molecule) non fixed##
        # if e_part != 0 :
        #     for i in range((s_part + p_part)//2,n_part//2):
        #         #print('E interaction')
        #         E_1 = 2*i
        #         E_2 = 2*i + 1
        #         x_r = x[E_1] - x[E_2]
        #         x_r = x_r - box_width * round(x_r / box_width)
        #         ff = - k_e * (x_r - l_e)
        #         f_x[E_1] = f_x[E_1] + ff
        #         f_x[E_2]= f_x[E_2] - ff
        #         en += 0.5 * k_e * (abs(x_r) - l_e)**2



        for i in range(n_part):
            f_drag_x = -gamma * vel_x[i]  # Dissipative Force
            #f_drag_y = -gamma * vel_y[i]
            #f_drag_z = -gamma * vel_z[i]
            f_brown_x = sigma_brown_divided_by_sqrt_dt * np.random.normal(0, 1)  # Brownian motion
            f_x[i] += f_brown_x
            f_x[i] += f_drag_x



        return en


    def subroutine_integrate(en):
        sumv = 0.0
        sumv2 = 0.0
        sumv2_x = 0.0
        ########################
        ## Verlet improved ##
        #########################
        vel_x_intermediate = vel_x[:]
        for index in range(2):
            x[index]= x[index] + vel_x[index] * dt + dt_2_over_2_mass * f_x[index]
            #x[index]= x[index] - round(x[index] / box_width) * box_width
            vel_x[index] += lambda_constant * dt * f_x[index]

        # if catalyst is rigid #
        # if e_part == 2:
        #     x[3] = x[2] + l_e
        #     x[3] = x[3] - round(x[3] / box_width)*box_width

        en = subroutine_force(en)

        for index in range(2):
            vel_x[index]= vel_x_intermediate[index] + (prev_f_x[index] + f_x[index]) * dt_over_2
            # verlet v_n+1 = v_n + 1/2 * dt * (a_n + a_n+1)
            v2 = vel_x[index] ** 2
            sumv += vel_x[index]
            sumv2 += v2





        #########################
        #Euler - Maruyama method
        #########################
        ## USING THE OVERDAMPED LANGEVIN EQUATION ##
        #x(t+∆t)=x(t)+μ(x,t)∆t+σ(x,t)η√∆t
        # en = subroutine_force(en)
        # for index in range(s_part+p_part):
        #     next_x[index] = x[index] + dt * f_x[index]*gamma_inv
        #     next_x[index]= next_x[index] - round(next_x[index] / box_width) * box_width
        #     #print(next_x[index],x[index])
        #     vel_x[index] = (next_x[index] - x[index])*inv_dt
        #     x[index] = next_x[index]
        #     v2 = vel_x[index] ** 2
        #     sumv += vel_x[index]
        #     sumv2 += v2



        T = sumv2 * mass_n_part_inv
        #T = sumv2

        Temp_list.append(T)
        #ec_list.append(0.5 * sumv2 * n_part_inv)
        e_tot = (en + 0.5 * sumv2) * n_part_inv
        #e_tot_list.append(e_tot)
        en = en * n_part_inv
        #ep_list.append(en)
        global time
        time += dt
        time_list.append(time)

        return T, e_tot, en


    en = subroutine_force(en)

    i = 0
    #while (not backward_reaction or not forward_reaction) and i < n_steps:
    #print(x[0],x[1],r)

    # ---------------------------
    #  LOOPING OVER THE N_STEPS
    # ---------------------------
    while i<n_steps:
        if i % 100 == 0:
            #print("hors loop " + str(r))
            r_t_file.write(str(r)+";"+str(time)+";"+str(en)+";"+str(Temp_average)+"\n")
        #r_list.append(str(r))
        #print(str(r))
            #print(i)

        # ---------------------------
        #     OVITO SIMULATION
        # ---------------------------
        # if i % 100 == 0:
        #     with open(str(nom_fichier)+".dat", 'a') as data_particle:
        #         data_particle.write("ITEM: TIMESTEP\n" + str(i) + "\nITEM: NUMBER OF ATOMS\n" + str(n_part) + "\n")
        #         data_particle.write("ITEM: BOX BOUNDS fff\n")
        #         data_particle.write(str(-box_width_over_2) + " " + str(box_width_over_2) + "\n")
        #         data_particle.write(str(-box_width_over_2) + " " + str(box_width_over_2) + "\n")
        #         data_particle.write(str(-box_width_over_2) + " " + str(box_width_over_2) + "\n")
        #         data_particle.write("ITEM: ATOMS radius x y z v_X v_Y v_Z Type\n")
        #         for j in range(n_part):
        #             data_particle.write(
        #                 str(radius) + " " + str(x[j]) + " " + str(y[j]) + " " + str(z[j]) + " " + str(vel_x[j]) + " " + str(
        #                     vel_y[j]) + " " + str(vel_z[j]) + " " + type[j] + "\n")
        #if i % 10000 == 0:
            #print(i)

            # ---------------------------
            # INTEGRATE THE MOTION EQUATION
            # ---------------------------
        T, e_tot, en = subroutine_integrate(en)
        i += 1
        section = Temp_list[:i+1]
        Temp_average = sum(section) / (len(section))
        Temp_average_list.append(Temp_average)


        # ---------------------------
        #       PLOT FUNCTIONS
        # ---------------------------

    #plt.plot(time_list,e_tot_list, label = "e_tot")
    #plt.plot(time_list,ec_list, label = "kinetic energy")
    #plt.plot(time_list,ep_list, label = "potential energy")
    #plt.plot(time_list,Temp_list, label = "Temperature")
    # plt.plot(time_list,Temp_average_x, label = "X Average Temperature")
    # plt.plot(time_list,Temp_average_y, label = "Y Average Temperature")
    # plt.plot(time_list,Temp_average_z, label = "Z Average Temperature")
    #plt.plot(time_list,Temp_average_list, label = "Average Temperature")
    #plt.plot(time_list,s_part_list, label = "S_part")
    #plt.plot(time_list,p_part_list, label = "P_part")
    # plt.plot(time_list,f_drag_z_list, label = "Z drag force")
    # plt.plot(time_list,f_drag_y_list, label = "Y drag force")
    # plt.plot(time_list,f_drag_x_list, label = "X drag force")
    # plt.plot(time_list,f_brown_z_list, label = "Z brown force")
    # plt.plot(time_list,f_brown_y_list, label = "Y brown force")
    # plt.plot(time_list,f_brown_x_list, label = "X brown force")
    #plt.plot(time_list,r_list,label = "distance")
    #plt.plot(r_list,ep_list, label = "potential energy")
    #plt.scatter(r_list,ep_list, label = "potential energy")
    #plt.legend()
    #plt.show()
    #plt.savefig(str(nom_fichier)+"_plotfig_ep(r)"+".pdf")
    # plt.xlim(0,8)
    #plt.ylim(0, n_part + 1)

    #ratio_file.write(";"+str(forward_reaction_time)+"\n")

#ratio_file.close()

    r_t_file.close()
