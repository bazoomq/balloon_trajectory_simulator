# /******************************************************************************

# This code simulates weather balloon ascent using the data of Standard Atmosphere. 

# *******************************************************************************/

from math import pow,sqrt,exp,atan,cos
import datetime
from astra_GFS import GFS_Handler
import pandas as pd 


#  Standard Atmosphere

Hatm = [
  0,
  500,
  1000,
  1500,
  2000,
  2500,
  3000,
  4000,
  5000,
  6000,
  7000,
  8000,
  9000,
  10000,
  11000,
  12000,
  14000,
  16000,
  18000,
  20000,
  24000,
  28000,
  32000,
  36000
 ]
Patm = [
  101330,
  95464,
  89877,
  84559,
  79499,
  74690,
  70123,
  61661,
  54052,
  47217,
  41106,
  35653,
  30801,
  26500,
  22700,
  19399,
  14170,
  10353,
  7565,
  5529,
  2971,
  1616,
  889,
  499
 ]
Tatm = [
  288.2,
  284.9,
  281.7,
  278.4,
  275.2,
  271.9,
  268.7,
  262.2,
  255.7,
  249.2,
  242.7,
  236.2,
  292.7,
  223.3,
  216.8,
  216.7,
  216.7,
  216.7,
  216.7,
  216.7,
  220.6,
  224.5,
  228.5,
  239.3
 ]
Beta = [
  0.000119267,
  0.000120614,
  0.000121985,
  0.00012341,
  0.000124796,
  0.000126191,
  0.000128599,
  0.000131705,
  0.000135193,
  0.0001386,
  0.000142321,
  0.000146286,
  0.000150402,
  0.00015478,
  0.000157143,
  0.000157047,
  0.000156925,
  0.000156872,
  0.000156763,
  0.000155277,
  0.000152236,
  0.000149403,
  0.000144373,
  0
 ]


# Standard Atmosphere 

# Input data 

def main():
    pi   = 3.1415926535    # Pi
    R_earth = 6371009      # m
    R_id = 8314.462        # gas ant                                   (J/K/mol*1000)
    pix  = pi / 180        # converts degree to radian
    pix_rev = 180 / pi
    Ksph = 4.0 * pi / 3    # sphere volume coefficient 
    g    = 9.80665         # free fall acceleration at release location     (m/s^2)
    k_Sp = 0.739668778     # coefficient of pumpkin shape area
    k_Vp = 1.21852421611857# coefficient of pumpkin shape volume
    k_S  = k_Sp * sqrt(pi) # coefficient of pumpkin shape area
    k_V  = 2.74581225      # coefficient of pumpkin shape volume
    k_p_rate = (5 / 3) / 100000 # converts pumpimg rate unit from l/min to m^3/s

    # Balloon and ballonet data
    ro_f = 920                   # density of the balloon film                    (kg/m^3)
    d_f  = 0.0000254             # thickness of the balloon film                  (m)
    R_b  = 6.122831782671        # radius of fully the inflated balloon pumpkin   (m)
    R_bl = 2                     # ballonet pumpkin shape radius                  (m)
                                 # 1/2 of the ballonet pumpkin meridian length    (m)
    S_b  = pi * R_b * pi * R_b     # balloon surface area                           (m^2)
    S_bl = pi * R_bl * pi * R_bl   # ballonet surface area                          (m^2)
    # Sb_drag                      # balloon drag area ( = pi*Rx^2 )                (m^2)
    # V_b                          # balloon volume                                 (m^3)
    # V_b_max                      # volume of fully inflated balloon pumpkin       (m^3)
    # V_bl                         # ballonet volume                                (m^3)
    # V_bl_max                     # volume of fully inflated ballonet pumpkin      (m^3)
    # m_b                          # balloon mass                                   (kg)
    # m_bl                         # ballonet mass                                  (kg)
    # Rx                           # current maximum radius of the partly inflated balloon's pumpkin-like top part (m) 
    # x_cone                       # radius of the bottom cone-shaped part of the partly inflated balloon          (m)
    # h_cone                       # height of the bottom cone-shaped part of the partly inflated balloon          (m)
    # l_cone                       # slant height of the bottom cone-shaped part of the partly inflated balloon    (m)
    # Cd_z_up                      # drag cefficient for ascend
    # Cd_z_down                    # drag cefficient for descent
    # Cd_z                         # drag cefficient along z_0
    Cd_z_max = 0.47                # max drag cefficient along z_0
    # Cd_xy                        # drag cefficient for horizontal motion
    # eps                          # = x/Rx 
    # eps3                         # = eps^3
    # eps0                         # lower limit of the search range for eps        (rad)    
    # eps1                         # uper  limit of the search range for eps        (rad) 
    # gama                         # semi-vertex angle of the bottom cone-shaped 
                                   #   part of the partly inflated balloon          (rad)
    gama_max = 31.4 * pix          # max gamma, beyond which Cd_z = 0.47            (deg*pix -> rad)
    eps_error = 0.0001             # accuracy of defining eps                       (%)
    # f_S                          # 2*So_b = f_S(eps)*R_b
    # f_V                          # V_b = f_V(eps)*R_b^3
    # V_b_eq                       # equation for eps
    # P_bl                         # air pressure in the ballonet                   (Pa)
    # dP                           # dP = P_gas - P_atm  balloon  diff pressure    (Pa)
    # dP_bl                        # dP = P_bl  - P_gas  ballonet diff pressure    (Pa)    
    dP_max    =  100               # balloon  max diff pressure = P_b-P_atm         (Pa)
    dP_max_bl = 7500               # ballonet max diff pressure = P_bl-P_gas        (Pa)
    m_bl_air   = 0.0               # current air mass in the ballonet               (kg)   0.9
    # m_bl_air_1                   # ballonet air mass at which P_gas = Patm        (kg)
    # m_bl_air_2                   # ballonet air mass at which P_bl = P_gas_max    (kg)
    # V_bl_air_1                   # ballonet volume at m_bl_air = m_bl_air_1       (m^3)
    # T_bl_air                     # air temperature in the ballonet                (K)
    # K_b_bl                       # auxiliary coefficient for defining V_bl
    diff_Pmode = 0                 # differential pressure mode values: 1/2/3/4
    # diff_1                       # flag for diff pressure mode 1  
    # diff_2                       # flag for diff pressure mode 2  
    # diff_3                       # flag for diff pressure mode 3
    # diff_4                       # flag for diff pressure mode 4
    # Covering net data (net with diamond-shaped cells)
    F_tendon_b =  4                # breaking strength of the balloon  tendon       (kg)
    F_tendon_bl = 20               # breaking strength of the ballonet tendon       (kg)
    k_sf       =  2                # safety factor
    # L_h                          # horizontal size of a mesh                      (m)
    # L_v                          # verticlal  size of a mesh                      (m)
    # n_h                          # number of cells along the equator
    # n_v                          # number of cells along the meridian    
    Ad         = 0.38              # aspect ratio of a diamond-shaped mesh (L_h/L_v)
    # F_mesh                       # force on the mesh film caused by diff pressure (kg)
    # n_meshes                     # number of net mehses
    # m_tendons                    # total mass of tendons of the net               (kg)

    # Payload data
    m_p  = 8                     # payload                                        (kg)
    S_p  = 0.16                  # payload drag area                              (m^2)
    Cd_p = 0.25                  # payload drag cefficient

    # Pump data
    N_pump          = 6             # number of pumps                                
    N_pumps_on      = 3             # number of pumps running at a time
    m_pump       = 0.5            # weight of one pump                             (kg)

    #     P_diff_pump  = 8000          # maximum differential pressure of the pump      (Pa)
    # m_pump_all                   # weight of all pumps                            (kg)
    pump_rate_V  = 600           # volumetric pumping rate of one pump            (l/min)
    pumping_rate_V = N_pumps_on * pump_rate_V * 0.001/60                
                                        # overall volumetric pumping rate                (m^3/s)
    # pumping_rate_m               # overall air mass pumping rate                  (kg/s)
    k_pump       = 0.0           # pump gain - due to more powerfull BLDC motor
    S_pump       = pi * 0.01 * 0.01  # pump outlet cross section area                 (m^2)
    pump_on_off     = 0             # 0 - off, 1 - on

    # LTA gas data
    # mu_gas                        # molar mass of the LTA gas                      (gram)
    k_He_H2 = 1                     # LTA gas: 1 - helium , 0 - hydrogen
    # xmu_gas                       # R/mu ratio for the LTA gas 
                                        # ( P = (ro_gas/mu)*R*T = (R/mu)*ro_gas*T )
    m_gas_0 = 0.569                   # initial LTA gass mass                          (kg)    # = 0.569
    # m_gas                        # working LTA gass mass                          (kg)    #  = 3.414515
    dT_gas  = 0.0                  # additional temperature of the LTA gas due to greenhouse effect
    # V_gas                        # volume of the LTA gas                          (m^3)
    # P_gas                        # pressure of the LTA gas                        (Pa)
    # T_gas                        # temperature of the LTA gas                     (Pa)
    # P_gas_max                    # max pressure of the LTA gas at V_bl = V_bl_max (Pa)

    # Atmospheric data and variables
    mu_air  = 28.966             # air molar mass                                 (gram)
    xmu_air = R_id / mu_air      # mu/R ratio for air    (  P = (rho/mu)RT  )
    # P_atm                        # air pressure at the current altitude           (Pa)
    # T_atm                        # air temperature at the current altitude        (K)
    # rho_atm                      # air density at the current altitude            (kg/m^3)
    # P_atm_hmax                   # air pressure at maximum altitude               (Pa)
    # T_atm_hmax                   # air temperature at at maximum altitude         (K)
    Vw_xy0  = 0                    # initial horizontal velocity of wind            (m/s)
    #   Vw_z0   = 0                # initial vertical velocity of wind              (m/s)
    #   Vw_z                       # current vertical velocity of wind              (m/s)
    # Vw_xy                        # current horizontal velocity of wind            (m/s)  
    dA_dh   = 0.4*pix              # wind azimuth change per km of altitude         (deg/km*pix -> rad/km)
    A0      = 50*pix               # initial wind azimuth at h0 (ref: North vector) (deg*pix -> rad)
    # A                            # current wind azimuth at h0 (ref: North vector) (rad)
    dVw_dh  = 0.77                 # change of Vw_xy per km of altitude increase    (m/s/km)        

    # Overload    
    # Pg0                          # previuos overload                              (g)
    # Pg                           # current overload                               (g)
    # Pg_max                       # maximum overload                               (g)
    # mg                           # current overload force                         (kg)
    # mg_max                       # maximum overload force                         (kg)

    # Declaration and initialization of variables 
    launch_mode = False            # True, if in launch mode, False - otherwise 
    pumping_on_off = False         # pumping mode: False - disabled, True - enabled
    t0   = 0                       # previous time                                  (s)
    t    = 0                       # current time                                   (s)
    # ad_x                         # air drag acceleration along x                  (m/s^2)
    # ad_y                         # air drag acceleration along y                  (m/s^2)
    # ad_z                         # air drag acceleration along z                  (m/s^2)
    vx_0 = 0                       # previous velocity along x                      (m/s)
    vy_0 = 0                       # previous velocity along y                      (m/s) 
    vz_0 = 0                       # previous velocity along z                      (m/s) 
    vx   = 0                       # current velocity along x                       (m/s)
    vy   = 0                       # current velocity along y                       (m/s) 
    vz   = 0                       # current velocity along z                       (m/s) 
    v_0  = 0                       # previous velocity                              (m/s)
    v    = 0                       # current velocity                               (m/s)
    h0   = 1000                    # 1120  24000  # initial altitude                 ( < 36000 )   (m)
    h_max= 25000                   # the max altitude to stop ascent  ( < 36000 )   (km)
    h_min= 0                       # the min altitude to stop descent ( < 36000 )   (km)
    x_0  = 0                       # previous x coordinate (East)                   (m)
    y_0  = 0                       # previous y coordinate (North)                  (m) 
    z_0  = h0                      # initial/previous z coordinate (altitude)       (m) 
    x    = 0                       # current x coordinate  ( -> East  )             (m)
    y    = 0                       # current y coordinate  ( -> North )             (m) 
    z    = h0                      # current z coordinate  (altitude)               (m) 
    r    = 0                       # current distance                               (m)
    # a_x                          # acceleration along x                           (m/s^2)
    # a_y                          # acceleration along y                           (m/s^2)
    # a_z                          # acceleration along z                           (m/s^2)
    # mu_                          # gas/air molar mass ratio
    dt   = 3                    # 0.001; numerical integration step                     (s)
    dt_pr = 60                      # cout step                                      (s)
    # Fa                           # buoyant force                                  (N)
    # Fg                           # gravitational force                            (N)
    # Fd_z                         # air drag force for vertical motion             (N)
    # m_net                        # overall net mass                   (kg)
    # m                            # overall mass                       (kg)
    dvx      = 0                   # x increnment                       (m)
    dvy      = 0                   # y increnment                       (m)
    dvz      = 0                   # z increnment (altitude)            (m)
    dv       = 0                   #  velocity increnment                (m/s)
    dx       = 0                   # x increnment                       (m)
    dy       = 0                   # y increnment                       (m)
    dz       = 0                   # z increnment (altitude)            (m)
    dh_plus  = 1000
    # n_osc       = 0
    # n                               # number of atmospհeric layers in the standard model 
    # i                               # index of the current atmospհeric layer
    # i_max                           # index of the atmospհeric layer for the max altitude
    j           = 1


    # # for Runge-Kutta method
    # u[5]    = 0,0,1./2,1./2,1
    # k_vy[5] = 0
    # k_vx[5] = 0
    # k_vz[5] = 0
    # k_y[5]  = 0
    # k_x[5]  = 0
    # k_z[5]  = 0
    # Vx_temp[5] = 0
    #      Vy_temp[5] = 0
    #      t_temp[5]  = 0
    #      x_temp[5]  = 0
    #      y_prev     =  0
    #      y_temp[5]  = 0
    #      m_temp[5]  = 0
    #      ms_temp[5] = 0
    #      W_temp[5]  = 0
    #      F_temp[5]  = 0
    #      k_Imp[5]   = 0

    Vw_xy = Vw_xy0                    
    A     = A0        
    So_b  = k_S * R_b 
    So_bl = k_S * R_bl 
    V_b_max  = k_Vp * pow(So_b, 3)
    V_bl_max = k_Vp * pow(So_bl, 3) 
    m_b  = ro_f * S_b  * d_f
    m_bl = ro_f * S_bl * d_f
    m_pump_all    = N_pump * m_pump
    m_net    = m_p + m_b + m_bl + m_pump_all # overall net mass                   (kg)
    
    if k_He_H2 == 1:
        mu_gas = 4      
    else:
        mu_gas = 2
    
    
    latitude_start, longitude_start = 40.238593, 44.505821
    lon = longitude_start
    lat = latitude_start
    
    # min_longitude, max_latitude = [longitude_start, longitude_start + 0.2], [latitude_start, latitude_start + 0.2]
    UTC_offset = 4
    myGFSlink = GFS_Handler(longitude_start, latitude_start, datetime.datetime.now() - datetime.timedelta(seconds=UTC_offset * 3600)-datetime.timedelta(days=7))
    myGFSlink.downloadForecast()    
    temperature, pressure, u_wind, v_wind = myGFSlink.interpolateData('t', 'p', 'd','s')
    
    mu_  = mu_gas/mu_air
    xmu_gas = R_id / mu_gas
    for n in range(0, 23):       
        if (h_max > Hatm[n]):
            i_max = n
        else:
            break
            
        
    
    T_atm_hmax = Tatm[i_max] + (h_max - Hatm[i_max]) * (Tatm[i_max + 1] - Tatm[i_max]) / (Hatm[i_max + 1] - Hatm[i_max])    
    P_atm_hmax = Patm[i_max] * exp(-Beta[i_max] * (h_max - Hatm[i_max]))
    # /*  cout << "i_max=" << i_max <<endl
    #     cout << "T_atm_hmax=" << T_atm_hmax <<endl
    #     cout << "P_atm_hmax=" << P_atm_hmax <<endl
    #     cout << "Hatm[i_max]=" << Hatm[i_max] <<endl
    #     cout << "Hatm[i_max+1]=" << Hatm[i_max+1] <<endl
    #     cout << "Tatm[i_max]=" << Tatm[i_max] <<endl
    #     cout << "z=" << z <<endl
    # */
    
    latlon_to_csv = pd.DataFrame(columns=['Longitude', 'Latitude', 'Altitude'])
    flight_info = pd.DataFrame(columns=['Altitude', 'vx', 'vy', 'mgas'])

    while (z <= h_max - 10000 and z > 0): # h_max + dh_plus
        
    #while (n_osc < 1 and z > 0)
    # TODO Check Hatm size it should be 24
        for n  in range(0, 24):
        
          if (z > Hatm[n]):
            
            i = n
           
        #T_atm = Tatm[i] + (z - Hatm[i]) * (Tatm[i + 1] - Tatm[i]) / (Hatm[i + 1] - Hatm[i])
        #P_atm = Patm[i] * exp(-Beta[i] * (z - Hatm[i]))
        T_atm = temperature(lat, lon, z, myGFSlink.getGFStime(datetime.datetime.now() - datetime.timedelta(seconds=UTC_offset * 3600)-datetime.timedelta(days=7)+datetime.timedelta(seconds=t))) 
        P_atm = pressure(lat, lon, z, myGFSlink.getGFStime(datetime.datetime.now() - datetime.timedelta(seconds=UTC_offset * 3600)-datetime.timedelta(days=7)+datetime.timedelta(seconds=t))) * 100
        
        rho_atm = P_atm / xmu_air / T_atm
        T_gas = T_atm + dT_gas
        
        if (z < h_min and launch_mode):
            
            m_gas = m_gas_0
            
        else:
            
            m_gas = V_b_max * P_atm_hmax / xmu_gas / (T_atm_hmax+dT_gas)
            
        
        P_gas_max = m_gas * xmu_gas * T_gas / (V_b_max - V_bl_max)
        
        T_bl_air   = T_atm

        m_bl_air_1 = P_atm * (mu_air/R_id/T_atm) * (V_b_max-(m_gas/mu_gas)*R_id*(T_atm+dT_gas)/P_atm)

        if(P_atm < P_gas_max):
            
            m_bl_air_2 = (m_gas / mu_) * (V_bl_max / (V_b_max - V_bl_max)) * T_gas / T_bl_air
        
        else:
            
            m_bl_air_2 = P_atm * V_bl_max / xmu_air / T_bl_air
        
        #cout << "m_bl_air_1=" << m_bl_air_1 <<endl
        #cout << "m_bl_air_2=" << m_bl_air_2 <<endl
        # /*
        # if (z > h_max) 
        #     cout <<" m_bl_air = "   
        #     cin >> m_bl_air   
        #     
        # */
        K_b_bl = (m_bl_air/m_gas)*(mu_gas/mu_air)*(T_atm/T_gas)
        
        diff_1 = False
        diff_2 = False
        diff_3 = False
        diff_4 = False
       
       #  IF(OR(AND(P_atm<P_gas_max,m_bl_air<mI_bl_air),AND(P_atm>=P_gas_max,m_bl_air<=mII_bl_air),AND(P_atm<P_gas_max,m_bl_air=mII_bl_air)),1,0)
       
        if (((P_atm < P_gas_max) and (m_bl_air <= m_bl_air_1)) or ((P_atm >= P_gas_max) and (m_bl_air <= m_bl_air_2)) or ((P_atm < P_gas_max) and (m_bl_air == m_bl_air_2))):
               
            diff_Pmode = 1
            diff_1 = True
         
            
        #  IF(AND(P_atm>=P_gas_max,m_bl_air>mII_bl_air),1,0)
            
        if (((P_atm >= P_gas_max) and (m_bl_air > m_bl_air_2))):
               
            diff_Pmode = 2
            diff_2 = True
         
            
        #  IF(AND(P_atm<P_gas_max,(m_bl_air-mI_bl_air)>0,(m_bl_air-mII_bl_air)<=0),1,0)
            
        if (((P_atm < P_gas_max) and ((m_bl_air - m_bl_air_1) > 0)) and ((m_bl_air-m_bl_air_2) <= 0)):
               
            diff_Pmode = 3
            diff_3 = True
          
            
        #  IF(AND(P_atm<P_gas_max,(m_bl_air-mI_bl_air)>0,(m_bl_air-mII_bl_air)<=0),1,0)
            
        if ((P_atm < P_gas_max) and (m_bl_air > m_bl_air_2)):
               
            diff_Pmode = 4
            diff_4 = True
          
            
        if (diff_1 + diff_2 + diff_3 + diff_4 != 1):
             
            print( "j=" , j )
            print( "P_atm=" , P_atm )
            print( "T_atm=" , T_atm )
            print( "P_gas_max=" , P_gas_max )
            print( "mu_=" , mu_ )
            print( "m_gas=" , m_gas )
            print( "V_b_max=" , V_b_max )
            print( "P_gas_max=" , P_gas_max )
            print( "m_bl_air=" , m_bl_air )
            print( "m_bl_air_1=" , m_bl_air_1 )
            print( "m_bl_air_2=" , m_bl_air_2 )
            print( "diff_Pmode=" , diff_Pmode )
            print( "diff_1=" , diff_1 )
            print( "diff_2=" , diff_2 )
            print( "diff_3=" , diff_3 )
            print( "diff_4=" , diff_4 )
            print( "pumping_rate_m=" , pumping_rate_m )
            print( "Mode definition failure")
            return
          

        if (diff_1):
                     
            V_bl  = (m_bl_air / mu_air) * (R_id *T_atm / P_atm)
            V_b   = ((m_gas / mu_gas) * (R_id * T_gas / P_atm) + V_bl) #TODO Check order of V_bl and V_b
            P_gas = P_atm
            P_bl  = P_atm
        
            
        if (diff_2):
             
            V_b   = ((m_gas / mu_gas) * (R_id * T_gas / P_atm) + V_bl_max)
            V_bl = V_bl_max
            P_gas = P_atm
            P_bl  = (m_bl_air / mu_air) * (R_id * T_atm / V_bl_max)
        
        
        if (diff_3):
             
            V_b   = V_b_max
            V_bl  = V_b_max * K_b_bl / (K_b_bl + 1)
            P_gas = (m_gas / mu_gas) * R_id * T_gas / (V_b_max - V_bl)
            if (V_bl == 0):
                P_bl = P_atm
            else: 
                P_bl = (m_bl_air / mu_air) * (R_id * T_atm / V_bl)
        

        if (diff_4):
             
            V_b   = V_b_max
            V_bl  = V_bl_max
            P_gas = P_gas_max
            P_bl  = (m_bl_air / mu_air) * R_id * T_atm / V_bl_max
        
            
        dP = P_gas - P_atm
        # /*
        # if (dP > max_dP) 
        #     print( "Balloon burst! Q -> @"
        #     exit(0)
        # 
        # */
        
        dP_bl = P_bl  - P_gas
        
        # if (dP_bl > max_dP_bl):
             
        #     print( "Ballonet burst! O -> @")
        #     exit(0)
        
           
           
        # Solve the following equation to get epsilon:
        # V_b = (pi/3)*(2+1/tan(eta)*(4/(pi+2/sin(eta)))^3*So^3
        eps0 = 0
        eps1 = 1
        while ((abs((eps1-eps0)/eps1))*100 >= eps_error):
            
            eps    = (eps0 + eps1) / 2
            eps3   = eps*eps*eps
            f_V    = (pi/3)*eps3*eps*eps/sqrt(1-eps*eps3) + k_V*(1-eps3/2)
            f_S    = eps/sqrt(1-eps*eps3) + k_S*(2-eps)
            #f_V    = (pi/3)*(1/eps3)/sqrt(eps*eps3-1) + k_V*(1-1/eps3/2)
            #f_S    = 1/eps/sqrt(eps*eps3-1) + k_S*(2-eps)
            V_b_eq = f_V*pow((2*So_b/f_S),3) - V_b
            
            # cout << "eps=" << eps <<endl
            if (V_b_eq < 0):
                eps1 = eps
            else:
                eps0 = eps
            # cout << V_b_eq << "\t" << eps << "\t" << eps0 << "\t" << eps1 <<endl
        
        vx = u_wind(lat, lon, z, myGFSlink.getGFStime(datetime.datetime.now() - datetime.timedelta(seconds=UTC_offset * 3600)-datetime.timedelta(days=7)+datetime.timedelta(seconds=t)))
        vy = v_wind(lat, lon, z, myGFSlink.getGFStime(datetime.datetime.now() - datetime.timedelta(seconds=UTC_offset * 3600)-datetime.timedelta(days=7)+datetime.timedelta(seconds=t)))
        
        Rx     = 2 * So_b / f_S 
        x_cone = eps * Rx
        gama   = atan(sqrt(1/(eps*eps3)-1))
        h_cone = x_cone/sqrt(1/(eps*eps3)-1)
        l_cone = sqrt(x_cone*x_cone+h_cone*h_cone)
        if (Rx < R_bl):
            Rx = R_bl  
        Sb_drag  = pi * Rx * Rx
              
        if (vz < 0):
            
            #Cd_z = Cd_z_max
            if (gama <= gama_max):
                
                Cd_z = 1 - pow(cos(gama),4) 
                #Cd_z = Cd_z_max
            
            else:
                Cd_z = Cd_z_max    
        
        else:
            Cd_z = Cd_z_max
        
        m    = m_net + m_gas + m_bl_air
        Fa   = rho_atm * V_b * g
        Fg   = m * g
        Fd_z = Cd_z * (rho_atm * vz * abs(vz) / 2) * Sb_drag
    
        a_z  = (Fa - Fg - Fd_z) / m
        # if t < 10:
            
        #     print(t, " ", a_z)
        
        dvz  = a_z * dt
        vz_0 = vz
        vz   = vz + dvz
        #if(vz*vz0 < 0) n_osc++
        dz   = ((vz + vz_0) / 2) * dt
        z    = z + dz
        
        dx = vx * dt
        dy = vy * dt

        # # from lat lon to x y
        # x = x_0 + (R_earth + z) * cos(latitude_start * pix) * (lon * pix - longitude_start * pix)
        # y = y_0 + (R_earth + z) * (lat * pix - latitude_start * pix)
        
        
        
        x = x + dx # u - vx 
        y = y + dy # v - vy    
        
        # from x y to  lon lat
        lon =  pix_rev * (longitude_start * pix + dx / (R_earth + z) / cos(latitude_start * pix))
        lat = pix_rev * (latitude_start * pix + dy / (R_earth + z))
        
        longitude_start = lon
        latitude_start = lat
        
        if ((P_bl - P_gas < dP_max_bl) and pumping_on_off):
            
            pumping_rate_m = N_pumps_on * sqrt(rho_atm * (rho_atm * pumping_rate_V*pumping_rate_V - 2*dP_bl*S_pump*S_pump/k_pump/k_pump))        
            m_bl_air  = m_bl_air + pumping_rate_m * dt
            pump_on_off = 1
        
        else:
            pump_on_off = 0
        
        
        # if ((z < h_target) and target_on_off):
        #     vent_on_off = 1
        #     venting_rate_m = ro_bl_air * S_vent 
            # pump_on_off * N_pumps_on * sqrt(rho_atm * (rho_atm * pumping_rate_V*pumping_rate_V - 2*(P_atm-P_bl)*S_pump*S_pump/k_pump/k_pump));        
            # m_bl_air  = m_bl_air + p
        Pg   = a_z / g
        mg   = Pg * m
        
        t = t + dt

        if ((t - t0) >= dt_pr):
            
            t0 = t
            #print(t , "\t" , z , "\t" , vz , "\t" , a_z , "\t" , Pg , "\t" , P_atm , "\t" , T_atm , "\t" , rho_atm , "\t" , m_bl_air , "\t" , diff_Pmode ,"\t" , dP , "\t" , dP_bl , "\t" , m_gas , "\t" , gama/pix ,  "\t" , Cd_z ,  "\t" , m ,  "\t" , pump_on_off , "\t" , Rx)
            print("Current Alt = ", z)
            print("Current Lon = ", lon)
            print("Current Lat = ", lat)
            new_row = {'Longitude':  lon, 'Latitude': lat, 'Altitude': z}
            info_row = {'Altitude': z, 'vx': vx, 'vy': vy, 'mgas': m_gas}
            latlon_to_csv = latlon_to_csv.append(new_row, ignore_index=True)
            flight_info = flight_info.append(info_row, ignore_index=True)
    
    latlon_to_csv.to_csv("flight2.csv", index=False)
    flight_info.to_csv("flight_info2.csv", index=False)
    return

if __name__ == "__main__":
    main()



# V

# 40.23822423742247,44.5079306954846
# 40.238809862507416,44.51110169829215
# 40.240032916754764,44.51530271713292
# 40.24146523010262,44.51938411995484
# 40.24372178933436,44.522827890893105