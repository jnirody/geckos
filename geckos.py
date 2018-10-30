'''
geckos.py 
Graphing and statistically analyzing gecko water walking
created: April 12, 2017
last updated: Aug 21, 2017
Judy Jinn
'''

'''                         #######################
#------------------------   ## ---   Set up  --- ##     ------------------------
                            #######################
'''
"""import every package you'll need"""
import os 
import re #used for natural sort function

import pandas as pd  #data frames
import math

import matplotlib
matplotlib.use("TkAgg") # Need to do this for tkinter to properly import and allow simulataneous matplotlib figures
from matplotlib import pyplot as plt
import matplotlib.animation as animation #animate plots, unused
import matplotlib.cm as cm
import matplotlib.patches as patches
import tkinter as Tk #opening files
from tkinter import filedialog
import numpy as np
import scipy.stats as stats
import scipy.io as sio
import scipy.stats as stats

# auto-detect width of terminal for displaying dataframes
pd.set_option('display.max_columns',0) 
plt.ion() # Solves that matplotlib hang problem
pd.options.mode.chained_assignment = None 


'''                                         #######################
#----------------------------------------   ## ---    MAIN   --- ##     -----------------------------------------
                                            #######################
'''

if __name__ == '__main__':
    path = os.getcwd()
    geckos_csv = pd.read_csv(path+"/scaled values.csv")
        
    soap_water_avg_period = geckos_csv.groupby('expt').mean()['period'].reset_index()
    geckos_csv['expt_avg_period'] = geckos_csv['expt']
    geckos_csv['expt_avg_period'] = geckos_csv['expt_avg_period'].map(
        {'soap': 0.104583, 'water': 0.088179})
    
    # for weird numbers that result if we use every individual period time
        # 0.09097 is average period for all (soap + water)
        # geckos_csv['period'] will use individual's periods
        # geckos_csv['expt_avg_period'] self explanatory   
    # Divide all by 2 for single stride
    stride = geckos_csv['period'] / 2
    geckos_csv['stride'] = geckos_csv['period'] / 2

    # Imput gecko weights
    geckos_csv['weight_g'] = geckos_csv['geckoID']    
    geckos_csv['weight_g'] = geckos_csv['weight_g'].map(
        {3: 5.4, 5: 5.8, 6: 4.89, 7: 5.71, 8: 6.9 , 9: 6.6, 10: 5.3}
        )
    

    
    # Transform weights in gram to mN and impulses per stride
    geckos_csv['weight_mN'] = geckos_csv['weight_g']/1000.0*9.8*1000
    
    geckos_csv['weight_imp'] = geckos_csv['weight_mN'] * stride


    
    # See excel sheet of equations for calculations
    # Finding the impulse for entire cycle

    # Gecko body measurements
    geckos_csv['fl_mm'] = geckos_csv['geckoID']   
    geckos_csv['fl_mm'] = geckos_csv['fl_mm'].map(
        {3: 19.9, 5: 19.3, 6: 19.3, 7: 21, 8: 19,\
        9: 17.5, 10: 19.3}
        ) 
    geckos_csv['ffoot_mm'] = geckos_csv['geckoID']   
    geckos_csv['ffoot_mm'] = geckos_csv['ffoot_mm'].map(
        {3: 6.2, 5: 6.65, 6: 7.0, 7: 6.5, 8: 6.2,\
        9: 7.4, 10: 6.6}
        )
    geckos_csv['bl_mm'] = geckos_csv['geckoID']   
    geckos_csv['bl_mm'] = geckos_csv['bl_mm'].map(
        {3: 24.1, 5: 23.38, 6: 22.8, 7: 23.0, 8: 19.2,\
        9: 26.9, 10: 24.3}
        )
    geckos_csv['bfoot_mm'] = geckos_csv['geckoID']   
    geckos_csv['bfoot_mm'] = geckos_csv['bfoot_mm'].map(
        {3: 7.8, 5: 8.583, 6: 8.1, 7: 9.5, 8: 7.9,\
        9: 10.3, 10: 7.9}
        )
    geckos_csv['tail_mm'] = geckos_csv['geckoID']   
    geckos_csv['tail_mm'] = geckos_csv['tail_mm'].map(
        {3: 45.65, 5: 52.425, 6: 63.8, 7: 39.2, 8: 48.1,\
        9: 58.6, 10: 59.2}
        )    
    geckos_csv['svl_mm'] = geckos_csv['geckoID']   
    geckos_csv['svl_mm'] = geckos_csv['svl_mm'].map(
        {3: 59.91, 5: 57.735, 6: 54.4, 7: 59.6, 8: 57.4,\
        9: 58.6, 10: 56.5}
        )
    geckos_csv['body_ht_mm'] = geckos_csv['geckoID']   
    geckos_csv['body_ht_mm'] = geckos_csv['body_ht_mm'].map(
        {3: 8.1, 5: 10.25, 6: 11.3, 7: 11.9, 8: 12.2,\
        9: 8.6, 10: 9.4}
        )
    
    # Add buoyancy of individual geckos
    rho = 999.97
    v_water = 0.00089
    
    
    # geckos_csv['buoyancy_N'] = geckos_csv['geckoID']
    # geckos_csv['buoyancy_N'] = geckos_csv['buoyancy_N'].map(
    #     {3: 0.023462094, 5: 0.00024, 6: 0.03894067, 7: 0.036983704, 8: 0.044893951,\
    #     9: 0.03387716, 10: 0.031308718}
    #     )
    # geckos_csv['buoyancy_mN'] =  geckos_csv['buoyancy_N']*1000.0

    # buoyancy calculated by modeling body as oval cylinder and tail as oval cone
    h1 = geckos_csv['svl_mm']-(geckos_csv['l_body_above_water']*10)+(0.5* (geckos_csv['body_ht_mm'])/ \
        (np.tan(np.deg2rad(180-geckos_csv['body_angle']))))
    h2 = geckos_csv['svl_mm']-(geckos_csv['l_body_above_water']*10)-(0.5* (geckos_csv['body_ht_mm'])/ \
        (np.tan(np.deg2rad(180-geckos_csv['body_angle']))))
        
    geckos_csv['body_volume'] = np.pi * (geckos_csv['body_ht_mm']/2/1000) * \
        (geckos_csv['body width']/2/1000) * (((h1+h2)/1000)/2)
        
    geckos_csv['tail_volume'] = (1/3)*np.pi*(geckos_csv['body_ht_mm']/2/1000)*(geckos_csv['body width']/2/1000) * \
        (geckos_csv['tail_mm']/1000)
    
    geckos_csv['buoyancy_mN'] = rho * 9.81 * (geckos_csv['body_volume']+geckos_csv['tail_volume'])*1000
    geckos_csv['buoyancy_imp'] = geckos_csv['buoyancy_mN'] * stride
    
    geckos_csv['min imp'] = geckos_csv['weight_mN']*stride
    
    
    # Imput surface tension of each gecko
    geckos_csv['surface_tension_value'] = geckos_csv['expt']
    geckos_csv['surface_tension_value'] = geckos_csv['surface_tension_value'].map(
        {'soap': 0.035, 'water': 0.07197})
    geckos_csv['surface_tension_N'] =  geckos_csv['surface_tension_value'] * np.pi /1000 * (geckos_csv['ffoot_mm']  + geckos_csv['bfoot_mm'])
    geckos_csv['surface_tension_mN'] =  geckos_csv['surface_tension_N']*1000.0
    
    geckos_csv['surface_tension_imp'] = geckos_csv['surface_tension_mN'] * stride
    
    # Calculate foot drag for each gecko

    geckos_csv['fl_drag'] = (rho * \
       (((v_water) * (geckos_csv['fl_mm']/1000))**(1/2)) * \
       ((geckos_csv['avg_strokevel_hor']/1000)**(3/2)) * (geckos_csv['fl_mm']/1000)) + \
       (rho * ((geckos_csv['avg_strokevel_hor']/1000)**2) * (geckos_csv['fl_mm']/1000))

    geckos_csv['bl_drag'] = (rho * \
       (((v_water) * (geckos_csv['bl_mm']/1000))**(1/2)) * \
       ((geckos_csv['avg_bstrokevel_hor']/1000)**(3/2)) * (geckos_csv['bl_mm']/1000)) + \
       (rho * ((geckos_csv['avg_bstrokevel_hor']/1000)**2)  * (geckos_csv['bl_mm']/1000))
       
    # Checking undulation and drag forces and drag
    rho = 999.97
    v_water = 0.00089
    body_len = geckos_csv['svl']/2.0/1000.0
    geckos_csv['body underwater depth'] = geckos_csv['body underwater depth']/100
    geckos_csv['tail underwater depth'] = geckos_csv['tail underwater depth']/100
    geckos_csv['tail length'] = geckos_csv['tail length']/1000
    geckos_csv['vel'] = (geckos_csv['vel']/1000)
    
    geckos_csv['body_thrust'] = rho * ((2*np.pi/(geckos_csv['per bmid']))**2) *      \
        ((geckos_csv['amp bmid'])**2) * body_len * geckos_csv['body underwater depth']
    
    geckos_csv['tail_thrust'] = rho * ((2*np.pi/(geckos_csv['per tmid']))**2) *      \
        ((geckos_csv['amp tmid'])**2) * geckos_csv['tail length'] * geckos_csv['tail underwater depth']
    geckos_csv['bodytail_thrust'] = geckos_csv['body_thrust']+geckos_csv['tail_thrust']
    
    geckos_csv['body_drag'] = (rho * \
       (v_water * (body_len**(1/2))) * \
       (geckos_csv['vel']**(3/2)) * geckos_csv['body underwater depth']) + \
       (rho * (geckos_csv['vel']**2) * body_len* geckos_csv['body underwater depth'])
    geckos_csv['tail_drag'] = (rho * \
       ((v_water * (geckos_csv['tail length']))**(1/2)) * \
       (geckos_csv['vel']**(3/2)) * geckos_csv['tail underwater depth']) + \
       (rho * (geckos_csv['vel']**2) * geckos_csv['tail length']* geckos_csv['tail underwater depth'])
    geckos_csv['bodytail_drag'] = geckos_csv['body_drag'] + geckos_csv['tail_drag']   
    
        
       
    # find total vertical and horizontal forces generated per full cycle
    geckos_csv['tot_vert_mN'] = geckos_csv['fl.slap.force']+                \
        geckos_csv['bl.slap.force'] +geckos_csv['bl.stroke.ver.force'] +      \
        geckos_csv['fl.stroke.ver.force'] + geckos_csv['buoyancy_mN'] +     \
        geckos_csv['surface_tension_mN']
    geckos_csv['leg_strokehor_mN'] = geckos_csv['bl.stroke.hor.force']+           \
        geckos_csv['fl.stroke.hor.force']
    geckos_csv['tot_hor_mN'] = geckos_csv['bl.stroke.hor.force']+           \
        geckos_csv['fl.stroke.hor.force'] + (geckos_csv['bodytail_thrust'])

    # Total drag with feet, body, and tail
    geckos_csv['total_drag_w_feet'] = geckos_csv['bodytail_drag'] +            \
        geckos_csv['bl_drag'] + geckos_csv['fl_drag']
    #
    # # Find contributions
    # geckos_csv['buoyancy/tot_vert_mN'] = geckos_csv['buoyancy_mN']/geckos_csv['tot_vert_mN']*100
    # geckos_csv['surface_tension/tot_vert_mN'] = geckos_csv['surface_tension_mN']/geckos_csv['tot_vert_mN']*100
    # geckos_csv['slap/tot_vert_mN'] = (geckos_csv['fl.slap.force']*2+geckos_csv['bl.slap.force']*2)/geckos_csv['tot_vert_mN']*100
    # geckos_csv['stroke_vert/tot_vert_mN'] = (geckos_csv['fl.stroke.ver.force']*2+geckos_csv['bl.stroke.ver.force']*2) / geckos_csv['tot_vert_mN']*100
    # geckos_csv['stroke_hor/tot_hor_mN'] = (geckos_csv['fl.stroke.hor.force']*2+geckos_csv['bl.stroke.hor.force']*2) / geckos_csv['tot_hor_mN']*100
    # geckos_csv['undulation/tot_hor_mN'] = geckos_csv['bodytail_thrust']/geckos_csv['tot_hor_mN']*100
    
    geckos_csv['thrust-drag_ratio'] = geckos_csv['tot_hor_mN']/geckos_csv['total_drag_w_feet']



    # Impulse of drag over a stride
    geckos_csv['fl_drag_imp'] = geckos_csv['fl_drag']*stride
    geckos_csv['bl_drag_imp'] = geckos_csv['bl_drag']*stride
    geckos_csv['body_drag_imp'] = geckos_csv['body_drag']*stride
    geckos_csv['tail_drag_imp'] = geckos_csv['tail_drag']*stride
    geckos_csv['bodytail_drag_imp'] =  geckos_csv['bodytail_drag'] * stride
    geckos_csv['tot_drag_w_feet_imp'] = geckos_csv['fl_drag_imp'] + geckos_csv['bl_drag_imp']+ \
        geckos_csv['body_drag_imp'] + geckos_csv['tail_drag_imp']
        
    # Impulses 
    geckos_csv['body_imp'] = geckos_csv['body_thrust']*stride
    geckos_csv['tail_imp'] = geckos_csv['tail_thrust']*stride
    geckos_csv['bodytail_imp'] =geckos_csv['bodytail_thrust']*stride
    
    geckos_csv['fl_bl_vert_imp'] = (geckos_csv['fl.slap.imp']) + (geckos_csv['fl.stroke.ver.imp']) +\
        (geckos_csv['bl.slap.imp']) + (geckos_csv['bl.stroke.ver.imp'])
    geckos_csv['leg_strokehor_imp'] = (geckos_csv['bl.stroke.hor.imp']+           \
        geckos_csv['fl.stroke.hor.imp'])
        
        
    geckos_csv['total_vert_imp'] = (geckos_csv['fl.slap.imp']) + (geckos_csv['fl.stroke.ver.imp']) +\
        (geckos_csv['bl.slap.imp']) + (geckos_csv['bl.stroke.ver.imp']) + geckos_csv['surface_tension_imp'] # + geckos_csv['buoyancy_imp']
    geckos_csv['total_hor_imp'] = geckos_csv['fl.stroke.hor.imp']+\
        geckos_csv['bl.stroke.hor.imp']+geckos_csv['body_imp']+geckos_csv['tail_imp']
    
    # Find contributions
    geckos_csv['buoyancy/tot_vert_imp'] = geckos_csv['buoyancy_imp']/geckos_csv['total_vert_imp']*100
    geckos_csv['surf_tension/tot_vert_imp'] = geckos_csv['surface_tension_imp']/geckos_csv['total_vert_imp']*100
    geckos_csv['slap/tot_vert_imp'] = ((geckos_csv['fl.slap.imp']) + (geckos_csv['bl.slap.imp']))/geckos_csv['total_vert_imp']*100
    geckos_csv['stroke_vert/tot_vert_imp'] = ((geckos_csv['fl.stroke.ver.imp']) + (geckos_csv['bl.stroke.ver.imp'])) / geckos_csv['total_vert_imp']*100
    geckos_csv['stroke_hor/tot_hor_imp'] = ((geckos_csv['fl.stroke.hor.imp']) + (geckos_csv['bl.stroke.hor.imp'])) / geckos_csv['total_hor_imp']*100
    geckos_csv['undulation/tot_hor_imp'] = (geckos_csv['body_imp']+geckos_csv['tail_imp'])/geckos_csv['total_hor_imp']*100
    geckos_csv['vert_imp/min_imp'] = geckos_csv['total_vert_imp']/geckos_csv['min imp']
    geckos_csv['slapstrokevert_imp/min_imp'] = ((geckos_csv['fl.slap.imp']) + (geckos_csv['fl.stroke.ver.imp']) +\
        (geckos_csv['bl.slap.imp']) + (geckos_csv['bl.stroke.ver.imp']))/geckos_csv['min imp']
    geckos_csv['slap/min_imp'] = ((geckos_csv['fl.slap.imp']) + (geckos_csv['bl.slap.imp']))/geckos_csv['min imp']*100
    geckos_csv['stroke/min_imp'] = ((geckos_csv['fl.stroke.ver.imp']) + (geckos_csv['bl.stroke.ver.imp']))/geckos_csv['min imp']*100
    geckos_csv['st/min_imp'] = geckos_csv['surface_tension_imp']/geckos_csv['min imp']*100
    
    
    # Bond and Weber
    geckos_csv['Bo'] = (999.97*9.81*geckos_csv['fl_mm']/1000)/(0.07197/(geckos_csv['ffoot_mm']/1000))
    geckos_csv['We'] = (999.97*(geckos_csv['vel']**2)*(geckos_csv['ffoot_mm']/1000))/0.07197
    

    # Misc calculations
    geckos_csv['freq bmid'] = (2*np.pi/geckos_csv['per bmid'])
    geckos_csv['freq tmid'] = (2*np.pi/geckos_csv['per tmid'])

    geckos_csv.to_csv("python_edited_geckos_all.csv", sep=',',index=False)
    
    
    # ------
    
    # Average the trials across individuals to find average per individual
    ind_gecko_avg = geckos_csv.groupby(['expt','geckoID']).mean().reset_index()
    
    ind_gecko_avg.to_csv("ind_gecko_avg.csv", sep=',',index=False)

    
    # Get overall average
    averages = ind_gecko_avg.groupby('expt').mean()
    averages = averages.rename(index={'soap': 'soap_avg', 'water': 'water_avg'})
    std = ind_gecko_avg.groupby('expt').std()
    std = std.rename(index={'soap': 'soap_std', 'water': 'water_std'})
    means = pd.concat([averages,std])
    means = means.reset_index().sort_values('expt')

    
    means.to_csv("means.csv", sep=',',index=False)
    
    # Label low height geckos to high height geckos for drag forces only in water group
    avg_headht_water = ind_gecko_avg[ind_gecko_avg['expt']=='water'].mean()['head.ht']
    ind_gecko_avg['head_group'] = np.where(ind_gecko_avg['head.ht']<avg_headht_water, 'low', 'high')
    # Get average drag forces for only the water group
    ind_gecko_avg_water =  ind_gecko_avg[ind_gecko_avg['expt']=='water']
    drag_hilo = ind_gecko_avg_water.groupby('head_group').mean()['bodytail_drag']
    print(drag_hilo)
    drag_hi = ind_gecko_avg_water[ind_gecko_avg_water['head_group']=='high']
    drag_lo = ind_gecko_avg_water[ind_gecko_avg_water['head_group']=='low']
    T_val, p_val = stats.ttest_ind(drag_hi['bodytail_drag'],drag_lo['bodytail_drag'])
    print('T test of head high vs head low: T=', T_val, 'p=', p_val)
    
    print('\n')
    
    print(ind_gecko_avg.groupby('expt').mean()['total_drag_w_feet'])
    water_drag = ind_gecko_avg[ind_gecko_avg['expt']=='water']
    soap_drag = ind_gecko_avg[ind_gecko_avg['expt']=='soap']
    T_val, p_val = stats.ttest_ind(water_drag['total_drag_w_feet'],soap_drag['total_drag_w_feet'])
    print('T test of soap total drag vs water total drag: T=', T_val, 'p=', p_val)
    
    stats.f_oneway(water_drag['total_drag_w_feet'],soap_drag['total_drag_w_feet'])
    
    print('\n')
    
    print(ind_gecko_avg.groupby('expt').mean()['tot_hor_mN'])
    water = ind_gecko_avg[ind_gecko_avg['expt']=='water']
    soap = ind_gecko_avg[ind_gecko_avg['expt']=='soap']
    T_val, p_val = stats.ttest_ind(water['tot_hor_mN'],soap['tot_hor_mN'])
    print('T test of soap total thrust vs water total thrust: T=', T_val, 'p=', p_val)
    
    water = ind_gecko_avg[ind_gecko_avg['expt']=='water']
    soap = ind_gecko_avg[ind_gecko_avg['expt']=='soap']
    stats.ttest_ind(water['per bmid'],soap['per bmid'])

    
    
    # # Rest of code is for graphing drag, thrust, drag/thurst ratio. Ignore
    # dragthurst = ind_gecko_avg.groupby('expt').mean()['bodytail_thrust'].reset_index()
    # dragthurst.columns=['expt','thrust']
    # dragthurst['thurst_std'] = (ind_gecko_avg.groupby('expt').std().reset_index())['bodytail_thrust']
    # dragthurst['drag'] = (ind_gecko_avg.groupby('expt').mean().reset_index())['bodytail_drag']
    # dragthurst['drag_std'] = (ind_gecko_avg.groupby('expt').std().reset_index())['bodytail_drag']
    # dragthurst['thurst/drag'] = dragthurst['thrust']/dragthurst['drag']
    
    fig = plt.figure(figsize=(7,6))
    ax = fig.gca()
    plt.style.use('ggplot')


    plt.bar([0,1], averages['tot_hor_mN'], yerr=std['tot_hor_mN'], color=['gold', 'royalblue'])

    plt.ylim(0, 200)
    plt.xticks([0,1], ['Surfactant', 'Water'], fontsize=25)
    plt.yticks(fontsize=25)
    plt.ylabel('Thurst (mN)', fontsize=25)
    ax.set_facecolor('white')
    plt.subplots_adjust(bottom=None, top=0.9, left=0.2)
    fig.savefig('thurst.png')
    plt.show()


    fig = plt.figure(figsize=(7,6))
    ax = fig.gca()
    plt.style.use('ggplot')

    # fig.suptitle('Ge', fontsize=20)

    ax.bar([0,1],averages['total_drag_w_feet'],yerr=std['total_drag_w_feet'], color=['gold', 'royalblue'])

    plt.ylim(0, 200)
    plt.xticks([0,1], ['Surfactant', 'Water'], fontsize=25)
    plt.yticks(fontsize=25)
    plt.ylabel('Drag (mN)', fontsize=25)
    ax.set_facecolor('white')
    plt.subplots_adjust(bottom=None, top=0.9, left=0.2)
    fig.savefig('drag.png')
    plt.show()


    fig = plt.figure(figsize=(7,6))
    ax = fig.gca()
    plt.style.use('ggplot')
    # fig.suptitle('Ge', fontsize=20)

    ax.bar([0,1],averages['thrust-drag_ratio'], yerr= std['thrust-drag_ratio'], color=['gold', 'royalblue'])

    plt.xticks([0,1], ['Surfactant', 'Water'], fontsize=25)
    plt.yticks(fontsize=25)
    plt.ylabel('Thurst/Drag Percentage', fontsize=25)
    ax.set_facecolor('white')
    plt.subplots_adjust(bottom=None, top=0.9, left=0.2)
    fig.savefig('thurstdrag_ratio.png')
    plt.show()

    plt.close('all')
    
    # box plots
    import seaborn as sns
    
    # thurst
    

    fig = plt.figure(figsize=(7,6))
    ax = fig.gca()
    
    ax = sns.boxplot(x="expt", y="tot_hor_mN", hue="expt", data=ind_gecko_avg, 
        linewidth=2.5)
    
    ax.set_xticklabels(['Surfactant', 'Water'], fontsize=25)
    plt.ylim(0, 270)
    plt.xlabel('', fontsize=25)
    plt.yticks(fontsize=25)
    plt.ylabel('Thurst (mN)', fontsize=25)
    plt.subplots_adjust(left=0.2)
    ax.legend_.remove()
    fig.savefig('thurst_box_inds.png')
    plt.show()
    
    
    fig = plt.figure(figsize=(7,6))
    ax = fig.gca()
    
    ax = sns.boxplot(x="expt", y='total_drag_w_feet', hue="expt", data=ind_gecko_avg, linewidth=2.5)
    
    ax.set_xticklabels(['Surfactant', 'Water'], fontsize=25)
    plt.ylim(0, 270)
    plt.xlabel('', fontsize=25)
    plt.yticks(fontsize=25)
    plt.ylabel('Drag (mN)', fontsize=25)
    plt.subplots_adjust(left=0.2)
    ax.legend_.remove()
    fig.savefig('drag_box_inds.png')
    plt.show()
    
    
    fig = plt.figure(figsize=(7,6))
    ax = fig.gca()
    
    ax = sns.boxplot(x="expt", y='thrust-drag_ratio', hue="expt", data=ind_gecko_avg, linewidth=2.5)
    
    ax.set_xticklabels(['Surfactant', 'Water'], fontsize=25)
    plt.xlabel('', fontsize=25)
    plt.yticks(fontsize=25)
    plt.ylabel('Thurst/Drag Percentage', fontsize=25)
    plt.subplots_adjust(left=0.2)
    ax.legend_.remove()
    fig.savefig('thurstdrag_ratio_box_inds.png')
    plt.show()
    
    plt.close('all')