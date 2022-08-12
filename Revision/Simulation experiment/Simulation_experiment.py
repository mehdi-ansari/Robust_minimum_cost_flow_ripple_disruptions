# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 13:11:48 2022

@author: mehdi
"""
import numpy as np
import pandas as pd
import ast

def Euclidean_distance(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

#Attack:
simulated_Attacks = pd.read_csv("1000attacks_data.csv") 

    
#Cost:
'''costfile = pd.read_csv("Cost_of_arcs.csv", header=None)
costfile = costfile.iloc[[5,25,45,65]]
costfile = costfile.reset_index(drop=True)'''

#Analysis:
resultFile = pd.read_csv("input.csv")
num_activeArcs = len(resultFile['activeArcs_record'].dropna())
Q = resultFile['radius'].dropna()
Q = Q[len(Q)-1]

df_dictionary = pd.DataFrame(columns=["arc",'C0','Yoptimal','location'])
for i in range(num_activeArcs):
    arc = resultFile['activeArcs_record'][i]
    df_dictionary = df_dictionary.append({
        'arc': int(arc),
        'C0': resultFile['C0a'][i],
        'Yoptimal': resultFile['Yoptimal_record'][i],
        'location': (resultFile['x'][i],resultFile['y'][i])}, ignore_index=True)

for i in range(len(simulated_Attacks)):
    scenario = simulated_Attacks.iloc[i]
    ds = ast.literal_eval(scenario['ripple_damage'])
    qs = [-0.001]
    for ripple in range(len(ds)):
        qs.append(round(Q*(ripple+1)/len(ds),2))    

    cost_and_flow = []
    for arc_index in range(len(df_dictionary)):
        Ca = df_dictionary['C0'][arc_index]
        flow = df_dictionary['Yoptimal'][arc_index]
        for attack in range(scenario['num_attacks']):
            epicenter = ast.literal_eval(scenario['attack_location'])[attack]
            for ripple in range(len(ds)):
                if(qs[ripple] < Euclidean_distance(epicenter, df_dictionary['location'][arc_index]) <= qs[ripple+1]):
                    Ca += ds[ripple]
                    break
        cost_and_flow.append(Ca*flow)
    df_dictionary[i] = cost_and_flow

total_costs = []
for i in range(len(simulated_Attacks)):
    total_costs.append(sum(df_dictionary[i]))
    
print(round(np.mean(total_costs),2))
print(round(np.std(total_costs),2))
print(round(np.min(total_costs),2))
print(round(np.max(total_costs),2))


'''activearcs_collection = []
for j in range(4):
    s = list(resultFile.iloc[j].dropna())
    activeArcs = s[14:]
    
    activearcs_collection.append(activeArcs)'''

    
'''
simulated_Attacks = pd.DataFrame(columns=['num_attacks', 'num_ripples', 'attack_location'
    numberofAttacks = np.random.randint(2,6), 'ripple_damage'])
for i in range(1000):
    numberofRipples = np.random.randint(5,10)
    DamageParameter = np.random.randint(20,61)
    simulated_Attacks = simulated_Attacks.append({'num_attacks': numberofAttacks,
                              'num_ripples': numberofRipples,
                              'attack_location': [(round(np.random.uniform(70),2),round(np.random.uniform(50),2)) for i in range(numberofAttacks)],
                              'ripple_damage': [round(DamageParameter/(1.5)**(k+1),2) for k in range(numberofRipples)]}, ignore_index=True)
    
simulated_Attacks.to_csv("1000attacks_data.csv", index=False)
'''