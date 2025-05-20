#!/usr/bin/env python
# coding: utf-8

'''
Code for simulating generic difference equation form.
Author: Srini Venkatramanan (srini@virginia.edu)
Last Modified: 20th May 2025
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse as ap
import re

parser = ap.ArgumentParser(description='Generic difference equation Simulator.')
parser.add_argument('--model', type=str,
                    help='Disease model specification')
parser.add_argument('--config', type=str,
                    help='Disease parameters')
parser.add_argument('--outfile', type=str,help='Output file')
parser.add_argument('--plotfile',type=str,help='Plot file')
args = parser.parse_args()


def custom_make_translation(text, translation):
    regex = re.compile('|'.join(map(re.escape, translation)))
    return regex.sub(lambda match: translation[match.group(0)], text)

# ### Read specifications
model = pd.read_csv(args.model,sep=':',comment="#",names=['key','value'],skip_blank_lines=True).set_index('key')
cfg = pd.read_csv(args.config,sep='=',comment="#",names=['key','value'],skip_blank_lines=True).set_index('key')

# ### Initialize StateArray
states = model.loc['States']['value'].split(';')
dim = len(states)
T = int(cfg.loc['Duration'].values[0])
X = np.zeros((dim,T))
state_dict = {}
c=0
for s in states:
    state_dict[s] = X[c]
    c+=1

# ### Set initial conditions
N = int(cfg.loc['N']['value'])
for s in states:
    try:
        state_dict[s][0] = int(cfg.loc[s]['value'])
    except:
        pass

# ### Read parameters
params = model.loc['Parameters']['value'].split(';')
param_dict = {}
for p in params:
    param_dict[p] = float(cfg.loc[p]['value'])
# ### Create model equations
equations = [x for x in model.index.values if x.startswith('T')]
trans_dict = {}

for s in states:
    trans_dict[s] = s
for e in equations:
    eqn,cat,par = model.loc[e]['value'].split(';')
    src = eqn.split('->')[0]
    dest = eqn.split('->')[1]
    
    par_term = ''
    if par in param_dict.keys():
        par_term = "param_dict['{}']".format(par)
    else:
        par_term = par
        map_dict = {}
        for k in param_dict.keys():
            if k in par_term:
                map_dict[k] = "param_dict['{}']".format(k)
        par_term = custom_make_translation(par_term,map_dict)
        #print(par_term)
            
        
    if cat=='': ## no catalyst  
        if src=='':
            trans_dict[dest] += " + {}".format(par_term)
        elif dest=='':
            trans_dict[src] += " - {}".format(par_term)
        else:
            trans_dict[src] += " - {}*{}".format(par_term,src)
            trans_dict[dest] += " + {}*{}".format(par_term,src)
    else:
        if src=='':
            trans_dict[dest] += " + {}*{}".format(par_term,cat)
        elif dest=='':
            trans_dict[src] += " - {}*{}".format(par_term,cat)
        else:
            trans_dict[src] += " - {}*{}*{}/N".format(par_term,src,cat)
            trans_dict[dest] += " + {}*{}*{}/N".format(par_term,src,cat)        
for s in states:
    for x in states:
        trans_dict[s] = trans_dict[s].replace(x,x+'[t]')
        
model_str = ''
for s in states:
    model_str+=(s+'[t+1] = '+trans_dict[s]+'\n')
model_str = model_str.strip()

print('**** Model Equations ****')
print(model_str.replace('param_dict','').replace("['",'').replace("']",''))
for s in states:
    model_str = model_str.replace(s,"state_dict['{}']".format(s))

# ### Create external update equations
externals = [x for x in model.index.values if x.startswith('X')]
ext_dict = {}
ext_str = ''
for ext in externals:
    eqn,fname = model.loc[ext].value.split(';')
    
    var = fname.split('.')[0]
    sched = pd.read_csv(fname,index_col=0,header=None,sep=' ',names=[var])
    ext_dict[ext] = sched.reindex(range(0,T)).fillna(0)[var].values
    
    src = eqn.split('->')[0]
    dest = eqn.split('->')[1]
    
    ext_str+= src+'[t] = '+src+'[t] - min('+src+"[t],ext_dict['{}'][t])".format(ext)+'\n'
    ext_str+= dest+'[t] = '+dest+'[t] + min('+src+"[t],ext_dict['{}'][t])".format(ext)+'\n'

print('**** External Updates ****')
print(ext_str.replace('ext_dict','').replace("['",'').replace("']",''))
for s in states:
    ext_str = ext_str.replace(s,"state_dict['{}']".format(s))

# ### Simulate the system and save output
for t in range(T-1):    
    exec(ext_str)
    exec(model_str)
out = pd.DataFrame(X.T,columns=states)
    
if args.outfile:
    out.to_csv(args.outfile)

if args.plotfile:
    plt.figure()
    out.plot()
    plt.savefig(args.plotfile,dpi=100)

