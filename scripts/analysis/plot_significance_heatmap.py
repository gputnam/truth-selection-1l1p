import os
import sys

from collections import OrderedDict
import itertools

import numpy as np
import matplotlib.pyplot as plt

import chi2

def main(background_dir, signal_dir, shower_dist, track_dist):
    selections = selections_list(track_dist, shower_dist) 
    significances = []

    for i in range(len(shower_dist) * len(track_dist)):
        background_name =  background_dir + "cov_bkg_%i" % i
        signal_name = signal_dir + "cov_sig_%i" % i
        enu_mx, enu_m, enu_ex, enu_e, enu_s, enu_ss, cov = \
            chi2.load(background_name + '/cov_all.root', signal_name + '/cov_all.root')
        eb = np.hstack((enu_m, enu_e))
        fcov = cov / np.einsum('i,j', eb, eb)
        em = enu_m
        ee = enu_e
        es = enu_s
    
        significances.append(
               chi2.significance(em, ee, es, fcov,
                   pot=66e19, syst=True, mc=False, eff_m=0.3, eff_e=0.3)
        )

    data = np.zeros(( len(track_dist) , len(shower_dist) ))
    for i,sig in enumerate(significances):
        t_ind = i / len(track_dist)
        s_ind = i % len(shower_dist)
        data[s_ind][t_ind] = sig
        print "NEXT"
        print t_ind, s_ind, sig
        print i, selections[i], track_dist[t_ind], shower_dist[s_ind]
    
    fig, ax = plt.subplots()
    
    heatmap = ax.pcolor(data) 
    ax.set_xticks([0.5,2.5,4.5])
    ax.set_xlabel("track energy distortion")
    ax.set_xticklabels([track_dist[0], track_dist[2], track_dist[4]])
    ax.set_yticks([0.5, 2.5, 4.5])
    ax.set_yticklabels([shower_dist[0], shower_dist[2], shower_dist[4]])
    ax.set_ylabel("shower energy distortion")
    plt.colorbar(heatmap)
    plt.show()

def selections_list(track, shower):
   s_len = len(shower)
   t_len = len(track)
   track_full = list(itertools.chain.from_iterable(itertools.repeat(x, s_len) for x in track))
   shower_full = shower * t_len
   return zip(track_full, shower_full)

if __name__ == "__main__":
    shower_dist = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    track_dist = [0,0.01,0.02,0.03,0.04,0.05]
    background_dir = sys.argv[1]
    signal_dir = sys.argv[2]
    main(background_dir, signal_dir, shower_dist, track_dist)


