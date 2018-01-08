import os
import sys

from collections import OrderedDict
import itertools

import numpy as np
import matplotlib.pyplot as plt

import chi2
import read_eff

def main(background_dir, signal_dir, shower_dist, track_dist, effs=None):
    selections = selections_list(track_dist, shower_dist) 
    significances = []

    for i in range(len(shower_dist) * len(track_dist)):
        background_name =  background_dir + "cov_bkg_%i" % i
        signal_name = signal_dir + "cov_sig_makross_%i" % i
        enu_mx, enu_m, enu_ex, enu_e, enu_s, enu_ss, cov = \
            chi2.load(background_name + '/cov_all.root', signal_name + '/cov_all.root')
        eb = np.hstack((enu_m, enu_e))
        fcov = cov / np.einsum('i,j', eb, eb)
        em = enu_m 
        ee = enu_e
        es = enu_s
    
        if effs is None:
            significances.append(
                   chi2.significance(em, ee, es, fcov,
                       pot=66e19, syst=True, mc=False, eff_m=0.3, eff_e=0.3)
            )
        else:
            significances.append(
                   chi2.significance(em, ee, es, fcov,
                       pot=66e19, syst=True, mc=False, eff_m=0.3, eff_e=0.3,
                       truth_eff_e=effs[i][0], truth_eff_s=effs[i][1], truth_eff_m=effs[i][2])
            )


    data = np.zeros(( len(track_dist) , len(shower_dist) ))
    for i,sig in enumerate(significances):
        t_ind = i / len(track_dist)
        s_ind = i % len(shower_dist)
        data[s_ind][t_ind] = sig
    
    print data
    fig, ax = plt.subplots()
    
    heatmap = ax.pcolor(data) 
    ax.set_xticks([0.5,2.5,4.5])
    ax.set_xlabel("track energy distortion")
    ax.set_xticklabels([track_dist[0], track_dist[2], track_dist[4]])
    ax.set_yticks([0.5, 2.5, 4.5])
    ax.set_yticklabels([shower_dist[0], shower_dist[2], shower_dist[4]])
    ax.set_ylabel("shower energy distortion")
    plt.colorbar(heatmap)
    fig.suptitle("Significance v. Energy Distortion")
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

    effs = []
    scratch_dir = "/pnfs/uboone/scratch/users/gputnam/smeared/"
    inc_to_nue_factor = 7863970000000000196608. / 509366499999999983616.

    for i in range(len(shower_dist) * len(track_dist)): 
        t_ind = i / len(track_dist)
        s_ind = i % len(shower_dist) 
        es_eff = read_eff.get_eff(47, scratch_dir + "smeared_signal_makross_weighted_5trial/", shower_dist[0], track_dist[0])[0]

        nue_data = read_eff.main(44, scratch_dir + "smeared_nue_5trial/", shower_dist[0], track_dist[0])
        inc_data = read_eff.main(1016, scratch_dir + "smeared_inclusive_slim_try2/", shower_dist[0], track_dist[0])

        em_eff = float(inc_data[4] + inc_data[5]) /float(inc_data[3])
        ee_eff = (float(nue_data[1] + nue_data[2]) + float(inc_data[1]+inc_data[2]) * inc_to_nue_factor) \
                 / (float(nue_data[0]) + float(inc_data[0])*inc_to_nue_factor)
        eff = (ee_eff,es_eff,em_eff)
        effs.append((ee_eff,es_eff,em_eff))

    main(background_dir, signal_dir, shower_dist, track_dist, effs)


