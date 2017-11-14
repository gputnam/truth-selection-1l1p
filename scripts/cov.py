'''Produce covariance matrices for selected events.

Usage:

  $ python cov.py file.root

The input is an analysis tree ROOT file, i.e. the output from
`selection.py`.

This script should be modified to use the desired combinations of
systematic parameters.

N.B. This is somewhat of an abuse of the `gallery_framework` framework,
as it just calls some C++ functions that run over an arbitary ROOT file
(not art ROOT input), for ease of use.

A. Mastbaum <mastbaum@uchicago.edu>, 2017/09
'''

import sys
import argparse
from glob import glob
from ROOT import galleryfmwk

# Groups of parameters, to be varied together
wg = {
    'genie': [
        'genie_AGKYpT_Genie',
        'genie_AGKYxF_Genie',
        'genie_DISAth_Genie',
        'genie_DISBth_Genie',
        'genie_DISCv1u_Genie',
        'genie_DISCv2u_Genie',
        'genie_FermiGasModelKf_Genie',
        'genie_FermiGasModelSf_Genie',
        'genie_FormZone_Genie',
        'genie_IntraNukeNabs_Genie',
        'genie_IntraNukeNcex_Genie',
        'genie_IntraNukeNel_Genie',
        'genie_IntraNukeNinel_Genie',
        'genie_IntraNukeNmfp_Genie',
        'genie_IntraNukeNpi_Genie',
        'genie_IntraNukePIabs_Genie',
        'genie_IntraNukePIcex_Genie',
        'genie_IntraNukePIel_Genie',
        'genie_IntraNukePIinel_Genie',
        'genie_IntraNukePImfp_Genie',
        'genie_IntraNukePIpi_Genie',
        'genie_NC_Genie',
        'genie_NonResRvbarp1pi_Genie',
        'genie_NonResRvbarp2pi_Genie',
        'genie_NonResRvp1pi_Genie',
        'genie_NonResRvp2pi_Genie',
        'genie_ResDecayEta_Genie',
        'genie_ResDecayGamma_Genie',
        'genie_ResDecayTheta_Genie',
        'genie_ccresAxial_Genie',
        'genie_ccresVector_Genie',
        'genie_cohMA_Genie',
        'genie_cohR0_Genie',
        'genie_ncelAxial_Genie',
        'genie_ncelEta_Genie',
        'genie_ncresAxial_Genie',
        'genie_ncresVector_Genie',
        'genie_qema_Genie',
        'genie_qevec_Genie',
    ],
    'flux': [
        'expskin_FluxUnisim',
        'horncurrent_FluxUnisim',
        'kminus_PrimaryHadronNormalization',
        'kplus_PrimaryHadronFeynmanScaling',
        'kzero_PrimaryHadronSanfordWang',
        'nucleoninexsec_FluxUnisim',
        'nucleonqexsec_FluxUnisim',
        'nucleontotxsec_FluxUnisim',
        'piminus_PrimaryHadronSWCentralSplineVariation',
        'pioninexsec_FluxUnisim',
        'pionqexsec_FluxUnisim',
        'piontotxsec_FluxUnisim',
        'piplus_PrimaryHadronSWCentralSplineVariation',
    ],
    'xsr': [
        'xsr_scc_Fv3_XSecRatio',
        'xsr_scc_Fa3_XSecRatio',
    ],
}

parser = argparse.ArgumentParser("Scale input files and produce covariances")
parser.add_argument("input_file")
parser.add_argument("-d", "--directory", default=".")
parser.add_argument("-s", "--signal", action="store_true")
args = parser.parse_args()

# A group with all parameters varied
wg['all'] = []
for v in wg.values():
    wg['all'] += v

# POT scale factors for the Monte Carlo statistics (run utils/pot.py
# to get the MC POT for the dataset.

"""
sfm = 5.0e19 / 2.81456e+19  # BNB inclusive (!CCnue) sample
if args.signal:
   sfe = 5.0e19 / 1.86828e+21  # Signal sample
else:
   sfe = 5.0e19 / 8.47443e+21  # CCnue sample
"""
# UPDATED VALUES:
sfm = 5.0e19 / 2.72261e+19  # BNB inclusive (!CCnue) sample
if args.signal:
   sfe = 5.0e19 / 1.668035e+21  # Signal sample
else:
   sfe = 5.0e19 / 7.86397e+21  # CCnue sample

# Produce output file with covariances for each group
for w, ws in wg.items():
    print w
    m = galleryfmwk.TSCovariance()
    m.SetInputFile(args.input_file)
    m.SetScaleFactorE(sfe)
    m.SetScaleFactorMu(sfm)
    for ww in ws:
        m.AddWeight(ww)
    m.SetOutputFile(args.directory + '/cov_%s.root' % w)
    m.init()
    m.analyze()
    del m

# Produce output file with covariances for each individual parameter
for w in wg['all']:
    print w
    m = galleryfmwk.TSCovariance()
    m.SetInputFile(sys.argv[1])
    m.SetScaleFactorE(sfe)
    m.SetScaleFactorMu(sfm)
    m.AddWeight(w)
    m.SetOutputFile(args.directory + '/cov_%s.root' % w)
    m.init()
    m.analyze()
    del m

