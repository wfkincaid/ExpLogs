"""
Nutation
========

A standard echo where the 90 time is varied so 
that we are able to see when the signal rotates through 90 to 
180 degrees.
"""

import pyspecdata as psd
import os
import SpinCore_pp
from SpinCore_pp import get_integer_sampling_intervals, save_data, prog_plen
from Instruments.XEPR_eth import xepr
from SpinCore_pp.ppg import run_spin_echo
from datetime import datetime
import numpy as np
from numpy import r_

my_exp_type = "ODNP_NMR_comp/nutation"
assert os.path.exists(psd.getDATADIR(exp_type=my_exp_type))
beta_range_s_sqrtW = np.linspace(0.5e-6, 100e-6, 20)
# {{{importing acquisition parameters
config_dict = SpinCore_pp.configuration("active.ini")
prog_p90_us = prog_plen(beta_range_s_sqrtW, config_dict)
(
    nPoints,
    config_dict["SW_kHz"],
    config_dict["acq_time_ms"],
) = get_integer_sampling_intervals(
    config_dict["SW_kHz"], config_dict["acq_time_ms"]
)
# }}}
# {{{add file saving parameters to config dict
config_dict["type"] = "nutation"
config_dict["date"] = datetime.now().strftime("%y%m%d")
config_dict["echo_counter"] += 1
# }}}
# {{{set phase cycling
ph1_cyc = r_[0, 2]
ph2_cyc = r_[0, 2]
nPhaseSteps = len(ph1_cyc) * len(ph2_cyc)
# }}}
# {{{let computer set field
input(
    "I'm assuming that you've tuned your probe to %f since that's what's in   "
    " your .ini file. Hit enter if this is true"
    % config_dict["carrierFreq_MHz"]
)
field_G = config_dict["carrierFreq_MHz"] / config_dict["gamma_eff_MHz_G"]
print(
    "Based on that, and the gamma_eff_MHz_G you have in your .ini file, I'm   "
    " setting the field to %f" % field_G
)
with xepr() as x:
    assert field_G < 3700, "are you crazy??? field is too high!"
    assert field_G > 3300, "are you crazy?? field is too low!"
    field_G = x.set_field(field_G)
    print("field set to ", field_G)
# }}}
# {{{check total points
total_pts = nPoints * nPhaseSteps
assert total_pts < 2**14, (
    "You are trying to acquire %d points (too many points) -- either change SW"
    "    or acq time so nPoints x nPhaseSteps is less than 16384\nyou could"
    " try    reducing the acq_time_ms to %f"
    % (total_pts, config_dict["acq_time_ms"] * 16384 / total_pts)
)
# }}}
data = None
for idx, beta in enumerate(beta_range_s_sqrtW):
    # Just loop over the 90 times and set the indirect axis at the end
    # just like how we perform and save IR data
    data = run_spin_echo(
        deadtime_us=config_dict["deadtime_us"],
        deblank_us=config_dict["deblank_us"],
        nScans=config_dict["nScans"],
        indirect_idx=idx,
        indirect_len=len(prog_p90_us),
        ph1_cyc=ph1_cyc,
        ph2_cyc=ph2_cyc,
        amplitude=config_dict["amplitude"],
        adcOffset=config_dict["adc_offset"],
        carrierFreq_MHz=config_dict["carrierFreq_MHz"],
        nPoints=nPoints,
        nEchoes=config_dict["nEchoes"],
        plen=beta,
        repetition_us=config_dict["repetition_us"],
        tau_us=config_dict["tau_us"],
        SW_kHz=config_dict["SW_kHz"],
        ret_data=data,
    )
data.rename("indirect", "beta")
data.setaxis("beta", beta_range_s_sqrtW).set_units("beta", "s√W")
data.set_prop("prog_p90_us", prog_p90_us)
# {{{ chunk and save data
data.chunk("t", ["ph2", "ph1", "t2"], [2, 2, -1])
data.setaxis("ph1", ph1_cyc / 4).setaxis("ph2", ph2_cyc / 4)
if config_dict["nScans"] > 1:
    data.setaxis("nScans", r_[0 : config_dict["nScans"]])
data.reorder(["nScans", "ph2", "ph1", "beta", "t2"])
data.set_units("t2", "s")
data.set_prop("postproc_type", "spincore_nutation_v6")
data.set_prop("coherence_pathway", {"ph1": +1, "ph2": -2})
data.set_prop("acq_params", config_dict.asdict())
config_dict = save_data(data, my_exp_type, config_dict, "echo")
config_dict.write()
