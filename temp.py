from pyspecdata import *
d = find_file('A174', exp_type='ODNP_NMR_comp/field_dependent')
d = find_file('A174', exp_type='ODNP_NMR_comp/field_dependent', expnum='field_2')
d = find_file('A174', exp_type='ODNP_NMR_comp/field_dependent', expno='field_2')
fl.image(d)
image(d)
figure()
d.ft('t2',shift=True)
d.ft('ph1', unitary=True)
fl.image(d)
image(d)
d.mean('nscans')
image(d)
%hist -f temp.py
