import optparse
import ConfigParser
from collections import OrderedDict

def runcfg():
    p = optparse.OptionParser()
    p.add_option('--eflag','-e',default=False,action='store_true',help='Creat eventxxx.cfg for event information')

def set_parse_event(eventname,obsnames):
    filename = 'event-%s.cfg'%eventname
    try:
        temp = open(filename,'rb')
        print 'There is an existing config file for this event!'
        command = raw_input('continue? (y/n) ')
        if command == 'n':
            return
    except IOError:
        print 'This is a new event'
    print 'Generating config file for this event'
    configfile = open(filename,'wb')
    p = ConfigParser.RawConfigParser(dict_type=OrderedDict)
    strwd = 13
### Configuration parameter for event information ###
    p.add_section('This is a configure file for sfitpy package')
    p.add_section('Event Info')
    if eventname[:2] == 'ob':
        team = 'OGLE'
        num = eventname[-4:]
    elif eventname[:2] == 'mb':
        team = 'MOA'
        num = eventname[-3:]
    year = '20'+eventname[2:4]
    formalname = team+'-'+year+'-BLG-'+num
    p.set('Event Info','name'.ljust(strwd),formalname)
    p.set('Event Info','ra_j2000'.ljust(strwd),'')
    p.set('Event Info','dec_j2000'.ljust(strwd),'')
    p.add_section('Data Files')
    p.set('Data Files','observatories'.ljust(strwd),obsnames)
    p.set('Data Files','bands'.ljust(strwd),'')
    obsnames = obsnames.split()
    nobs = len(obsnames)
    for i in range(nobs):
        p.set('Data Files',obsnames[i].ljust(strwd),'')
    p.add_section('Source Info')
    p.set('Source Info','gamma_i','0.44')
    p.set('Source Info','gamma_v','0.72')
    p.set('Source Info','gamma_h','0.26')
    p.set('Source Info','gamma_l','0.')
### Parameters for the fitting ###
    p.add_section('Reference Parameters')
    strwd = 7
    p.set('Reference Parameters','tbegin'.ljust(strwd),'6000.')
    p.set('Reference Parameters','talert'.ljust(strwd),'8000.')
    p.set('Reference Parameters','t0par'.ljust(strwd),'7200.')
    p.set('Reference Parameters','t0base'.ljust(strwd),'')
    p.set('Reference Parameters','mjd'.ljust(strwd),'True')
    p.add_section('Fitting Parameters')
    strwd = 5
    p.set('Fitting Parameters','t0'.ljust(strwd),'')
    p.set('Fitting Parameters','u0'.ljust(strwd),'')
    p.set('Fitting Parameters','te'.ljust(strwd),'')
    p.set('Fitting Parameters','rho'.ljust(strwd),'')
    p.set('Fitting Parameters','pien'.ljust(strwd),'')
    p.set('Fitting Parameters','piee'.ljust(strwd),'')
    p.add_section('Flux Parameters')
    strwd = 15
    if 'spitzer' in obsnames:
        p.set('Flux Parameters','I_minus_L'.ljust(strwd),'')
        p.set('Flux Parameters','use_color'.ljust(strwd),'False')
    for i in range(nobs):
        p.set('Flux Parameters',('(fs,fb)_'+obsnames[i]).ljust(strwd),'')
    p.add_section('Error Rescalings')
    for i in range(nobs):
        p.set('Error Rescalings',('errfac_'+obsnames[i]).ljust(strwd),'1.0')
    p.write(configfile)
    configfile.close()
    return

if __name__ == '__main__':
    eventname = raw_input('Please give the event name (eg, ob140124): ')
    obsnames = raw_input('Please give the observatory names (eg, ogle spitzer): ')
#    eventname = 'ob140124'
#    obsnames = 'ogle spitzer'
    set_parse_event(eventname,obsnames)
#    set_parse_sfit(eventname)
