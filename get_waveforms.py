# Create a class object that get stations and waveforms

import hvsrpy
import sigpropy 
from obspy.clients.fdsn import Client 


def site_f0(ew, ns,vt, windowlength):

        ###########################
    '''Data for methods to interact with each other'''
    '''Parameters to create hvsr object'''

    # windowlength = 300

    # Boolean to control whether Butterworth filter is applied. 
    # Geopsy does not apply a bandpass filter.
    filter_bool = False        
    # Low-cut frequency for bandpass filter.
    filter_flow = 0.1                  
    # High-cut frequency for bandpass filter.
    filter_fhigh = 30                   
    # Filter order.
    filter_order = 5

    # Width of cosine taper {0. - 1.}. Geopsy default of 0.05 is equal to 0.1 -> 0.1 is recommended
    width = 0.1


    """FREQUENCY DOMAIN SETTING"""

    # Konno and Ohmachi smoothing constant. 40 is recommended.
    bandwidth = 40

    # Minimum frequency after resampling
    resample_fmin = 0.1 
    # Maximum frequency after resampling
    resample_fmax = 50
    # Number of frequencies after resampling
    resample_fnum = 200
    # Type of resampling {'log', 'linear'}
    resample_type = 'log'

    # Upper and lower frequency limits to restrict peak selection. To use the entire range use `None`.
    peak_f_lower = None
    peak_f_upper = None

    """HVSR Settings"""
    bp_filter = {"flag":filter_bool, "flow":filter_flow, "fhigh":filter_fhigh, "order":filter_order}
    resampling = {"minf":resample_fmin, "maxf":resample_fmax, "nf":resample_fnum, "res_type":resample_type}

    # Method for combining horizontal components {"squared-average", "geometric-mean", "single-azimuth"}.
    # Geopsy's default is "squared-average" -> "geometric-mean" is recommended.
    method = "geometric-mean"
    # If method="single-azimuth", set azimuth in degree clock-wise from north. If method!="single-azimuth", value is ignored.
    azimuth = 0

    # Boolean to control whether frequency domain rejection proposed by Cox et al. (2020) is applied.
    # Geopsy does not offer this functionality.
    rejection_bool = True
    # Number of standard deviations to consider during rejection. Smaller values will reject more windows -> 2 is recommended.
    n = 2
    # Maximum number of iterations to perform for rejection -> 50 is recommended.
    max_iterations = 50

    # Distribution of f0 {"lognormal", "normal"}. Geopsy default "normal" -> "lognormal" is recommended.
    distribution_f0 = "lognormal"
    # Distribution of mean curve {"lognormal", "normal"}. Geopsy default "lognormal" -> "lognormal" is recommended.
    distribution_mc = "lognormal"



    hvsr = hvsrpy.Sensor3c(ew=ew, ns=ns, vt=vt, meta=None)

    response = hvsr.hv(windowlength=windowlength, bp_filter=bp_filter, taper_width = width,
    bandwidth = bandwidth, resampling=resampling, method=method, f_low = None, f_high = None,
    azimuth = None)

    return response.mean_f0_frq(distribution_f0), response.mean_curve(distribution='lognormal'), response.frq


'''Functions to retrie waveforms from FDSN client service'''
'''A list of channels is required to retrieve east, north, and vertical component'''

def get_waves(client, starttime, endtime, net, sta, loc, channel):

    # import sigpropy
    # from obspy.clients.fdsn import Client

    net, sta, loc = str(net), str(sta), str(loc)
    
    client = Client(client)
    
    if isinstance(channel, list):
        pass
    else:
        print('channel has to be a list for E, N, and Z components')
        
        
        
    try:
        st_e = client.get_waveforms(network = net, station = sta, location = loc, channel = channel[0], starttime=starttime,
                            endtime=endtime, attach_response=True)

        st_n = client.get_waveforms(network = net, station = sta, location = loc, channel = channel[1], starttime=starttime,
                            endtime=endtime, attach_response=True)

        st_z = client.get_waveforms(network = net, station = sta, location = loc, channel = channel[2], starttime=starttime,
                            endtime=endtime, attach_response=True)

        st_e.remove_response(output = 'ACC')
        st_n.remove_response(output = 'ACC')
        st_z.remove_response(output = 'ACC')

        st_e.merge()
        st_n.merge()
        st_z.merge()

        ew = sigpropy.TimeSeries(st_e.traces[0].data, dt = st_e[0].stats.delta)
        ns = sigpropy.TimeSeries(st_n.traces[0].data, dt = st_n[0].stats.delta)
        vt = sigpropy.TimeSeries(st_z.traces[0].data, dt = st_z[0].stats.delta)
        
    except (RuntimeError, TypeError, NameError):
        pass
            
        
    return ew, ns, vt


