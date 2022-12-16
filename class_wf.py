# Create a class object that get stations and waveforms
class eq_data:
    import hvsrpy
    import sigpropy
    

    def __init__(self, client, t, starttime, endtime, net, sta, loc, channel):
        self.client = client
        self.t = t
        self.starttime = starttime
        self.endtime = endtime
        self.net = net
        self.sta = sta
        self.loc = loc
        self.channel = channel 
        waves = self.get_waveforms()
        
    
    
    def site_f0(self):
        # Return ew, ns, vt components to waves var
        # Decompose them into time series objects
        ew = waves[0]
        ns = waves[1]
        vt = waves[2]


        hvsr = hvsrpy.Sensor3c(ew=ew, ns=ns, vt=vt, meta=None)

        response = hvsr.hv(windowlength=windowlength, bp_filter=bp_filer, taper_width = width,
        bandwidth = bandwidth, resampling=resampling, method=method, f_low = None, f_high = None,
        azimuth = None)

        return response.mean_f0_fre(distribution_f0)


    '''Functions to retrie waveforms from FDSN client service'''
    '''A list of channels is required to retrieve east, north, and vertical component'''

    def get_waveforms(self):
        # return class values into local function
        client_0 = self.client 
        t = self.t
        starttime = self.starttime 
        endtime=self.endtime 
        net=self.net
        sta=self.sta
        loc=self.loc
        channel=self.channel
        
        
        
        
        
        
        from obspy.clients.fdsn import Client
        
        
        net, sta, loc = str(net), str(sta), str(loc)
        
        client = Client(client_0)
        
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

            st_e.merge()
            st_n.merge()
            st_z.merge()

            ew = sigpropy.TimeSeries(st_e.traces[0].data, dt = st_e[0].stats.delta)
            ns = sigpropy.TimeSeries(st_n.traces[0].data, dt = st_n[0].stats.delta)
            vt = sigpropy.TimeSeries(st_z.traces[0].data, dt = st_z[0].stats.delta)
            
        except (RuntimeError, TypeError, NameError):
    
                
            
        return ew, ns, vt

















    ###########################
    '''Data for methods to interact with each other'''
    '''Parameters to create hvsr object'''

    windowlength = 300

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
