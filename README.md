# hvsr_analysis
Fundamental frequency of a site using eHVSR and mHVSR


# Background
### This project goal is to conduct Horizontal to Vertical Spectral Ratio (HVSR) analysis proposed by Nogoshi & Igarashi (1971) and Nakamura (1989) to estimate fundemental frequency <(f_{0})> of a site condition. There have been existing literature to validate the authenticity of f_{0} to its physics-based accuracy. We hope that we could develope a relationship between basin amplification in seismic ground motion models using f_{0} derived from HVSR

* The project utilizes HVSRPy package developed by John Vantessel. Please refer to their documentation for more information [HVSRPy](https://github.com/jpvantassel/hvsrpy)

* To eradicate the effect of altitude in seismic recordings, we decide to choose 10 stations which is denoted as **ground** from [Strong Motion Database](https://www.strongmotioncenter.org/)
and 5 earthquakes with magnitudes of 4 or greater

* We are under development to establish a framework for our analysis. Our first step is to gather data using `obspy` and analyze data using `hvsrpy`. Our data is generated at 2 stages: before instrument response and after instrument reponse


# HVSR Mean Curves for 5 earthquakes, 10 stations. Time duration includes 2-hr before and 2-hr after earthquakes (Categorized as microtremor HVSR or mHVSR)

