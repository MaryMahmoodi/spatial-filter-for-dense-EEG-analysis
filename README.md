# spatial-filter-for-dense-EEG-analysis

spatial-filter-for-dense-EEG-analysis

* Spatial filters are time-geometry dependent coefficients, applied to the dense EEG signals


The aim of an spatial filter is to increase the EEG spatial resolution for better diagnosis of brain-maps, and to make a virtual sensor for extraction of neural sources on each area of the brain, as is the case for the 
# LCMV beamformer,
shown in the paper A robust beamforming approach for early detection of readiness potential with application to brain-computer interface systems.

* The functionon scalef.m creates the surface spatial filters in 3 polar directions, called 
# Laplacian 
as the representative of "radial" neural sources and "tangential" electric field as the representative of tangential direction of neural sources. 
The output matrix, applied to the EEG, will distinctively increase the sensitivity of the EEG to 
the activity of radial and tangential neural sources.


* Also, the function optimummomentum_lcmv_sekihara2015.m applies an optimization to the lcmv beamformer in order to acquire the optimum dipole source direction from any brain location. written according to the book entitled "electromagnetic brain imaging", also "adaptive spatial filters for electromagnetic brain imaging" by  Sekihara et.al. 

* Then, the source signal is extracted by applying the lcmv beamformer to the surface spatial filtered EEG in the optimum direction with maximum source power.


