# spatial-filter-for-dense-EEG-analysis

spatial-filter-for-dense-EEG-analysis

Spatial filters are time-geometry dependent coefficients, applied to the dense EEG signals


The aim of an spatial filter is to increase the EEG spatial resolution for better diagnosis of brain-maps, and to make a virtual sensor for extraction of neural sources on each area of the brain, as is the case for the LCMV beamformer, shown in the paper A robust beamforming approach for early detection of readiness potential with application to brain-computer interface systems.

The functionon scalef.m creates the surface spatial filters in 3 polar directions, called Laplacian as the representative of radial neural sources and tangential electric field as the representative of tangential direction of neural sources. 
The output matrix, applied to the EEG, will distinctively increase the sensitivity of the EEG to 
the activity of radial and tangential neural sources.


Also, the code contains an optimization function on the lcmv beamformer for acquiring the optimum dipole source direction from any brain location,

Then, the source is extracted by applying the lcmv beamformer to the surface spatial filtered EEG in the optum direction with maximum source power.


