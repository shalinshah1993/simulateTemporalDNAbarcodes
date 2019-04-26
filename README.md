# (simulate) Temporal DNA barcodes
These software scripts are a part of the temporal DNA barcoding framework. It is a single-molecule imaging techqniue that uses time-domain to encode information for optical multiplexing. Although, the scripts were developed to our framework, they can be used for simulating any smFISH type where kinetics of single-molecule are required to be simulated.

If you find our scripts helpful, please cite our paper(s):

[1] `Shalin Shah, Abhishek Dubey, and John Reif. "Improved optical multiplexing with temporal DNA barcodes". ACS Synthetic Biology (2019) DOI: 10.1021/acssynbio.9b00010` [[PDF]](https://pubs.acs.org/doi/10.1021/acssynbio.9b00010)

[2] `Shalin Shah, Abhishek Dubey, and John Reif. "Programming temporal DNA barcodes for single-molecule fingerprinting". Nano Letters (2019) DOI: 10.1021/acs.nanolett.9b00590` [[PDF]](https://pubs.acs.org/doi/10.1021/acs.nanolett.9b00590)

[3] `Shalin Shah, and John Reif. "Temporal DNA Barcodes: A Time-Based Approach for Single-Molecule Imaging." International Conference on DNA Computing and Molecular Programming. Springer, Cham, 2018. DOI: 10.1007/978-3-030-00030-1_5` [[PDF]](https://link.springer.com/content/pdf/10.1007%2F978-3-030-00030-1_5.pdf)

## CTMC simulation scripts

A set of MATLAB scripts to program and simulate each parameter presented in the references [1, 2] above. Briefly, we used MATLAB’s SimBiology toolbox to represent our Markov models as chemical reactions and simulated them using SSA algorithm. For a better estimation of the experimental data, we added Gaussian noise to the simulated state chains since the combined effect of shot noise, dark noise, and all other detector noises are usually approximated by a Gaussian distribution. The unbinding rate constants for 7, 8, 9, and 10 nt were 10 ms, 60 ms, 550 ms, and 9 s, respectively. The binding rate constants for all the simulation were 1E6 /M– /s. More details on rate constants can be obtained from the papers. Random samples of on-time, off-time, double-blink time, and other higher-level on-times were calculated from the simulated state-chain and were binned into discrete bins and counts per bin. The observed distribution is then fitted with exponential probability density functions to estimate the mean parameters. The simulation experiments repeat several times to generate the parameter scatter plots.

- `t1_ssdna_devices_params_est.m`: tune the length of DNA to tune its kinetics and temporal barcode
- `t2_dbdomain_devices_params_est.m`: tune the number of reporter domains to tune the temporal barcode
- `t3_domainhide_devices_params_est.m`: tune the length of a multi-domain DNA hairpin to tune its kinetics and temporal barcode
- `t4_hairpin_devices_params_est.m`: tune the length of hairpin stem to change the dark time of temporal barcode

### Dependencies
- SimBiology - https://www.mathworks.com/products/simbiology.html
- Generate maximally perceptually-distinct colors - https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
