# DARK_CURRENT
This is a repository of macros to analyse the dark current data of SiPMs and ARAPUCAs

The macros DC_read_adc_50us_inv.C and DC_read_adc_filt_50us_inv.C have been tested to work with the data that is provided by the configuration files. It is unclear how well the format of these macros is compatible with the IR02_Software that we use for the rest of the analysis.

Goals:
- Format the macros to read adc files (currently oscy.)
- Make use of the functions and cuts available for IR02_Software
- Use the fit to analyse the data (implement it in the same macro)
