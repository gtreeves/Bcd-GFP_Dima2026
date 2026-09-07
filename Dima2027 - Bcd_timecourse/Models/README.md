\### fit\_MWC\_global

This function fits the MWC model to experimental genome-wide dose/response map. Parameter screening is performed by systematic exploration within ranges of the input parameters: n, L, kn, ko. The capacity (ltot) is determined analytically. The function saves the best-fit parameter set and the SSE for all tested parameter combinations. 



\### plot\_enhancerstate

This function uses the parameters determined from the MWC model to simulate the fraction of active nuclei assuming it is equal to the probability that the enhancer is in the active state. Plots the minimum and maximum values predicted across all parameter sets, along with the best-fit prediction. Three cases are considered in which the promoter can be competent for transcription: ‘anybound’, ‘atleastonebound’, and ‘fullybound’ corresponding to Eq. (2-4). 



\### fit\_kon0\_tau

This function uses the parameters determined from the MWC model to fit transcriptional dynamics such as fraction active nuclei and transcriptional onset time using the multistate promoter model and estimate the values of kon0, tau, and t0 as reported in Supplementary Figure 6. Three cases are considered in which the promoter can be competent for transcription: ‘anybound’, ‘atleastonebound’, and ‘fullybound’ corresponding to Eq. (2-4). 











