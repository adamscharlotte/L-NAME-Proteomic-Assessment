Proteomic assessment of C57BL/6 hippocampi after non-selective pharmacological inhibition of nitric oxide synthase activity
========

In a previous study we found that long-term inhibition of all NOS isoforms led to impaired visuospatial learning and memory in C57BL/6 mice <https://doi.org/10.3390/biomedicines9121905>. As a follow-up we investigated the progressive molecular signatures in hippocampal tissue. Non-selective nitric oxide synthase inhibition was induced by pharmacological treatment of male C57BL/6 mice from the age of 8 weeks onward with 0.5 mg/mL N(G)-nitro-L-arginine methyl ester (L-NAME) in drinking water for 2 weeks, 8 weeks, and 16 weeks. At sacrifice, the hippocampus of one cerebral hemisphere were collected and immediately frozen in liquid nitrogen. Tissues were homogenized, proteins were actracted, samples were digested with trypsin, and peptides were purified and injected for LC-MS/MS analysis on an Ultimate 3000 RSLCnano system (Thermo) in line connected to a Q Exactive HF Biopharma mass spectrometer (Thermo).

Spectral identification was performed with MaxQuant (version 1.6.11.0) using the Andromeda search engine with default search settings including a false discovery rate set at 1% on PSM and peptide level. The spectra were searched against the mouse proteins in the Swiss-Prot Reference Proteome database (database release version of June 2019 containing 22,282 mouse protein sequences). The mass tolerance for precursor and fragment ions were set to 4.5 and 20 ppm, respectively. Enzyme specificity was set as C-terminal to arginine and lysine, also allowing cleavage at proline bonds with a maximum of two missed cleavages. Variable modifications were set to oxidation of methionine residues and acetylation of protein N-termini, and carbamidomethylation of cysteine residues was set as fixed modification.

The peptides.txt output was analyzed further with MSqRob2 <https://github.com/statOmics/msqrob2>. We performed data normalisation, protein inference, protein quantification, and looked for differentially expressed proteins (DEPs) at 2 weeks, 8 weeks, and 16 weeks. Additionally we performed a 2-way ANOVA.

The specific settings used to compile the lists of DEPs can be found in the R Markdown scripts <https://github.com/adamscharlotte/L-NAME-Proteomic-Assessment/tree/main/MSqRob2/Rmd>.

Contact
-------

For more information you can send an email to <charlotte.adams@uantwerpen.be>.
