
# __Population genomic analyses reveal geographic structure in Rhizoctonia solani AG1-IA isolates associated with different crops in the USA and the Caribbean__
by Juanita Gil, Kensy Rodriguez-Herrera, Vanina Castroagudin, Felipe Dalla Lana, Xin-Gen Zhou, Sara Thomas-Sharma, Terry Spurlock, Jim Correll, Pierre Gladieux and Alejandro Rojas

GitHub repository for data analysis for publication.

## Data
Raw sequence reads are deposited on NCBI in the SRA under BioProject ID PRJNA1292868 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1292868). 

__Abstract__

Soilborne fungi threaten food security, but limited knowledge of their population biology, including genetic variability and population structure, hinders the development of effective strategies to prevent crop losses. _Rhizoctonia solani_ AG-1 (Basidiomycota: Agaricomycetes) is a significant soilborne pathogen worldwide, divided into subgroups based on host range and molecular diversity. _Rhizoctonia solani_ AG1-IA causes sheath blight in rice and aerial blight in soybean, two devastating diseases in these economically important crops. Rice and soybean are often used in rotation, increasing inoculum in the fields over cropping seasons. Monitoring genetic variability and structure in _R. solani_ AG1-IA is critical to understand the population biology of the pathogen and aid in disease management practices, including screening and developing new resistant cultivars. A total of 145 isolates of _R. solani_ collected between 1993 and 2022 from different hosts and states in the USA were sequenced. The population genetic structure was inferred with clustering approaches based on approximately two million biallelic single nucleotide polymorphisms. _Rhizoctonia solani_ AG1-IA showed relatively high genetic diversity and little departure from mutation-drift equilibrium, pointing to an ancient, widespread, indigenous pathogen in the USA. While populations had a clear geographic structure, they lacked host specialization, suggesting dispersal is mainly distance limited. Shared ancestry between populations and the discovery of clonal lineages, however, indicated recent connections between geographic areas. Our work provides valuable insights into the evolutionary history and population biology of _R. solani_ AG1-IA, offering a foundation for developing targeted management strategies and resistant crop varieties.

# Organization of this repository

# Scripts:
For most of the open source sotware used in this study, please refer to the corresponding manual for details on how to use. A few example scripts are provided, but these mostly correspond to custom scripts for data processing and analysis.

runMapping_bowtie2.sh 

## Variants_genotyping_filtering 
- detection of variants and genotyping: runVariantCallingBcftools.sh
- annotation and filtering: runAnnotateVCF.sh, filter_vcf_stdin.py

## Pop_structure
Assessment of population structure
- generate fasta file from population vcf for SplitsTree: vcf2fastadiploid.py
- calculate ancestry coefficients: ancestry.R
- generate plots: pophelper.R

## Clone_correction
- Assess clonality in the population and identify MLGs: poppr.R

## Selection
  
Linkage disequilibrium, recombination, and genome scan for selective sweeps
Over-representation analysis of genes in regions under positive selection


