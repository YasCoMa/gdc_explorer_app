# GDC Exploratory Data Analysis Portal for Coverage & Demographic subgroup stratification

This is a web application that allows exploring basic metrics of the projects data provided by the GDC portal (hosted by National Institute of Health (NIH)) for further deep analysis.

Check it out at https://yascoma.github.io/gdc_explorer_app/ 

## Summary
In the first tab, the focus is assessing the size of projects cohort (number of cases) and the number of open access files in each of them. Once choosing a specific project, it shows details about the stratification of its cases according to the demographic variables gender, race and ethnicity. Finally, it shows the data coverage for these samples for each data category (biospecimen, clinical, somatic structural variant, etc) and specific experimental strategy used to acquire the files (WXS, WGS, ATAC-Seq, etc).

In the second tab, for all the projects belonging to the TCGA program, we processed the clinical xml files looking for main metrics of interest that are used for later analysis in the published articles. And by choosing one project and one of the three demographic variable (race, gender, ethnicity), it is possible to check four main analysis sections: 
1) The proportions of the disease pathological stage according to each subgroup (ex.: female, male); 
2) The most prescribed drugs in the treatments for each subgroup; 
3) The distribution of the following continuous features across the subgroups:
	- Age at diagnosis
	- Qty. followup events
	- Qty. drugs administration
	- Qty. radiation prescriptions
	- Duration of treatment by drugs
	- Duration of treatment by radiation
	- Dosis used in treatment by drugs
	- Dosis used in treatment by radiation
4) The traditional Kaplan-Meier plot that models the survival probability across time.

In the third tab, the user may explore the mutational landscape by subgroup stratification containing a brief overview of the mutation annotations such as structural impact, the consequence in the transcription, the pathogenicity classification and the most affected genes.

In the fourth tab, the single nucleotide variation was integrated with the CIVIC db latest release of approved evidence for relations between mutations and drugs, regarding predictive entries that support statements asserting sensitivity/resistance. The user may explore for a given project, the number of cases for each mutation found in the project that matched civic entries, and check exactly the cases information that have the specific mutation, and the prescribed drugs (when available in the clinical metadata) for the cases. This part of the app was designed to be a demonstration of a simple dashboard to aid molecular tumor boards.

## Backend data processing
The data consumed by the rendering functions in the second tab can be reproduced using the script found in ./processing_scripts/process_data.py. It stores in the folder./data_processed.
There are also separated handlers to process and reorganize data for specific types of analysis, such as the handle\_bridgePharmacogenomics, that integrates data from the June, 2026 release of civic db with the mutations found in a specific chosen TCGA project. 
This script requires few packages as dependency, and the recipe to create a conda environment can be found in environment.yml.