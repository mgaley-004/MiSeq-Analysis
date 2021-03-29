## HTS Analysis Repository

A home for R, mothur, and some python based solutions for analysing your HTS datasets.

### Tutorials/SOPs

Follow these standard workflows for a jumping off point.

- [Aligning and Classifying Sequences using mothur](https://mgaley-004.github.io/MiSeq-Analysis/Tutorials/AligningSequences.html) (formerly FastQ to Mothur)
  - Start here. Learn how to use the command line program `mothur` to align and identify DNA sequences from several regions of 16S.
- [Data Organization Throughout Your Project](https://mgaley-004.github.io/MiSeq-Analysis/Tutorials/DataOrganization.html)
  - Advice on setting up your metadata
  - Moving files to and from MSI
  - MSI workspace organization
  - R workspace organization
- [Read Scaling, Normalization, Transformation, Rarefactions and Diversity](https://mgaley-004.github.io/MiSeq-Analysis/Tutorials/rarefactions.html)
  - How to import your mothur output files into R for analysis
  - basics of phyloseq
  - OTU Richness, evenness, and coverage
  - The difference between proportional and relative abundance
  - Creating an even depth dataset
- [Ordinations and Redundancy Analysis](https://mgaley-004.github.io/MiSeq-Analysis/Tutorials/ordinationRDA.html)
  - Formerly known as beta diversity analysis
  - The differences between distance metrics
    - Bray Curtis, Jaccard, Manhattan, etc.
  - Mapping categorical and scalar values onto your microbial community
    - Principle components of analysis (PCA)
    - Principle coordinates of analysis (PCoA)
    - Non Metric Dimensional Scaling (NMDS)
  - Redundancy Analysis to connect environmental variables to your data
  - Minimum Spanning Trees for visualizing small differences
  - Cluster identification
- [Describing Major Lineages and Differential Abundance](https://mgaley-004.github.io/MiSeq-Analysis/Tutorials/MajorTaxonomicGroups.html)
  - Who is in  your samples?
  - Summarizing taxonomic differences at varying levels
  - Making stacked bar charts that aren't lies
  - Visualizations with heat maps
  - Statistically measuring differences in abundance
  - Using ALDEx2 to detect significant differences
- [FAPROTAX and Functions](https://mgaley-004.github.io/MiSeq-Analysis/Tutorials/faprotax_demo.html)
  - What are the microbes doing?
  - Assign functions
  - Visualization of FAPROTAX results
  - Make your own function tables
- [Source Sink Dynamics](https://mgaley-004.github.io/MiSeq-Analysis/Tutorials/SourceSink.html)
  - Set up Sourcetracker with phyloseq data
  - Visualize sourcetracker results
  - Case study

### General Help

- [Installing Python](https://mgaley-004.github.io/MiSeq-Analysis/Help/Python.html)
- [SLURM computing](https://mgaley-004.github.io/MiSeq-Analysis/Tutorials/SLURM/SlurmTransition.html) [(PDF)](https://mgaley-004.github.io/MiSeq-Analysis/Help/meeting121820.pdf)
  - Download metadata for this tutorial [here](https://github.com/mgaley-004/MiSeq-Analysis/tree/main/Help).
  - This tutorial is slightly outdated. MSI has retired login nodes for both `sftp` and `ssh` access. Now connect directly to msi using `ssh mesabi.msi.umn.edu -l yourx500` or `ssh mangi.msi.umn.edu -l yourx500` and only login once!
 
### Other Helpful Packages and Tutorials

- Having trouble in mothur? Check the [SOP](https://mothur.org/wiki/miseq_sop/) and the [wiki](https://mothur.org/wiki/mothur_manual/)
- Compositional analysis is complicated. The ALDEx2 package, used here frequently, is really helpful and the Gloor lab has extensive documentation and a helpful [workshop](https://github.com/ggloor/CoDa_microbiome_tutorial/wiki)
- I'm a big fan of [mixOmics](http://mixomics.org/). For fancier PCA options, check it out.


### Helpful Publications

- Miranda's long list of favorite bioinformatics articles will someday become available.

### Authors
Miranda Galey [@mgaley-004](https://github.com/mgaley-004/)

Contact at galey004-at-umn.edu

[Learn more about my work](https://miranda-galey.owlstown.net/)
