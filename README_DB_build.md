## Taxonomy database flavours

#### Snakefile_db_build + config_DB.yaml

The rule in "Snakefile_db_build" select the genomes to create the database according to the flavour specified in the configuration file "config_DB.yaml".

The possible flavours are:

flavour_main: "vanilla" or "hires"
  - vanilla: includes the GTDB rs202 species cluster representatives (original GTDB, not derepG)
  - hires: includes the GTDB rs202 genomes dereplicated with derepG

flavour_sec: "prok", "organelle_euk", "organelle_euk_virus"
  - vanilla_organelles_euk: includes vanilla + ncbi-organelles + coarse euk
  - vanilla_organelles_euk_virus: includes vanilla + ncbi-organelles + coarse euk + CheckV v0.6 (representative genomes)

  - hires_organelle_euk: includes GTDB rs202 - derepG + ncbi-organelles + mikroeuk + coarse-macro-euk
  - hires-organelle_euk_virus: includes GTDB rs202 - derepG + ncbi-organelles + mikroeuk + coarse-macro-euk + CheckV v0.6 (derepG)


If specified in the config_DB.yaml, the flavours include also the custom selection of prokaryotic, eukaryotic and viral genomes. The lists are in the folder assets/.

At the moment this module creates a folder with the flavour name (ex: "vanilla/") with the following structure:
vanilla/
  select_accessions.txt
  select_taxonomy.txt
  genomes/


#### Snakefile_kraken2 + config_kraken2.yaml

Rules and configuration to build a kraken2 database based on the flavour selection.

This module create the database inside the previously created flavour folder:

vanilla/
  select_accessions.txt
  select_taxonomy.txt
  genomes/
  taxonomy/
  db/
    kraken2/
    library/
    *.k2d files
