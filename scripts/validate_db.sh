# database validation

### simulate metaG reads from genomes not in the database
# select from assembly_taxonomy files accnos which are not in derep_assembly_taxonomy
# inlcude human as well (to make sure that this contamination will be recognized)
# select 2 random genomes per phylum
cd /vol/cloud/christiane/NCBI_taxdb_integration/Testing4
mkdir Simulation

PART="fungi invertebrate plant protozoa vertebrate_mammalian vertebrate_other viral"
cut -f1 tax_combined/gtdb_derep_taxonomy.txt | grep -v -F -f - gtdb/metadata/gtdb_taxonomy.txt > Simulation/non_db_genomes_taxonomy.txt
for part in $(echo "$PART")
do
  cut -f1 tax_combined/${part}"_derep_taxonomy.txt" | grep -v -F -f - ${part}/assembly_taxonomy.txt >> Simulation/non_db_genomes_taxonomy.txt
done

cd Simulation
sort -t$'\t' -k2,2 non_db_genomes_taxonomy.txt | sort -t$'\t' -u -k2,2 --merge > non_db_genomes_taxonomy_nr.txt
rm sim_genomes_taxonomy.txt
cut -f2 non_db_genomes_taxonomy_nr.txt | cut -d';' -f-2 | sort -t$'\t' | uniq > phylum_list.txt
while read line
do
  grep "${line}" non_db_genomes_taxonomy_nr.txt | shuf -n 2 >> sim_genomes_taxonomy.txt
done < phylum_list.txt
mv sim_genomes_taxonomy.txt tmp
sort -t$'\t' -k2,2 tmp | uniq > sim_genomes_taxonomy.txt # shuf may duplicate lines...
rm tmp

HUMAN=$(grep "s__Homo sapiens" sim_genomes_taxonomy.txt)
if [[ $HUMAN == "" ]]
then
  grep "s__Homo sapiens" non_db_genomes_taxonomy.txt | shuf -n 1 >> sim_genomes_taxonomy.txt
fi

cut -f2 ../gtdb/metadata/gtdb_download_info.txt | cat - ../*/assembly_url_genomic.txt | grep -F -f <(cut -f1 sim_genomes_taxonomy.txt) - > sim_genomes_download.txt
sed 's/.*\///' sim_genomes_download.txt | cut -d'_' -f2 | grep -v -F -f - sim_genomes_taxonomy.txt | cut -f1 | cut -d'_' -f2 | sed 's/\.[0-9]\+//' | grep -F -f - <(cut -f2 ../gtdb/metadata/gtdb_download_info.txt) >> sim_genomes_download.txt
mkdir genomes
aria2c -i sim_genomes_download.txt -c -l links.log --dir genomes --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads=28
gunzip genomes/*.gz

# metaG simulator
# https://github.com/nick-youngblut/MGSIM
# install dependencies first, then run pip install command
# maybe better to use setup.py?

# format genome table
# Important: no duplicate taxa allowed in genome table
paste <(sed 's/.*\//genomes\//' sim_genomes_download.txt | sed 's/\.gz//') <(sed 's/.*\//genomes\//' sim_genomes_download.txt | cut -d'_' -f2 | sed 's/\.[0-9]\+//') | sort -t$'\t' -k2,2 > tmp1
cut -f1 sim_genomes_taxonomy.txt | cut -d'_' -f2 | sed 's/\.[0-9]\+//' | paste sim_genomes_taxonomy.txt - | sort -t$'\t' -k3,3 > tmp2
diff <(cut -f2 tmp1) <(cut -f3 tmp2)
echo -e "Accnos\tTaxon\tFasta" > sim_genomes_table.txt
cut -f2,1 tmp2 | paste - <(cut -f1 tmp1) >> sim_genomes_table.txt
rm tmp1 tmp2

# run community simulation
MGSIM communities sim_genomes_table.txt sim_out
MGSIM reads --sr-seq-depth 1e7 --art-paired -n 16 sim_genomes_table.txt sim_out_abund.txt sim_out_reads


### create a database subset with low resolution for euks and high resolution for microbes
cat ../tax_combined/gtdb_derep_taxonomy.txt ../tax_combined/viral_derep_taxonomy.txt ../tax_combined/fungi_derep_taxonomy.txt ../tax_combined/protozoa_derep_taxonomy.txt | cut -f1 > db_subset_micro.accnos
grep -v "p__Streptophyta" ../tax_combined/plant_derep_taxonomy.txt | cut -f1 >> db_subset_micro.accnos
grep -F -f db_subset_micro.accnos ../tax_combined/derep_taxonomy_combined.txt > db_subset_micro_taxonomy.txt

grep "p__Streptophyta" ../tax_combined/plant_derep_taxonomy.txt | cut -f1 > db_subset_macro.accnos
cut -f1 ../tax_combined/invertebrate_derep_taxonomy.txt >> db_subset_macro.accnos
cut -f1 ../tax_combined/vertebrate_mammalian_derep_taxonomy.txt >> db_subset_macro.accnos
cut -f1 ../tax_combined/vertebrate_other_derep_taxonomy.txt >> db_subset_macro.accnos
grep -F -f db_subset_macro.accnos ../tax_combined/derep_taxonomy_combined.txt > db_subset_macro_taxonomy.txt

# estimate sizes
for i in $(seq 2 7)
do 
  cut -f1 db_subset_macro_taxonomy.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5,9 | sed -e 's/\.\.\/kraken2_genomes\/genome_files\///' -e 's/\.fa//' | sort -k2,2 | paste - <(sort -k1,1 db_subset_macro_taxonomy.txt | cut -f2 | cut -d';' -f-$i) | sort -t$'\t' -k3,3 -k1,1gr | sort -t$'\t' -u -k3,3 --merge | cut -f1 | paste -sd+ | bc -l
done > tmp
cut -f1 db_subset_macro_taxonomy.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5 | paste -sd+ | bc -l >> tmp
paste ../Summary_stats/tmp_taxlevels.txt tmp > db_subset_macro_sizes.txt
rm tmp
cut -f1 db_subset_micro_taxonomy.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5 | paste -sd+ | bc -l
# 321.5 GB micro + 173.5 GB macro (class level): too large without memory limitation
# only one representative per species for non-viral microbes
grep -v "d__Virus" db_subset_micro_taxonomy.txt > db_subset_micro_taxonomy_nv.txt
cut -f1 db_subset_micro_taxonomy_nv.txt | sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5,9 | sed -e 's/\.\.\/kraken2_genomes\/genome_files\///' -e 's/\.fa//' | sort -k2,2 | paste - <(sort -k1,1 db_subset_micro_taxonomy_nv.txt | cut -f2) | sort -t$'\t' -k3,3 -k1,1gr | sort -t$'\t' -u -k3,3 --merge | cut -f1 | paste -sd+ | bc -l
grep "d__Virus" db_subset_micro_taxonomy.txt > db_subset_micro_taxonomy_vi.txt
cut -f1 db_subset_micro_taxonomy_vi.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5 | paste -sd+ | bc -l
# 232.6 GB micro: that should work with class-level macro

# select largest genome per class for macro
cut -f1 db_subset_macro_taxonomy.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5,9 | sed -e 's/\.\.\/kraken2_genomes\/genome_files\///' -e 's/\.fa//' | sort -k2,2 | paste - <(sort -k1,1 db_subset_macro_taxonomy.txt | cut -f2 | cut -d';' -f-3) | sort -t$'\t' -k3,3 -k1,1gr | sort -t$'\t' -u -k3,3 --merge | cut -f2 | grep -F -f - db_subset_macro_taxonomy.txt > db_subset_class_macro_taxonomy.txt

# select largest genome per species for non-viral micro
cut -f1 db_subset_micro_taxonomy_nv.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5,9 | sed -e 's/\.\.\/kraken2_genomes\/genome_files\///' -e 's/\.fa//' | sort -k2,2 | paste - <(sort -k1,1 db_subset_micro_taxonomy_nv.txt | cut -f2) | sort -t$'\t' -k3,3 -k1,1gr | sort -t$'\t' -u -k3,3 --merge | cut -f2 | grep -F -f - db_subset_micro_taxonomy_nv.txt > db_subset_species_micro_taxonomy_nv.txt

# combine final tax tables
cat db_subset_class_macro_taxonomy.txt db_subset_species_micro_taxonomy_nv.txt db_subset_micro_taxonomy_vi.txt > db_subset_v1_taxonomy.txt
cut -f1 db_subset_v1_taxonomy.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | xargs ls -l | sed 's/\s\+/\t/g' | cut -f5 | paste -sd+ | bc -l

# prepare db directory (already has UniVec_Core)
cp -r ../kraken2_db/ ./kraken2_db1

# mask low complexity and build prelim_map
mkdir kraken2_db1/library/class/ # IMPORTANT: don't call directory 'added'
cut -f1 db_subset_v1_taxonomy.txt |  sed -e 's/^/\.\.\/kraken2_genomes\/genome_files\//' -e 's/$/\.fa/' | parallel -j28 'dustmasker -in {} -outfmt fasta' | sed -e '/^>/!s/[a-z]/x/g' >> kraken2_db1/library/class/library.fna
LC_ALL=C grep '^>' kraken2_db1/library/class/library.fna | sed 's/^>//' > kraken2_db1/library/class/tmp.accnos
NSEQ=$(wc -l kraken2_db1/library/added/tmp.accnos | cut -d' ' -f1)
printf 'TAXID\n%.0s' $(seq 1 $NSEQ) | paste - kraken2_db1/library/class/tmp.accnos | paste - <(cut -d'|' -f3 kraken2_db1/library/class/tmp.accnos) > kraken2_db1/library/class/prelim_map.txt
rm kraken2_db1/library/class/tmp.accnos

# build database
# use screen
scontrol -o show nodes
salloc --cpus-per-task=28 --nodelist="bibigrid-worker1-3-jdzmrjtottq9ctr" bash
conda activate source-tracking
srun /usr/bin/time -v kraken2-build --build --threads 26 --db kraken2_db1 --kmer-len 31 --minimizer-len 25 --minimizer-spaces 5 --max-db-size 460000000000
#Creating sequence ID to taxonomy ID map (step 1)...
#Sequence ID to taxonomy ID map complete. [50.440s]
#Estimating required capacity (step 2)...
#Estimated hash table requirement: 281012681872 bytes
#Capacity estimation complete. [3h19m7.301s]
#Building database files (step 3)...
#Taxonomy parsed and converted.
#CHT created with 17 bits reserved for taxid.
#Completed processing of 75247056 sequences, 432105004922 bp
#Writing data to disk...  complete.
#Database files completed. [9h27m13.416s]
#Database construction complete. [Total: 12h47m11.299s]
#        Command being timed: "kraken2-build --build --threads 26 --db kraken2_db1 --kmer-len 31 --minimizer-len 25 --minimizer-spaces 5 --max-db-size 460000000000"
#        User time (seconds): 757244.79
#        System time (seconds): 2446.90
#        Percent of CPU this job got: 1650%
#        Elapsed (wall clock) time (h:mm:ss or m:ss): 12:47:11
#        Average shared text size (kbytes): 0
#        Average unshared data size (kbytes): 0
#        Average stack size (kbytes): 0
#        Average total size (kbytes): 0
#        Maximum resident set size (kbytes): 349795148
#        Average resident set size (kbytes): 0
#        Major (requiring I/O) page faults: 54
#        Minor (reclaiming a frame) page faults: 150870106
#        Voluntary context switches: 57123202
#        Involuntary context switches: 29074168
#        Swaps: 0
#        File system inputs: 1758394488
#        File system outputs: 566472936
#        Socket messages sent: 0
#        Socket messages received: 0
#        Signals delivered: 0
#        Page size (bytes): 4096
#        Exit status: 0

# run classification
/usr/bin/time -v kraken2 --db kraken2_db1 --threads 10 --report sim_v1.kreport --output sim_v1.kraken --paired sim_out_reads/1/R1.fq sim_out_reads/1/R2.fq
#Loading database information... done.
#9674824 sequences (2902.45 Mbp) processed in 132.503s (4381.0 Kseq/m, 1314.29 Mbp/m).
#  9674310 sequences classified (99.99%)
#  514 sequences unclassified (0.01%)
#        Command being timed: "kraken2 --db kraken2_db1 --threads 10 --report sim_v1.kreport --output sim_v1.kraken --paired sim_out_reads/1/R1.fq sim_out_reads/1/R2.fq"
#        User time (seconds): 497.27
#        System time (seconds): 372.01
#        Percent of CPU this job got: 29%
#        Elapsed (wall clock) time (h:mm:ss or m:ss): 49:18.13
#        Average shared text size (kbytes): 0
#        Average unshared data size (kbytes): 0
#        Average stack size (kbytes): 0
#        Average total size (kbytes): 0
#        Maximum resident set size (kbytes): 275129108
#        Average resident set size (kbytes): 0
#        Major (requiring I/O) page faults: 0
#        Minor (reclaiming a frame) page faults: 69200750
#        Voluntary context switches: 1221554
#        Involuntary context switches: 29415
#        Swaps: 0
#        File system inputs: 450498208
#        File system outputs: 6367728
#        Socket messages sent: 0
#        Socket messages received: 0
#        Signals delivered: 0
#        Page size (bytes): 4096
#        Exit status: 0



# with memory limitation (70% of input file size) 
