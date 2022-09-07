rule format_fasta_kraken2:
	input:
		kraken2_select = config["rdir"] + "/" + config["db_name"] + "/select_accessions.txt",
		linked = config["rdir"] + "/" + config["db_name"] + "/genomes/done"
	output:
		library = config["rdir"] + "/" + config["db_name"] + "/tmp/library.fna"
	params:
		script = config["wdir"] + "/scripts/format_kraken2.sh",
		gendir = config["rdir"] + "/" + config["db_name"] + "/genomes"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	threads: config["parallel_threads"]
	shell:
		"""
		cut -f1,3 {input.kraken2_select} | parallel -j{threads} '{params.script} {{}} {params.gendir}' | sed -e '/^>/!s/[a-z]/x/g' >> {output.library}
		"""

rule prelim_map:
	input:
		fasta = config["rdir"] + "/" + config["db_name"] + "/tmp/library.fna"
	output:
		map = config["rdir"] + "/" + config["db_name"] + "/tmp/prelim_map.txt"
	params:
		libdir = config["rdir"] + "/" + config["db_name"] + "/tmp"
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	shell:
		"""
		LC_ALL=C grep '^>' {input.fasta} | sed 's/^>//' > {params.libdir}/tmp.accnos
		NSEQ=$(wc -l {params.libdir}/tmp.accnos | cut -d' ' -f1)
		printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp.accnos | paste - <(cut -d'|' -f3 {params.libdir}/tmp.accnos) > {output.map}
		rm {params.libdir}/tmp.accnos
		"""

if config["kingdoms"]:
	rule format_taxid_conterminator:
		input:
			names = config["rdir"] + "/DB_taxonomy/names.dmp"
		output:
			kstring = config["rdir"] + "/" + config["db_name"] + "/decontamination/conterminator_string.txt",
			xstring = config["rdir"] + "/" + config["db_name"] + "/decontamination/conterminator_blacklist.txt"
		params:
			script = config["wdir"] + "/scripts/get_kingdoms_conterminator.R",
			kingdoms = config["kingdoms"],
			taxdir = config["rdir"] + "/DB_taxonomy/"
		conda:
			config["wdir"] + "/envs/r.yaml"
		shell:
			"""
			# parse taxid string for conterminator kingdoms parameter
			{params.script} -t {params.taxdir} -k "{params.kingdoms}" -s "{params.taxdir}/accessionTaxa.sql" -o {output.kstring} -x {output.xstring}
			"""

	rule detect_contamination:
		input:
			fasta = config["rdir"] + "/" + config["db_name"] + "/tmp/library.fna",
			map = config["rdir"] + "/" + config["db_name"] + "/tmp/prelim_map.txt",
			names = config["rdir"] + "/DB_taxonomy/names.dmp",
			kstring = config["rdir"] + "/" + config["db_name"] + "/decontamination/conterminator_string.txt",
			xstring = config["rdir"] + "/" + config["db_name"] + "/decontamination/conterminator_blacklist.txt"
		output:
			cmap = config["rdir"] + "/" + config["db_name"] + "/decontamination/cmap.txt",
			contam = config["rdir"] + "/" + config["db_name"] + "/decontamination/db_conterm_prediction"
		params:
			tmpdir = config["rdir"] + "/" + config["db_name"] + "/decontamination/tmp",
			taxdir = config["rdir"] + "/DB_taxonomy/",
			prefix = config["rdir"] + "/" + config["db_name"] + "/decontamination/db",
			cmem = config["cmem"]
		threads: config["parallel_threads"]
		log:
			config["rdir"] + "/logs/run_conterminator.log"
		shell:
			"""
			# prepare fasta header mapping file for conterminator
			cut -f2,3 {input.map} > {output.cmap}
			# run conterminator
			KSTR=$(cat {input.kstring})
			XSTR=$(cat {input.xstring})
			conterminator dna {input.fasta} {output.cmap} {params.prefix} {params.tmpdir} --mask-lower-case 1 --ncbi-tax-dump {params.taxdir} --threads {threads} --split-memory-limit {params.cmem} --blacklist $XSTR --kingdoms $KSTR &>> {log}
			"""

	rule filter_contamination:
		input:
			contam = config["rdir"] + "/" + config["db_name"] + "/decontamination/db_conterm_prediction",
			fasta = config["rdir"] + "/" + config["db_name"] + "/tmp/library.fna",
			map = config["rdir"] + "/" + config["db_name"] + "/tmp/prelim_map.txt"
		output:
			contam_filt = config["rdir"] + "/" + config["db_name"] + "/decontamination/db_conterm_prediction_filt",
			id_contam = config["rdir"] + "/" + config["db_name"] + "/decontamination/contam_id.accnos",
			fasta_noncontam = config["rdir"] + "/" + config["db_name"] + "/library/selection/library.fna",
			fasta_contam = config["rdir"] + "/" + config["db_name"] + "/decontamination/library_contam.fna",
			map_noncontam = config["rdir"] + "/" + config["db_name"] + "/library/selection/prelim_map.txt"
		conda:
			config["wdir"] + "/envs/bbmap.yaml"
		log:
			config["rdir"] + "/logs/contam_filter.log"
		shell:
			"""
			awk -v FS="\\t" -v OFS="\\t" '$5 >= 0 && $6 >= 0' {input.contam} > {output.contam_filt}
			cut -f2 {output.contam_filt} | sort | uniq > {output.id_contam}
			filterbyname.sh in={input.fasta} out={output.fasta_contam} names={output.id_contam} include=t ow=t &>> {log}
			filterbyname.sh in={input.fasta} out={output.fasta_noncontam} names={output.id_contam} include=f ow=t &>> {log}
			grep -v -F -f {output.id_contam} {input.map} > {output.map_noncontam}
			"""

	rule remove_contamination:
		input:
			contam = config["rdir"] + "/" + config["db_name"] + "/decontamination/db_conterm_prediction_filt",
			fasta_contam = config["rdir"] + "/" + config["db_name"] + "/decontamination/library_contam.fna",
			map_noncontam = config["rdir"] + "/" + config["db_name"] + "/library/selection/prelim_map.txt",
			fasta_noncontam = config["rdir"] + "/" + config["db_name"] + "/library/selection/library.fna"
		output:
			cleaned_fasta = config["rdir"] + "/" + config["db_name"] + "/decontamination/cleaned.fna",
			cleaned_map = config["rdir"] + "/" + config["db_name"] + "/decontamination/cleaned_map.txt"
		params:
			script = config["wdir"] + "/scripts/remove_contam_contigs.R",
			contam_dir = config["rdir"] + "/" + config["db_name"] + "/decontamination"
		conda:
			config["wdir"] + "/envs/r.yaml"
		log:
			config["rdir"] + "/logs/contam_cleaning.log"
		shell:
			"""
			{params.script} -i {input.fasta_contam} -c {input.contam} -o {output.cleaned_fasta} &>> {log}
			LC_ALL=C grep '^>' {output.cleaned_fasta} | sed 's/^>//' > "{params.contam_dir}/tmp.accnos"
			NSEQ=$(wc -l "{params.contam_dir}/tmp.accnos" | cut -d' ' -f1)
			printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - "{params.contam_dir}/tmp.accnos" | paste - <(cut -d'|' -f3 "{params.contam_dir}/tmp.accnos") > {output.cleaned_map}
			rm "{params.contam_dir}/tmp.accnos"
			cat {output.cleaned_map} >> {input.map_noncontam}
			cat {output.cleaned_fasta} >> {input.fasta_noncontam}
			# to be implemented later: remove tmp folder in krakendb to save disk space
			"""

else:
	rule collect_db_files:
		input:
			fasta = config["rdir"] + "/" + config["db_name"] + "/tmp/library.fna",
			map = config["rdir"] + "/" + config["db_name"] + "/tmp/prelim_map.txt"
		output:
			fasta = config["rdir"] + "/" + config["db_name"] + "/library/selection/library.fna",
			map = config["rdir"] + "/" + config["db_name"] + "/library/selection/prelim_map.txt"
		shell:
			"""
			mv {input.fasta} {output.fasta}
			mv {input.map} {output.map}
			"""

rule copy_taxonomy:
	input:
		names = config["rdir"] + "/DB_taxonomy/names.dmp"
	output:
		names = config["rdir"] + "/" + config["db_name"] + "/taxonomy/names.dmp"
	params:
		taxdir = config["rdir"] + "/DB_taxonomy/",
		kraken_dir = config["rdir"] + "/" + config["db_name"],
		kraken_tax = config["rdir"] + "/" + config["db_name"] + "/taxonomy/"
	shell:
		"""
		mkdir -p {params.kraken_dir}
		cp -r {params.taxdir}/*.dmp {params.kraken_tax}
		"""

# to avoid ftp issue, recreate kraken2 code for adding UniVec files
# https://github.com/DerrickWood/kraken2/blob/561cc73fababe1dfd996e553e36ea1aff5642ef8/scripts/download_genomic_library.sh#L102-L117
if config["univec"]:
	rule add_univec:
		output:
			fasta = config["rdir"] + "/" + config["db_name"] + "/library/" + config["univec"] + "/library.fna",
			map = config["rdir"] + "/" + config["db_name"] + "/library/" + config["univec"] + "/prelim_map.txt"
		params:
			ncbi_server = config["ncbi_server"],
			uv_name = config["univec"],
			libdir = config["rdir"] + "/" + config["db_name"] + "/library/" + config["univec"]
		conda:
			config["wdir"] + "/envs/kraken2.yaml"
		log:
			config["rdir"] + "/logs/add_krakendb_univec.log"
		shell:
			"""
			wget -O "{params.libdir}/tmp.fna" "{params.ncbi_server}/pub/UniVec/{params.uv_name}"
			# choosing random artificial taxid (this taxid must not exist elsewhere in the database)
			sed -i 's/^>/>kraken:taxid|123456789|/' "{params.libdir}/tmp.fna"
			dustmasker -in "{params.libdir}/tmp.fna" -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' > {output.fasta}
			rm "{params.libdir}/tmp.fna"
			grep '^>' {output.fasta} | sed 's/^>//' > {params.libdir}/tmp.accnos
			NSEQ=$(wc -l {params.libdir}/tmp.accnos | cut -d' ' -f1)
			printf 'TAXID\\n%.0s' $(seq 1 $NSEQ) | paste - {params.libdir}/tmp.accnos | paste - <(cut -d'|' -f2 {params.libdir}/tmp.accnos) > {output.map}
			rm {params.libdir}/tmp.accnos
			"""

rule build_kraken2_db:
	input:
		univec_fasta = config["rdir"] + "/" + config["db_name"] + "/library/" + config["univec"] + "/library.fna" if config["univec"] else [],
		univec_map = config["rdir"] + "/" + config["db_name"] + "/library/" + config["univec"] + "/prelim_map.txt" if config["univec"] else [],
		lib_fasta = config["rdir"] + "/" + config["db_name"] + "/library/selection/library.fna",
		lib_map = config["rdir"] + "/" + config["db_name"] + "/library/selection/prelim_map.txt",
		names = config["rdir"] + "/" + config["db_name"] + "/taxonomy/names.dmp"
	output:
		hash = config["rdir"] + "/" + config["db_name"] + "/hash.k2d",
		opts = config["rdir"] + "/" + config["db_name"] + "/opts.k2d",
		map  = config["rdir"] + "/" + config["db_name"] + "/seqid2taxid.map",
		taxo = config["rdir"] + "/" + config["db_name"] + "/taxo.k2d"
	params:
		dbdir = config["rdir"] + "/" + config["db_name"],
		kmer_len = config["kmer_len"],
		min_len = config["minimizer_len"],
		min_spaces = config["minimizer_spaces"],
		max_dbsize = config["max_dbsize"]
	threads: config["krakenbuild_threads"]
	conda:
		config["wdir"] + "/envs/kraken2.yaml"
	log:
		config["rdir"] + "/logs/" + config["db_name"] + "_build_db.log"
	shell:
		"""
		echo "kmer: {params.kmer_len}" &>> {log}
		echo "minimizer length: {params.min_len}" &>> {log}
		echo "minimizer spaces: {params.min_spaces}" &>> {log}
		kraken2-build --build --threads {threads} --db {params.dbdir} --kmer-len {params.kmer_len} --minimizer-len {params.min_len} --minimizer-spaces {params.min_spaces} --max-db-size {params.max_dbsize} &>> {log}
	"""
