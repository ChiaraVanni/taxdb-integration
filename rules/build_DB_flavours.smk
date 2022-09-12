## Build vanilla DB
# vanilla: GTDB (using the species representatives)
if config["flavour_main"] == "vanilla":
	rule collect_vanilla_prok:
		input:
			gtdb_reps = config["rdir"] + "/gtdb/metadata/gtdb_reps_tax.txt",
			gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
		output:
			pro_select = config["rdir"] + "/" + config["db_name"] + "/gtdb_select_accessions.txt"
		params:
			gendir = config["rdir"] + "/gtdb/reps_genomes",
			outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
		shell:
			"""
			mkdir -p {params.outdir}
			find {params.gendir} -type f -name '*.gz' | while read line
			do
			  ln -sf "$line" {params.outdir}
			done
			cut -f1 {input.gtdb_reps} | grep -F -f - {input.gen2taxid} > {output.pro_select}
			"""

	if config["custom_gtdb_post_derep"] != "n":
		rule add_custom_vanilla_prok:
			input:
				gtdb_reps = config["rdir"] + "/gtdb/metadata/gtdb_reps_tax.txt",
				gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt",
				tax_added = config["rdir"] + "/tax_combined/pro_custom_post_derep_taxonomy.txt"
			output:
				tax_select = config["rdir"] + "/" + config["db_name"] + "/custom_pro_select_taxonomy.txt",
				acc_select = config["rdir"] + "/" + config["db_name"] + "/custom_pro_select_accessions.txt"
			params:
				script = config["wdir"] + "/scripts/coarse_selection_custom.R",
				rank = "species",
				nmax = 1,
				add = config["custom_gtdb_post_derep"],
				gendir = config["rdir"] + "/derep_combined",
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			conda:
				config["wdir"] + "/envs/r.yaml"
			log:
				config["rdir"] + "/logs/select_coarse_custom_gtdb.log"
			shell:
				"""
				{params.script} -t {input.gtdb_reps} -c {params.add} -r {params.rank} -n {params.nmax} -o {output.tax_select} &>>{log}
				mkdir -p {params.outdir}
				find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {output.tax_select}) | while read line
				do
					ln -sf "$line" {params.outdir}
				done
				cut -f1 {output.tax_select} | grep -F -f - {input.gen2taxid} > {output.acc_select}
				"""
# If the secondary flavour is set to "prok", check that all GTDB were collected
	if config["flavour_sec"] == "prok":
		rule check_vanilla_prok:
			input:
				select_gtdb = config["rdir"] + "/" + config["db_name"] + "/gtdb_select_accessions.txt",
				custom_pro = config["rdir"] + "/" + config["db_name"] + "/custom_pro_select_accessions.txt" if config["custom_gtdb_post_derep"] != "n" else []
			output:
				tax = config["rdir"] + "/" + config["db_name"] + "/select_taxonomy.txt",
				select = config["rdir"] + "/" + config["db_name"] + "/select_accessions.txt",
				checked = config["rdir"] + "/" + config["db_name"] + "/genomes/done"
			params:
				gtdb_tax = config["rdir"] + "/gtdb/metadata/gtdb_reps_tax.txt",
				custom_pro_tax = config["rdir"] + "/" + config["db_name"] + "/custom_pro_select_taxonomy.txt" if config["custom_gtdb_post_derep"] != "n" else [],
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			shell:
				"""
				cat {input.select_gtdb} {input.custom_pro} > {output.select}
				cat {params.gtdb_tax} {params.custom_pro_tax} > {output.tax}

				if [[ $(cat {output.select} | wc -l) == $(find {params.outdir} -name '*.gz' | wc -l) ]]; then
					touch {output.checked}
				fi
				rm {input.select_gtdb} {input.custom_pro} {params.custom_pro_tax}
				"""

# If the secondary flavour is not set to "prok", get the genomes for the organelles
	if config["flavour_sec"] != "prok":
		rule collect_organelle:
			input:
				tax_organelle = config["rdir"] + "/tax_combined/organelle_derep_taxonomy.txt",
				gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
			output:
				organelle_select = config["rdir"] + "/" + config["db_name"] + "/organelle_select_accessions.txt",
				linked = config["rdir"] + "/" + config["db_name"] + "/genomes/organelle_done"
			params:
				gendir = config["rdir"] + "/derep_combined",
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			shell:
				"""
				cat {input.tax_organelle} | cut -f1 | grep -F -f - {input.gen2taxid} > {output.organelle_select}
				mkdir -p {params.outdir}
				find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {output.organelle_select}) | while read line
				do
					ln -sf "$line" {params.outdir}
				done
				touch {output.linked}
				"""
# Also add the eukaryotes in coarse-mode
		rule select_euk_coarse:
			input:
				tax = config["rdir"] + "/tax_combined/{library_name}_derep_taxonomy.txt",
				meta = config["rdir"] + "/{library_name}/metadata/genome_metadata.txt",
				gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
			output:
				tax_select = config["rdir"] + "/" + config["db_name"] + "/{library_name}_select_taxonomy.txt",
				euk_select = config["rdir"] + "/" + config["db_name"] + "/{library_name}_select_accessions.txt"
			params:
				script = config["wdir"] + "/scripts/coarse_ncbi_selection.R",
				rank = lambda wildcards: config["rank_coarse"][wildcards.library_name],
				sublineage = lambda wildcards: config["coarse_sublineage"][wildcards.library_name],
				rank_sublineage = lambda wildcards: config["rank_sublineage"][wildcards.library_name],
				nmax = config["nmax_coarse"]
			conda:
				config["wdir"] + "/envs/r.yaml"
			log:
				config["rdir"] + "/logs/select_coarse_{library_name}.log"
			shell:
				"""
				{params.script} -t {input.tax} -m {input.meta} -r {params.rank} -s "{params.sublineage}" -R {params.rank_sublineage} -n {params.nmax} -o {output.tax_select} &>>{log}
				cut -f1 {output.tax_select} | grep -F -f - {input.gen2taxid} > {output.euk_select}
				"""

		rule collect_euk_coarse:
			input:
				euk_select = config["rdir"] + "/" + config["db_name"] + "/{library_name}_select_accessions.txt"
			output:
				linked = config["rdir"] + "/" + config["db_name"] + "/genomes/{library_name}_done"
			params:
				gendir = config["rdir"] + "/derep_combined",
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			shell:
				"""
				mkdir -p {params.outdir}
				find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {input.euk_select}) | while read line
				do
					ln -sf "$line" {params.outdir}
				done
				touch {output.linked}
				"""

		if config["custom_ncbi_post_derep"] != "n":
			rule add_custom_ncbi_coarse:
				input:
					tax_select = expand(config["rdir"] + "/" + config["db_name"] + "/{library_name}_select_taxonomy.txt", library_name = LIBRARY_NAME),
					gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt",
					tax_added = config["rdir"] + "/tax_combined/euk_custom_post_derep_taxonomy.txt"
				output:
					custom_select = config["rdir"] + "/" + config["db_name"] + "/custom_euk_select_taxonomy.txt",
					euk_select = config["rdir"] + "/" + config["db_name"] + "/custom_euk_select_accessions.txt"
				params:
					script = config["wdir"] + "/scripts/coarse_selection_custom.R",
					rank = config["rank_custom"],
					nmax = config["nmax_coarse"],
					dbdir = config["rdir"] + "/" + config["db_name"],
					add = config["custom_ncbi_post_derep"],
					gendir = config["rdir"] + "/derep_combined",
					outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
				conda:
					config["wdir"] + "/envs/r.yaml"
				log:
					config["rdir"] + "/logs/select_coarse_custom_ncbi.log"
				shell:
					"""
					cat {input.tax_select} > "{params.dbdir}/tmp_ncbi_select_taxonomy.txt"
					{params.script} -t "{params.dbdir}/tmp_ncbi_select_taxonomy.txt" -c {params.add} -r {params.rank} -n {params.nmax} -o {output.custom_select} &>>{log}
					mkdir -p {params.outdir}
					find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {output.custom_select}) | while read line
					do
						ln -sf "$line" {params.outdir}
					done
					cut -f1 {output.custom_select} | grep -F -f - {input.gen2taxid} > {output.euk_select}
					rm "{params.dbdir}/tmp_ncbi_select_taxonomy.txt"
					"""

# If the secondary flavour is set to "organelle_euk", check the collection of GTDB, organelles and eukaryotes
	if config["flavour_sec"] == "organelle_euk":
		rule check_vanilla_organelle_euk:
			input:
				select_gtdb = config["rdir"] + "/" + config["db_name"] + "/gtdb_select_accessions.txt",
				select_euk = expand(config["rdir"] + "/" + config["db_name"] + "/{library_name}_select_accessions.txt", library_name = LIBRARY_NAME),
				euk_linked = expand(config["rdir"] + "/" + config["db_name"] + "/genomes/{library_name}_done", library_name = LIBRARY_NAME),
				select_organelle = config["rdir"] + "/" + config["db_name"] + "/organelle_select_accessions.txt"
			output:
				flav_dir = config["rdir"] + "/" + config["db_name"],
				tax = config["rdir"] + "/" + config["db_name"] + "/select_taxonomy.txt",
				select = config["rdir"] + "/" + config["db_name"] + "/select_accessions.txt",
				checked = config["rdir"] + "/" + config["db_name"] + "/genomes/done"
			params:
				custom_euk = config["rdir"] + "/" + config["db_name"] + "/custom_euk_select_accessions.txt" if config["custom_ncbi_post_derep"] != "n" else [],
				custom_pro = config["rdir"] + "/" + config["db_name"] + "/custom_pro_select_accessions.txt" if config["custom_gtdb_post_derep"] != "n" else [],
				gtdb_tax = config["rdir"] + "/gtdb/metadata/gtdb_reps_tax.txt",
				organelle_tax = config["rdir"] + "/tax_combined/organelle_derep_taxonomy.txt",
				coarse_euk_tax = config["rdir"] + "/" + config["db_name"] + "/{library_name}_select_taxonomy.txt",
				custom_pro_tax = config["rdir"] + "/" + config["db_name"] + "/custom_pro_select_taxonomy.txt" if config["custom_gtdb_post_derep"] != "n" else [],
				custom_euk_tax = config["rdir"] + "/" + config["db_name"] + "/custom_euk_select_taxonomy.txt" if config["custom_ncbi_post_derep"] != "n" else [],
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			shell:
				"""
				cat {input.select_gtdb} {input.select_euk} \
					{input.select_organelle} {params.custom_euk} {params.custom_pro} > {output.select}
				cat {params.gtdb_tax} {params.organelle_tax} {params.coarse_euk_tax} \
					{params.custom_pro_tax} {params.custom_euk_tax} > {output.tax}

				if [[ $(cat {output.select} | wc -l) == $(find {params.outdir} -name '*.gz' | wc -l) ]]; then
					touch {output.checked}
				fi
				rm {params.flav_dir}/*_select_accessions.txt {params.flav_dir}/*_select_taxonomy.txt
				"""

# If the secondary flavour is set to "organelle_euk_virus" (last option), get the viral genomes (checkv representatives)
	if config["flavour_sec"] == "organelle_euk_virus":
		rule collect_checkv_coarse:
			input:
				reps_tax = config["rdir"] + "/checkv/checkv_reps_taxonomy.txt",
				reps_fna = config["rdir"] + "/checkv/checkv_reps.fna",
				gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
			output:
				vir_select = config["rdir"] + "/" + config["db_name"] + "/checkv_select_accessions.txt",
				linked = config["rdir"] + "/" + config["db_name"] + "/genomes/checkv_done"
			params:
				gendir = config["rdir"] + "/checkv/reps_genomes",
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			conda:
				config["wdir"] + "/envs/parallel.yaml"
			threads: config["parallel_threads"]
			shell:
				"""
				mkdir -p {params.gendir}
				if [[ $(find {params.gendir} -type f -name '*.gz' | wc -l) != $(grep -c '^>' {input.reps_fna}) ]]; then
					cd {params.gendir}
					cat {input.reps_fna} | awk '{{ if (substr($0, 1, 1)==">") {{filename=(substr($0,2) ".fa")}} print $0 > filename }}'
					find . -type f -name '*.fa' | parallel -j {threads} gzip {{}}
				fi
				cut -f1 {input.reps_tax} | grep -F -f - {input.gen2taxid} > {output.vir_select}

				mkdir -p {params.outdir}
				find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {output.vir_select}) | while read line
				do
					ln -sf "$line" {params.outdir}
				done
				touch {output.linked}
				"""

		if config["custom_checkv_post_derep"] != "n":
			rule add_custom_checkv_coarse:
				input:
					reps_tax = config["rdir"] + "/checkv/checkv_reps_taxonomy.txt",
					gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt",
					tax_added = config["rdir"] + "/tax_combined/vir_custom_post_derep_taxonomy.txt"
				output:
					custom_select = config["rdir"] + "/" + config["db_name"] + "/custom_vir_select_taxonomy.txt",
					vir_select = config["rdir"] + "/" + config["db_name"] + "/custom_vir_select_accessions.txt"
				params:
					script = config["wdir"] + "/scripts/coarse_selection_custom.R",
					rank = "species",
					nmax = 1,
					add = config["custom_checkv_post_derep"],
					gendir = config["rdir"] + "/derep_combined",
					outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
				conda:
					config["wdir"] + "/envs/r.yaml"
				log:
					config["rdir"] + "/logs/select_coarse_custom_vir.log"
				shell:
					"""
					{params.script} -t {input.reps_tax} -c {params.add} -r {params.rank} -n {params.nmax} -o {output.custom_select} &>>{log}
					mkdir -p {params.outdir}
					find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {output.custom_select}) | while read line
					do
						ln -sf "$line" {params.outdir}
					done
					cut -f1 {output.custom_select} | grep -F -f - {input.gen2taxid} > {output.vir_select}
					"""

# check the collection of GTDB, organelles, eukaryotes and viruses
		rule check_vanilla_organelle_euk_virus:
			input:
				select_checkv = config["rdir"] + "/" + config["db_name"] + "/checkv_select_accessions.txt",
				select_gtdb = config["rdir"] + "/" + config["db_name"] + "/gtdb_select_accessions.txt",
				select_euk = expand(config["rdir"] + "/" + config["db_name"] + "/{library_name}_select_accessions.txt", library_name = LIBRARY_NAME),
				euk_linked = expand(config["rdir"] + "/" + config["db_name"] + "/genomes/{library_name}_done", library_name = LIBRARY_NAME),
				select_organelle = config["rdir"] + "/" + config["db_name"] + "/organelle_select_accessions.txt"				
			output:
				tax = config["rdir"] + "/" + config["db_name"] + "/select_taxonomy.txt",
				select = config["rdir"] + "/" + config["db_name"] + "/select_accessions.txt",
				checked = config["rdir"] + "/" + config["db_name"] + "/genomes/done"
			params:
				custom_euk = config["rdir"] + "/" + config["db_name"] + "/custom_euk_select_accessions.txt" if config["custom_ncbi_post_derep"] != "n" else [],
				custom_pro = config["rdir"] + "/" + config["db_name"] + "/custom_pro_select_accessions.txt" if config["custom_gtdb_post_derep"] != "n" else [],
				custom_vir = config["rdir"] + "/" + config["db_name"] + "/custom_vir_select_accessions.txt" if config["custom_checkv_post_derep"] != "n" else [],
				gtdb_tax = config["rdir"] + "/gtdb/metadata/gtdb_reps_tax.txt",
				organelle_tax = config["rdir"] + "/tax_combined/organelle_derep_taxonomy.txt",
				checkv_tax = config["rdir"] + "/checkv/checkv_reps_taxonomy.txt",
				coarse_euk_tax = config["rdir"] + "/" + config["db_name"] + "/{library_name}_select_taxonomy.txt",
				custom_pro_tax = config["rdir"] + "/" + config["db_name"] + "/custom_pro_select_taxonomy.txt" if config["custom_gtdb_post_derep"] != "n" else [],
				custom_euk_tax = config["rdir"] + "/" + config["db_name"] + "/custom_euk_select_taxonomy.txt" if config["custom_ncbi_post_derep"] != "n" else [],
				custom_vir_tax = config["rdir"] + "/" + config["db_name"] + "/custom_vir_select_taxonomy.txt" if config["custom_checkv_post_derep"] != "n" else [],
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			shell:
				"""
				cat {input.select_checkv} {input.select_gtdb} \
					{input.select_euk} {input.select_organelle} \
					{params.custom_euk} {params.custom_pro} {params.custom_vir} > {output.select}
				cat {params.gtdb_tax} {params.organelle_tax} {params.checkv_tax} \
					{params.coarse_euk_tax} {params.custom_pro_tax} \
					{params.custom_euk_tax} {params.custom_vir_tax} > {output.tax}

				if [[ $(cat {output.select} | wc -l) == $(find {params.outdir} -name '*.gz' | wc -l) ]]; then
			    	touch {output.checked}
				fi
				"""

##################################################################################################################################################
## Build the high resolution DB
# Basic hires is GTDB using the results/genomes from derepG
if config["flavour_main"] == "hires":
	rule select_genomes_highres_prok:
		input:
			gtdb = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
			added_nuc_gtdb = config["rdir"] + "/tax_combined/pro_custom_post_derep_taxonomy.txt" if config["custom_gtdb_post_derep"] != "n" else [],
			gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
		output:
			pro_select = config["rdir"] + "/" + config["db_name"] + "/pro_select_accessions.txt"
		shell:
			"""
			cat {input.gtdb} {input.added_nuc_gtdb} | cut -f1 | grep -F -f - {input.gen2taxid} > {output.pro_select}
			"""

	rule collect_genomes_highres_prok:
		input:
			pro_select = config["rdir"] + "/" + config["db_name"] + "/pro_select_accessions.txt"
		output:
			linked = config["rdir"] + "/" + config["db_name"] + "/genomes/done"
		params:
			gendir = config["rdir"] + "/derep_combined",
			outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
		shell:
			"""
			mkdir -p {params.outdir}
			find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {input.pro_select}) | while read line
			do
				ln -sf "$line" {params.outdir}
			done
			if [[ $(cat {input.pro_select} | wc -l) == $(find {params.outdir} -name '*.gz' | wc -l) ]]
			then
				touch {output.linked}
			fi
			"""

# If the secondary flavour is set to "prok", check that all GTDB were collected
	if config["flavour_sec"] == "prok":
		rule check_hires_prok:
			input:
				pro_select = config["rdir"] + "/" + config["db_name"] + "/pro_select_accessions.txt",
			output:
				tax = config["rdir"] + "/" + config["db_name"] + "/select_taxonomy.txt",
				select = config["rdir"] + "/" + config["db_name"] + "/select_accessions.txt",
				checked = config["rdir"] + "/" + config["db_name"] + "/genomes/done"
			params:
				flav_dir = config["rdir"] + "/" + config["db_name"],
				gtdb_tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
				custom_pro_tax = config["rdir"] + "/tax_combined/pro_custom_post_derep_taxonomy.txt" if config["custom_gtdb_post_derep"] != "n" else [],
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			shell:
				"""
				cat {input.pro_select} > {output.select}
				cat {params.gtdb_tax} {params.custom_pro_tax} > {output.tax}
				if [[ $(cat {output.select} | wc -l) == $(find {params.outdir} -name '*.gz' | wc -l) ]]; then
					touch {output.checked}
				fi

				rm {params.flav_dir}/*_select_accessions.txt {params.flav_dir}/*_select_taxonomy.txt
				"""

# If the secondary flavour is not set to "prok", get the accessions/genomes for the organelles (results from derepG)
	if config["flavour_sec"] != "prok":
		rule collect_organelle:
			input:
				tax_organelle = config["rdir"] + "/tax_combined/organelle_derep_taxonomy.txt",
				gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
			output:
				organelle_select = config["rdir"] + "/" + config["db_name"] + "/organelle_select_accessions.txt",
				linked = config["rdir"] + "/" + config["db_name"] + "/genomes/organelle_done"
			params:
				gendir = config["rdir"] + "/derep_combined",
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			shell:
				"""
				cat {input.tax_organelle} | cut -f1 | grep -F -f - {input.gen2taxid} > {output.organelle_select}
				mkdir -p {params.outdir}
				find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {output.organelle_select}) | while read line
				do
					ln -sf "$line" {params.outdir}
				done
				touch {output.linked}
				"""
# Also collect the micro-eukaryotes accessions/genomes (results from derepG)
		rule collect_euk_micro:
			input:
				tax = config["rdir"] + "/tax_combined/{library_micro}_derep_taxonomy.txt",
				gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
			output:
				mikro_select = config["rdir"] + "/" + config["db_name"] + "/{library_micro}_select_accessions.txt",
				linked = config["rdir"] + "/" + config["db_name"] + "/genomes/{library_micro}_done"
			params:
				gendir = config["rdir"] + "/derep_combined",
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			shell:
				"""
				cut -f1 {input.tax} | grep -F -f - {input.gen2taxid} > {output.mikro_select}
				mkdir -p {params.outdir}
				find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {output.mikro_select}) | while read line
				do
					ln -sf "$line" {params.outdir}
				done
				touch {output.linked}
				"""
# And collect the macro-eukaryotes accessions/genomes (coarse-mode)
		rule select_euk_macro:
			input:
				tax = config["rdir"] + "/tax_combined/{library_macro}_derep_taxonomy.txt",
				meta = config["rdir"] + "/{library_macro}/metadata/genome_metadata.txt",
				gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
			output:
				tax_select = config["rdir"] + "/" + config["db_name"] + "/{library_macro}_select_taxonomy.txt",
				macro_select = config["rdir"] + "/" + config["db_name"] + "/{library_macro}_select_accessions.txt"
			params:
				script = config["wdir"] + "/scripts/coarse_ncbi_selection.R",
				rank = lambda wildcards: config["rank_macro"][wildcards.library_macro],
				sublineage = lambda wildcards: config["macro_sublineage"][wildcards.library_macro],
				rank_sublineage = lambda wildcards: config["macro_rank_sublineage"][wildcards.library_macro],
				nmax = config["nmax_macro"]
			conda:
				config["wdir"] + "/envs/r.yaml"
			log:
				config["rdir"] + "/logs/select_microeuk_{library_macro}.log"
			shell:
				"""
				{params.script} -t {input.tax} -m {input.meta} -r {params.rank} -s "{params.sublineage}" -R {params.rank_sublineage} -n {params.nmax} -o {output.tax_select} &>>{log}
				cut -f1 {output.tax_select} | grep -F -f - {input.gen2taxid} > {output.macro_select}
				"""

		rule collect_euk_macro:
			input:
				macro_select = config["rdir"] + "/" + config["db_name"] + "/{library_macro}_select_accessions.txt"
			output:
				linked = config["rdir"] + "/" + config["db_name"] + "/genomes/{library_macro}_done"
			params:
				gendir = config["rdir"] + "/derep_combined",
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			shell:
				"""
				mkdir -p {params.outdir}
				find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {input.macro_select}) | while read line
				do
					ln -sf "$line" {params.outdir}
				done
				touch {output.linked}
				"""

		if config["custom_ncbi_post_derep"] != "n":
			rule add_custom_ncbi_euk:
				input:
					tax_select = expand(config["rdir"] + "/" + config["db_name"] + "/{library_name}_select_taxonomy.txt", library_name = LIBRARY_NAME),
					gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt",
					tax_added = config["rdir"] + "/tax_combined/euk_custom_post_derep_taxonomy.txt"
				output:
					custom_select = config["rdir"] + "/" + config["db_name"] + "/custom_euk_select_taxonomy.txt",
					euk_select = config["rdir"] + "/" + config["db_name"] + "/custom_euk_select_accessions.txt"
				params:
					script = config["wdir"] + "/scripts/coarse_selection_custom.R",
					rank = config["rank_custom"],
					nmax = config["nmax_coarse"],
					dbdir = config["rdir"] + "/" + config["db_name"],
					add = config["custom_ncbi_post_derep"],
					gendir = config["rdir"] + "/derep_combined",
					outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
				conda:
					config["wdir"] + "/envs/r.yaml"
				log:
					config["rdir"] + "/logs/select_hires_custom_ncbi.log"
				shell:
					"""
					cat {input.tax_select} > "{params.dbdir}/tmp_ncbi_select_taxonomy.txt"
					{params.script} -t "{params.dbdir}/tmp_ncbi_select_taxonomy.txt" -c {params.add} -r {params.rank} -n {params.nmax} -o {output.custom_select} &>>{log}
					mkdir -p {params.outdir}
					find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {output.custom_select}) | while read line
					do
						ln -sf "$line" {params.outdir}
					done
					cut -f1 {output.custom_select} | grep -F -f - {input.gen2taxid} > {output.euk_select}
					rm "{params.dbdir}/tmp_ncbi_select_taxonomy.txt"
					"""

# If the secondary flavour is set to "organelle_euk", check the collection of GTDB, organelles and eukaryotes
	if config["flavour_sec"] == "organelle_euk":
		rule check_hires_organelle_euk:
			input:
				select_gtdb = config["rdir"] + "/" + config["db_name"] + "/gtdb_select_accessions.txt",
				custom_pro = config["rdir"] + "/" + config["db_name"] + "/custom_pro_select_accessions.txt" if config["custom_gtdb_post_derep"] != "n" else [],
				select_macro = expand(config["rdir"] + "/" + config["db_name"] + "/{library_macro}_select_accessions.txt", library_macro = LIBRARY_MACRO),
				select_micro = expand(config["rdir"] + "/" + config["db_name"] + "/{library_micro}_select_accessions.txt", library_micro = LIBRARY_MICRO),
				select_organelle = config["rdir"] + "/" + config["db_name"] + "/organelle_select_accessions.txt",
				macro_linked = expand(config["rdir"] + "/" + config["db_name"] + "/genomes/{library_macro}_done", library_macro = LIBRARY_MACRO),
				micro_linked = expand(config["rdir"] + "/" + config["db_name"] + "/genomes/{library_micro}_done", library_micro = LIBRARY_MICRO),
				custom_euk = config["rdir"] + "/" + config["db_name"] + "/custom_euk_select_accessions.txt" if config["custom_ncbi_post_derep"] != "n" else []
			output:
				tax = config["rdir"] + "/" + config["db_name"] + "/select_taxonomy.txt",
				select = config["rdir"] + "/" + config["db_name"] + "/select_accessions.txt",
				checked = config["rdir"] + "/" + config["db_name"] + "/genomes/done"
			params:
				flav_dir = config["rdir"] + "/" + config["db_name"],
				gtdb_tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
				organelle_tax = config["rdir"] + "/tax_combined/organelle_derep_taxonomy.txt",
				microeuk_tax = config["rdir"] + "/tax_combined/{library_micro}_derep_taxonomy.txt",
				macro_euk_tax = config["rdir"] + "/" + config["db_name"] + "/{library_macro}_select_taxonomy.txt",
				custom_pro_tax = config["rdir"] + "/tax_combined/pro_custom_post_derep_taxonomy.txt" if config["custom_gtdb_post_derep"] != "n" else [],
				custom_euk_tax = config["rdir"] + "/" + config["db_name"] + "/custom_euk_select_taxonomy.txt" if config["custom_ncbi_post_derep"] != "n" else [],
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			shell:
				"""
				cat {input.select_gtdb} {input.custom_pro} \
					{input.select_macro} {input.select_micro} \
					{input.select_organelle} {input.custom_euk} > {output.select}

				cat {params.gtdb_tax} {params.organelle_tax} \
					{params.microeuk_tax} {params.macroeuk_tax} \
					{params.custom_pro_tax} {params.custom_euk_tax}  > {output.tax}

				if [[ $(cat {output.select} | wc -l) == $(find {params.outdir} -name '*.gz' | wc -l) ]]; then
					touch {output.checked}
				fi

				rm {params.flav_dir}/*_select_accessions.txt {params.flav_dir}/*_select_taxonomy.txt
				"""

# If the secondary flavour is set to "organelle_euk_virus", collect the viral accessions/genomes (results from derepG)
	if config["flavour_sec"] == "organelle_euk_virus":
		rule select_viruses:
			input:
				checkv = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt",
				added_nuc_checkv = config["rdir"] + "/tax_combined/vir_custom_post_derep_taxonomy.txt" if config["custom_checkv_post_derep"] != "n" else [],
				gen2taxid = config["rdir"] + "/tax_combined/full_genome2taxid.txt"
			output:
				vir_select = config["rdir"] + "/" + config["db_name"] + "/vir_select_accessions.txt"
			shell:
				"""
				cat {input.checkv} {input.added_nuc_checkv} | cut -f1 | grep -F -f - {input.gen2taxid} > {output.vir_select}
				"""

		rule collect_viruses:
			input:
				vir_select = config["rdir"] + "/" + config["db_name"] + "/vir_select_accessions.txt"
			output:
				linked = config["rdir"] + "/" + config["db_name"] + "/genomes/pro_done"
			params:
				gendir = config["rdir"] + "/derep_combined",
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			shell:
				"""
				mkdir -p {params.outdir}
				find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {input.vir_select}) | while read line
				do
					ln -sf "$line" {params.outdir}
				done
				touch {output.linked}
				"""

# check the collection of GTDB, organelles, eukaryotes and viruses
		rule check_hires_organelle_euk_virus:
			input:
				select_gtdb = config["rdir"] + "/" + config["db_name"] + "/gtdb_select_accessions.txt",
				custom_pro = config["rdir"] + "/" + config["db_name"] + "/custom_pro_select_accessions.txt" if config["custom_gtdb_post_derep"] != "n" else [],
				select_macro = expand(config["rdir"] + "/" + config["db_name"] + "/{library_macro}_select_accessions.txt", library_macro = LIBRARY_MACRO),
				select_micro = expand(config["rdir"] + "/" + config["db_name"] + "/{library_micro}_select_accessions.txt", library_micro = LIBRARY_MICRO),
				select_organelle = config["rdir"] + "/" + config["db_name"] + "/organelle_select_accessions.txt",
				macro_linked = expand(config["rdir"] + "/" + config["db_name"] + "/genomes/{library_macro}_done", library_macro = LIBRARY_MACRO),
				micro_linked = expand(config["rdir"] + "/" + config["db_name"] + "/genomes/{library_micro}_done", library_micro = LIBRARY_MICRO),
				custom_euk = config["rdir"] + "/" + config["db_name"] + "/custom_euk_select_accessions.txt" if config["custom_ncbi_post_derep"] != "n" else [],
				select_vir = config["rdir"] + "/" + config["db_name"] + "/vir_select_accessions.txt"
			output:
				tax = config["rdir"] + "/" + config["db_name"] + "/select_taxonomy.txt",
				select = config["rdir"] + "/" + config["db_name"] + "/select_accessions.txt",
				checked = config["rdir"] + "/" + config["db_name"] + "/genomes/done"
			params:
				flav_dir = config["rdir"] + "/" + config["db_name"],
				gtdb_tax = config["rdir"] + "/tax_combined/gtdb_derep_taxonomy.txt",
				organelle_tax = config["rdir"] + "/tax_combined/organelle_derep_taxonomy.txt",
				microeuk_tax = config["rdir"] + "/tax_combined/{library_micro}_derep_taxonomy.txt",
				macro_euk_tax = config["rdir"] + "/" + config["db_name"] + "/{library_macro}_select_taxonomy.txt",
				checkv_tax = config["rdir"] + "/tax_combined/checkv_derep_taxonomy.txt",
				custom_pro_tax = config["rdir"] + "/tax_combined/pro_custom_post_derep_taxonomy.txt" if config["custom_gtdb_post_derep"] != "n" else [],
				custom_euk_tax = config["rdir"] + "/" + config["db_name"] + "/custom_euk_select_taxonomy.txt" if config["custom_ncbi_post_derep"] != "n" else [],
				custom_vir_tax = config["rdir"] + "/tax_combined/vir_custom_post_derep_taxonomy.txt" if config["custom_checkv_post_derep"] != "n" else [],
				outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
			shell:
				"""
				cat {input.select_gtdb} {input.custom_pro} \
					{input.select_macro} {input.select_micro} \
					{input.select_organelle} {input.custom_euk} \
					{input.select_vir} > {output.select}

				cat {params.gtdb_tax} {params.organelle_tax} {params.checkv_tax} \
					{params.microeuk_tax} {params.macroeuk_tax} \
					{params.custom_pro_tax} {params.custom_euk_tax} {params.custom_vir_tax} > {output.tax}

				if [[ $(cat {output.select} | wc -l) == $(find {params.outdir} -name '*.gz' | wc -l) ]]; then
					touch {output.checked}
				fi

				rm {params.flav_dir}/*_select_accessions.txt {params.flav_dir}/*_select_taxonomy.txt
				"""

if config["flavour_main"] == "user":
	rule collect_genomes:
		input:
			custom_select = config["flavour_custom"]
		output:
			linked = config["rdir"] + "/" + config["db_name"] + "/genomes/done"
		params:
			gendir = config["rdir"] + "/derep_combined",
			outdir = config["rdir"] + "/" + config["db_name"] + "/genomes/"
		shell:
			"""
			mkdir -p {params.outdir}
			find {params.gendir} -name '*.gz' | grep -F -f <(cut -f1 {input.custom_select}) | while read line
			do
			  ln -sf "$line" {params.outdir}
			done
			if [[ $(cat {input.custom_select | wc -l}) == $(find {params.outdir} -name '*.gz' | wc -l) ]]
			then
			  touch {output.linked}
			fi
			"""
