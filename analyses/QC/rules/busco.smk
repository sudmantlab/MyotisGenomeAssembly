rule busco:
	version: "5.2.2"  # augustus takes way too long - use METAEUK
	input:
		fasta = "data/genomes/{genome}.fa"
	params:
		prefix = "{genome}",
		busco_path = "data/funannotate_database",
		busco_set = "mammalia",
		augustus_species = "homo_sapiens",
		wd = os.getcwd(),
	threads: 20
	singularity:
		"docker://ezlabgva/busco:v5.2.2_cv1"
	log:
		stdout = "logs/BUSCO/{genome}.stdout.txt",
		stderr = "logs/BUSCO/{genome}.stderr.txt"
	output:
		buscos = directory("output/BUSCO/{genome}"),
	shell:
		"""
		echo -e "\n$(date)\tStarting on host: $(hostname) ...\n"
		
		PROJDIR="$(pwd -P)"

		if [[ ! -d results/{params.prefix}/BUSCO ]]
		then
			mkdir -p output/BUSCO/{wildcards.genome}
		fi
		cd output/BUSCO/

		if [[ ! -d tmp ]]
		then
			mkdir tmp
		fi

		#cp -rf /usr/share/augustus/config tmp/config
		#AUGUSTUS_CONFIG_PATH=$(pwd)/tmp/config

                echo "$(pwd -P)"

                #busco --list-datasets

		#run BUSCO
		busco -f \\
		  --in $PROJDIR/{input.fasta} \\
		  --out {wildcards.genome} \\
		  -l $PROJDIR/{params.busco_path}/{params.busco_set} \\
		  --mode genome \\
		  -c {threads} \\
		  1> $PROJDIR/{log.stdout} \\
		  2> $PROJDIR/{log.stderr}
		
		rm -r tmp

		echo -e "\n$(date)\tFinished!\n"
		"""
