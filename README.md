# About

This is currently working only on the BYU RC cluster, with the 20kb probe sequences, starting from the filtered ortholog files (i.e. alignment is the first step).

# Download

	git clone https://github.com/paytoncarter14/odonata-nextflow.git
	cd odonata-nextflow

# Create run directory and taxa.txt file

Create an empty folder in the `runs` directory, then create a `taxa.txt` file with the identifiers that will be searched in the locus files. For example, we'll make a run named `chlorocyphidae` that will include all of the family Chlorocyphidae, all of the genus _Calopteryx_, and the single sample `GEODE7884_Euphaeidae_Euphaea_guerini`:

	mkdir runs/chlorocyphidae && cd $_
	nano taxa.txt

Then the content of taxa.txt (the search is case sensitive):

	Chlorocyphidae
	Calopterygidae_Calopteryx
	GEODE7884_Euphaeidae_Euphaea_guerini

# Activate Nextflow module

	module load jdk/1.20
	module load nextflow/23.10

If you want Nextflow to load automatically when you login, also run `module save`

# Start your Nextflow run

From your run directory:

	nextflow run ../../main.nf

When the pipeline finished, your tree file and all intermediate files will be in the `output` folder in your run directory.

If your terminal exits for any reason, the pipeline will stop. You can resume the pipeline with cached output files:

	nextflow run -resume ../../main.nf

To ensure your pipeline will not be stopped if your terminal exits, you can start a [screen](https://www.gnu.org/software/screen/) session and submit the Nextflow command to the Slurm cluster:

	screen
 	srun --time=72:00:00 --mem=8G nextflow run ../../main.nf

# More runs

To use the pipeline with another taxa sampling, simply create another run folder and `taxa.txt` file:

    odonata-nextflow
    |-- 20kb_locus_list.txt
    |-- README.md
    |-- main.nf
    |-- nextflow.config
    `-- runs
        |-- chlorocyphidae
        |   |-- taxa.txt
        |   `-- output
        |       |-- tree.treefile
        |       `-- ...
        `-- another-run
            |-- taxa.txt
            `-- output
                |-- tree.treefile
                `-- ...


