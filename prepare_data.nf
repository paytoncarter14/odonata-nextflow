/* process init_db {
    shell: '''
    rm -f stats.db
    echo "create table sample (sample_id text, total_reads_after_filter int, pct_q30_after_filter float, pct_gc_after_filter float, pct_duplication float, pct_all_bases_filtered float, pct_adapter_bases_filtered float);"
    '''
} */

process check_files {
    input: tuple val(file_id), val(sample_id)
    output: val reads

    exec:
    if (!sample_id.matches("[a-zA-Z0-9_]+")) {
        error("Sample IDs in key.txt can have only alphanumeric characters and underscores. Offending value: ${sample_id}")
    }
    r1 = "${launchDir}/fastq/${file_id}_R1.fastq.gz"
    if (!file(r1).exists()) {
        error("${r1} not found.")
    }
    r2 = "${launchDir}/fastq/${file_id}_R2.fastq.gz"
    if (!file(r2).exists()) {
        error("${r2} not found.")
    }
    reads = [sample_id, [r1, r2]]
}

process init_db {
    input: val reads

    shell: '''
    rm -f stats.db
    echo "" | sqlite3 stats.db
    '''
}


process fastp {
    input: val reads
    output: tuple path("*.fastp.R1.fastq.gz"), path("*.fastp.R2.fastq.gz"), path("*.fastp.json"), path("*.fastp.html")

    executor "slurm"
    time {1.hour * task.attempt}
    memory {8.GB * task.attempt}
    cpus 3
    
    shell: '''
    !{projectDir}/bin/fastp \
    -w !{task.cpus} \
    -i !{reads[1][0]} \
    -I !{reads[1][1]} \
    -o !{reads[0]}.fastp.R1.fastq.gz \
    -O !{reads[0]}.fastp.R2.fastq.gz \
    -j !{reads[0]}.fastp.json \
    -h !{reads[0]}.fastp.html \
    -g
    '''
}

process fastp_to_db {
    input: path(jsons)

    script: """
    #!${projectDir}/bin/venv/bin/python3

    import glob, orjson, sqlite3, os
    import pandas as pd

    out = {}

    for file in glob.glob("*.json"):
        sample_id = os.path.basename(file).split('.')[0]
        with open(file, 'r') as f:
            json = orjson.loads(f.read())
        out_dict = {}
        out_dict['total_reads_after_filter'] = json['summary']['after_filtering']['total_reads']
        out_dict['pct_q30_after_filter'] = json['summary']['after_filtering']['q30_rate']
        out_dict['pct_gc_after_filter'] = json['summary']['after_filtering']['gc_content']
        out_dict['pct_duplication'] = json['duplication']['rate']
        out_dict['pct_all_bases_filtered'] = 1 - (json['summary']['after_filtering']['total_bases'] / json['summary']['before_filtering']['total_bases'])
        out_dict['pct_adapter_bases_filtered'] = json['adapter_cutting']['adapter_trimmed_bases'] / json['summary']['before_filtering']['total_bases']
        out[sample_id] = out_dict
    pd.DataFrame(out).transpose().to_sql("fastp", sqlite3.connect("${launchDir}/stats.db"), if_exists="replace")
    """
}

workflow {
    Channel.fromPath(launchDir / "key.txt").splitCsv(sep: "\t") | check_files | fastp | map {it[2]} | collect | fastp_to_db
}