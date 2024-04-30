taxa_sql = "${launchDir}/taxa.sql"
filtered_ortholog_dir = '/nobackup/scratch/grp/grp_geode/filtered_orthologs'
locus_list_20kb = '/nobackup/scratch/grp/grp_geode/payton/nextflow/20kb_locus_list.txt'

process sql_to_sample_ids {
    output: path('taxa_sample_ids.txt')

    executor 'local'

    shell: '''sqlite3 !{projectDir}/geode.db < !{taxa_sql} > taxa_sample_ids.txt'''
}

process make_locus_list {
    output: path('locus_list.txt')
    shell: '''ls !{filtered_ortholog_dir}/probe_orthologs/*.fasta > locus_list.txt'''
}

process select_taxa_from_filtered_orthologs {
    input:
        path('taxa_sample_ids.txt')
        val(locus)
    output: path('L*.fasta')

    cpus 1
    memory 16.GB
    time 15.minute
    errorStrategy 'ignore'

    shell: '''grep --no-group-separator -h -A1 -f taxa_sample_ids.txt !{filtered_ortholog_dir}/probe_orthologs/!{locus}.fasta > $(basename !{locus}).fasta'''
}

process mafft_alignment {
    input: path(locus)
    output: path 'L*.aligned.fasta'

    cpus 4
    memory {16.GB * task.attempt}
    time {4.hour * task.attempt}

    shell: '''!{projectDir}/bin/mafft-linux64/mafft.bat --maxiterate 1000 --localpair --thread !{task.cpus} --adjustdirectionaccurately !{locus} | sed "s/|/__/g" > !{locus.simpleName}.aligned.fasta'''
}

process aliscore {
    input: path(locus)
    output: tuple path(locus), path("*.fasta.profiles.svg"), path("*.fasta.tre"), path("*.fasta_List_l_all.txt"), path("*.fasta_Profile_l_all.txt")

    errorStrategy { task.exitStatus == 255 ? 'ignore' : 'retry' } // 255 means not enough taxa for tree reconstruction

    cpus 1
    memory {16.GB * task.attempt}
    time {4.hour * task.attempt}

    shell: '''Aliscore.02.2.pl -i !{locus}'''
}

process alicut {
    input: tuple path(locus), path(profiles_svg), path(tre), path(list_txt), path(profile_txt)
    output: tuple path("*.alicut.fasta"), path("*.alicut_info.xls")

    cpus 1
    memory {8.GB * task.attempt}
    time {1.hour * task.attempt}

    shell: '''
        if [[ -z $(grep '[^[:space:]]' *.fasta_List_l_all.txt) ]]; then
            ln -s !{locus} !{locus.simpleName}.alicut.fasta
            touch !{locus.simpleName}.alicut_info.xls
        else
            ALICUT_V2.31.pl -s
            mv ALICUT_info.xls !{locus.simpleName}.alicut_info.xls
            mv ALICUT_L* !{locus.simpleName}.alicut.fasta
        fi
        sed -i -e 's/_R_//g' -e 's/^>L[0-9]\\{1,4\\}_/>/g' -e '/^>/s/_comp[0-9]*.*//' !{locus.simpleName}.alicut.fasta
    '''
}

process fasconcat {
    input: path(alicut_fastas)
    output: tuple path('FcC_info.xls'), path('partition_def.txt'), path('FcC_smatrix.fas')

    cpus 1
    memory {64.GB * task.attempt}
    time {2.hour * task.attempt}

    shell: '''
    FASconCAT_v1.11.pl -i -s
    extract_partition_file_from_fcc.py FcC_info.xls
    '''
}

process iqtree_model_selection {
    input: tuple path(fcc_info), path(partition), path(fcc_smatrix)
    output: tuple path(fcc_smatrix), path('partition.best_scheme.nex')

    cpus 8
    memory 64.GB
    time 72.hour
    shell: 'iqtree2 -s FcC_smatrix.fas -spp partition_def.txt -nt AUTO -safe -pre partition -m TESTMERGEONLY -mset GTR+G'
}

process iqtree_tree_search {
    input: tuple path(fcc_smatrix), path(best_scheme_nex)
    output: path("tree.*")

    cpus 8
    memory 64.GB
    time 72.hour

    shell: 'iqtree2 -s FcC_smatrix.fas -spp partition.best_scheme.nex -nt AUTO -safe -pre tree -m MFP -bb 1000 -bnni'
}

workflow {
    sql_to_sample_ids_c = sql_to_sample_ids()
    make_locus_list_c = Channel.fromPath(locus_list_20kb).splitText().map{it -> it.trim()}
    select_taxa_from_filtered_orthologs_c = select_taxa_from_filtered_orthologs(sql_to_sample_ids_c, make_locus_list_c)
    mafft_alignment_c = mafft_alignment(select_taxa_from_filtered_orthologs_c)
    aliscore_c = aliscore(mafft_alignment_c)
    alicut_c = alicut(aliscore_c)
    fasconcat_c = fasconcat(alicut_c.map{it[0]}.collect())
    //fasconcat_c = fasconcat(mafft_alignment_c.collect())
    iqtree_model_selection_c = iqtree_model_selection(fasconcat_c)
    iqtree_tree_search_c = iqtree_tree_search(iqtree_model_selection_c)
    
}