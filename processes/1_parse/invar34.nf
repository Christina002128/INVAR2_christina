// tapasSetting is a closure that can be used inside more than one process that
// gives the full string for the TAPAS setting.

tapasSetting =
{
    return "${params.ERROR_SUPPRESSION_NAME}_BQ_${params.BASEQ}.MQ_${params.MAPQ}"
}

process createMutationsTable
{
    memory '8g'
    cpus 1
    time '1h'

    publishDir 'mutations', mode: 'link'

    input:
        path mutationsFile
        path tumourMutationsFile
        path layoutFile

    output:
        path 'mutation_table.all.rds', emit: "allMutationsFile", optional: true
        path 'mutation_table.filtered.rds', emit: "filteredMutationsFile"

    shell:

        """
        Rscript --vanilla "!{projectDir}/R/1_parse/createMutationsTable.R" \
            --mutations="!{mutationsFile}" \
            --tumour-mutations="!{tumourMutationsFile}" \
            --layout="!{layoutFile}" \
            --tapas="!{tapasSetting}" \
            --cosmic-threshold=!{params.cosmic_threshold} \
            --mqsb-threshold=!{params.individual_MQSB_threshold} \
            --max-depth=!{params.max_depth} \
            --min-ref-depth=!{params.min_ref_depth} \
            --alt-alleles-threshold=!{params.n_alt_alleles_threshold} \
            --minor-alt-allele-threshold=!{params.minor_alt_allele_threshold}
        """
}

process offTargetErrorRates
{
    memory '4g'
    cpus 2
    time '1h'

    publishDir 'off_target', mode: 'link'

    input:
        path mutationsFile
        path layoutFile

    output:
        path 'locus_error_rates.off_target.rds', emit: 'locusErrorRates'
        path 'mutation_table.error_rates.cosmic.rds', emit: "cosmicErrorRates"
        path 'mutation_table.error_rates.no_cosmic.rds', emit: "noCosmicErrorRates"
        path '*.tsv', optional: true

    shell:
        """
        Rscript --vanilla "!{projectDir}/R/1_parse/offTargetErrorRates.R" \
            --mutations="!{mutationsFile}" \
            --layout="!{layoutFile}" \
            --control-proportion=!{params.proportion_of_controls} \
            --max-background-allele-frequency=!{params.max_background_mean_allele_frequency} \
            !{params.is_bloodspot ? "--bloodspot" : ""}
        """
}

process createOnTargetMutationsTable
{
    memory '4g'
    cpus 1
    time '1h'

    input:
        path mutationsFile
        path tumourMutationsFile
        path layoutFile
        path errorRatesFile

    output:
        path 'mutation_table.on_target.all.rds', emit: "onTargetMutationsFile"
        path '*.tsv', optional: true

    shell:
        """
        Rscript --vanilla "!{projectDir}/R/1_parse/createOnTargetMutationsTable.R" \
            --mutations="!{mutationsFile}" \
            --tumour-mutations="!{tumourMutationsFile}" \
            --layout="!{layoutFile}" \
            --error-rates="!{errorRatesFile}"
        """
}

process onTargetErrorRatesAndFilter
{
    memory '4g'
    cpus 1
    time '1h'

    publishDir 'on_target', mode: 'link'

    input:
        path mutationsFile
        path tumourMutationsFile
        path layoutFile

    output:
        path 'locus_error_rates.on_target.rds', emit: 'locusErrorRates'
        path 'mutation_table.on_target.rds', emit: "onTargetMutationsFile"
        path '*.tsv', optional: true
        path '*.pdf', optional: true

    shell:
        """
        Rscript --vanilla "!{projectDir}/R/1_parse/onTargetErrorRatesAndFilter.R" \
            --mutations="!{mutationsFile}" \
            --layout="!{layoutFile}" \
            --study="!{params.STUDY_ID}" \
            --tapas="!{tapasSetting}" \
            --cosmic-threshold=!{params.cosmic_threshold} \
            --control-proportion=!{params.proportion_of_controls} \
            --max-background-allele-frequency=!{params.max_background_mean_allele_frequency} \
            !{params.is_bloodspot ? "--bloodspot" : ""} \
            --allele-frequency-threshold=!{params.allele_frequency_threshold}
        """
}

workflow invar34
{
    take:
        mutationsChannel
        tumourMutationsChannel
        layoutChannel

    main:
        createMutationsTable(mutationsChannel, tumourMutationsChannel, layoutChannel)

        offTargetErrorRates(createMutationsTable.out.filteredMutationsFile, layoutChannel)

        createOnTargetMutationsTable(
            createMutationsTable.out.filteredMutationsFile,
            tumourMutationsChannel,
            layoutChannel,
            offTargetErrorRates.out.noCosmicErrorRates)

        onTargetErrorRatesAndFilter(createOnTargetMutationsTable.out.onTargetMutationsFile,
                                    tumourMutationsChannel,
                                    layoutChannel)

    emit:
        onTargetMutationFile = onTargetErrorRatesAndFilter.out.onTargetMutationsFile
}
