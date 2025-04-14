#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include modules
include { EXTRACT_INFO } from '../modules/extract_info'
include { GATHER_INFO_LISTS } from '../modules/gather_info_lists'
include { RUN_R_SELECT } from '../modules/run_r_select'
include { FILTER_VCF } from '../modules/filter_vcf'
include { MERGE_VCF } from '../modules/merge_vcf'

workflow POST_IMPUTATION {
    take:
    ch_imputedVCFs // [[meta.id, meta.cohort, meta.chr, meta.target_cohort], imputed_vcf]
    ch_ordered_cohorts // List of cohorts in order from CSV

    main:

    ch_extractedINFO = EXTRACT_INFO(ch_imputedVCFs)
    
    ch_extractedINFO
        .map { meta, info_file -> [ meta.cohort, meta.target_cohort, [meta, info_file] ] } 
        .groupTuple(by: [0]) // [ 'A', ['X', 'X', 'Y'], [ [meta1, file1], [meta2, file2], [meta3, file3] ] ]
        .map { cohort, target_cohort, meta_info_list ->
            def info_paths = meta_info_list.collect { it[1] } // Extract just the paths
            def target_cohort_gather = target_cohort.unique() 
            return [ cohort, target_cohort_gather, info_paths ]
        }
        .set { ch_gather_input } 
    
    // ch_gather_input.view { "Gather input: ${it}" }

    ch_sumINFO = GATHER_INFO_LISTS(ch_gather_input, ch_ordered_cohorts)

    // Output of GATHER_INFO_LISTS is chr-independent, it just need cohort and target cohort, but we would keep chr_gather for the sake of pairing in filter step
    
    ch_listRout = RUN_R_SELECT(ch_sumINFO)
    
    ch_filter_files_typed = ch_listRout.filterlists.map { cohort, target_cohort_gather, f_list -> [cohort, target_cohort_gather, f_list, 'filterout'] }
    ch_keep_files_typed   = ch_listRout.keeplists.map { cohort, target_cohort_gather, k_list -> [cohort, target_cohort_gather, k_list, 'keeplist'] }
    ch_all_filter_lists = ch_filter_files_typed.mix(ch_keep_files_typed)
    
    // But the output of RUN_R_SELECT is chr-dependent

    // Now handle the filter lists with cohort matching
    ch_filter_lists_for_join = ch_all_filter_lists
        .flatMap { cohort, target_cohort_gather, list_paths, list_type ->
            // Handle the case where list_paths might be a collection or a single path
            def paths = list_paths instanceof Collection ? list_paths : [list_paths]
            def results = []
            
            paths.each { path ->
                def basename = path.getBaseName()
                // def matchedCohort = null

                def chrMatch = basename =~ /chr(\d+)/
                def matchedChr = chrMatch ? chrMatch[0][1].toString() : null
 
                def matchedCohort = target_cohort_gather.find { targetCohort -> 
                    basename.contains(targetCohort) 
                }
                

                if (matchedCohort && matchedChr) {
                    log.info "Found match: cohort=${cohort}, targetCohort=${matchedCohort}, chr=${matchedChr}, file=${path.getName()}"
                    results.add([cohort, matchedCohort, matchedChr, path, list_type])
                } else {
                    log.warn "Could not match either or both cohort and chr in filter list: ${basename}"
                }
            }
            
            return results
        }
    
    // ch_filter_lists_for_join.view{ "Filter lists for join: ${it}" }  //['A', 'X', file for X]

    ch_imputed_vcfs_for_join = ch_imputedVCFs.map { meta, vcf -> 
        [ meta.cohort,
          meta.target_cohort,
          meta.chr.toString(),
          meta, vcf ] 
    }

    // Join imputed VCFs with their corresponding filter list using both imputed cohort and target cohort
    ch_to_filter = ch_imputed_vcfs_for_join
        .join(ch_filter_lists_for_join, by: [0, 1, 2])
        .map { cohort, target_cohort, chr, meta, vcf, filter_path, filter_type -> 
            [ meta, vcf, filter_path ]
        }
    // ch_to_filter.view { "VCF to filter: ${it}" }

    FILTER_VCF(ch_to_filter)
    
    FILTER_VCF.out.filtered_vcf
        .collect(flat: false) { meta, vcf -> 
            def new_meta = [
                cohort: meta.cohort,
                chr: meta.chr,
                id: meta.id
            ]    
            return [new_meta, vcf] 
        }
        .flatMap()
        .groupTuple()        
        .set { ch_to_merge }
        
    MERGE_VCF(ch_to_merge)

    emit:
    postimputed_vcfs = MERGE_VCF.out.merged_vcf
    postimputed_idx = MERGE_VCF.out.merged_idx

}