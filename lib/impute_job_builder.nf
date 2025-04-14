// utils/impute_job_builder.nf

def pairImputeJobs(ch_phasedVCFs, ch_ref_raw) {
        // Create a channel with reference panels by cohort and chromosome
        def ch_ref_panels = ch_ref_raw
            .map { meta, panel -> [meta.cohort, meta.chr, meta, panel] }
        
        // Prepare phased VCFs for combination
        def ch_phased_for_combo = ch_phasedVCFs
            .map { meta, vcf, index -> [meta.cohort, meta.chr, meta, vcf, index] }
        
        // Create all valid combinations of phased VCFs and reference panels
        def ch_impute_jobs = ch_phased_for_combo
            .combine(ch_ref_panels, by: 1)  // Join by chromosome first
            .filter { chr, phased_cohort, phased_meta, phased_vcf, phased_index, ref_cohort, ref_meta, ref_panel ->
                // Only impute when cohorts are different
                phased_cohort != ref_cohort
            }
            .map { chr, phased_cohort, phased_meta, phased_vcf, phased_index, ref_cohort, ref_meta, ref_panel ->
                // Create a new meta with target_cohort information required by post_imputation
                def newMeta = phased_meta.clone()
                newMeta.target_cohort = ref_cohort
                
                // Return format that matches IMPUTE_MINIMAC4 input: [meta, vcf, index, panel]
                [newMeta, phased_vcf, phased_index, ref_panel]
            }


        return ch_impute_jobs
}       