        
def gmapFinder(tool, gmap) {

    def pattern = null
    
    if (tool == "shapeit5") {
        pattern = /chr(\d+|X|X_par1|X_par2)\.b37\.gmap\.gz/

    } else if (tool == "beagle5") {
        pattern = /plink\.chr(\d+|X|X_par1|X_par2)\.GRCh37\.map/
    }
    else {
        error "Unsupported tool: ${tool}. Supported tools are 'shapeit5' and 'beagle5'."
    }

    def ref_map_chr = Channel.fromPath(["${gmap}/*"], type : 'file')
                            .map{ file ->
                            def match = file.name =~ pattern
                            if (match) {
                                def chromosome = match[0][1]
                                return [chromosome, file]
                            } else {
                                error "File name does not match expected pattern: ${file.name}"
                            }
                            }
                            .filter { it[0] in params.chromosomes*.toString() }
    return ref_map_chr
}