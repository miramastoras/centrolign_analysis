version 1.0

import "./extract_hors.wdl" as extract_hors_wf

workflow extract_hor_sequence {
    input {
        File horArrayBed
        File asmToRefPaf
        File assemblyFasta

        String sampleID
    }

    call extract_hors_wf.extract_hor_sequence as extract_seq {
        input:
            horArrayBed  = horArrayBed,
            asmToRefPaf  = asmToRefPaf,
            assemblyFasta = assemblyFasta,
            sampleID = sampleID
    }

   output {
        Array[File] horFastas =extract_seq.horFastas
    }
}
