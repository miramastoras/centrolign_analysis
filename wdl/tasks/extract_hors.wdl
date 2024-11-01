version 1.0

workflow extract_hors {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "extracts HOR arrays for each chromosome in input assembly"
    }

    input {
        File assemblyFasta
        File asmToRefPaf

        File AsHorSFBedFile
        File CenSatBedFile

        String sampleName

        String dockerImage="miramastoras/centromere_scripts:v0.1"
    }

    ## use alignment to CHM13 to assign HORs to chromosomes
    ## filters out discontiguous HORs, cross-chromosomal groups not all fully resolved

    call locate_hors {
        input:
            CenSatBedFile=CenSatBedFile,
            asmToRefPaf=asmToRefPaf,
            assemblyFasta=assemblyFasta,
            AsHorSFBedFile=AsHorSFBedFile,
            sampleID=sampleName
    }

    output {

    }
}

task locate_hors {
    input {
      File CenSatBedFile
      File asmToRefPaf
      File assemblyFasta
      File AsHorSFBedFile

      String sampleID
      String dockerImage
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        ASM_FILENAME=$(basename -- "~{assemblyFasta}")

        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} ~{sampleID}.fasta.gz
            gunzip -f ~{sampleID}.fasta.gz
        else
            cp ~{assemblyFasta} ~{sampleID}.fasta
        fi

        samtools faidx ~{sampleID}.fasta

        python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/locate_hors_from_censat.py \
            -c ~{CenSatBedFile} \
            -a ~{AsHorSFBedFile} \
            -p ~{asmToRefPaf} \
            -f ~{sampleID}.fasta \
            > ~{sampleID}_hor_arrays.bed
    >>>
    output {
        File horArrayBed=glob("*hor_arrays.bed")[0]
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
}
