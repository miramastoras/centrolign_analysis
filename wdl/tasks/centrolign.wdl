version 1.0

workflow runCentrolign {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Run centrolign for multiple samples, on a single chromosome"
    }
    call centrolign

    output {
        File centrolignGFA=
    }
}

task centrolign {
    input {
      File CenSatBedFile
      File asmToRefPaf
      File assemblyFasta
      File AsHorSFBedFile

      String sampleID
      String dockerImage="miramastoras/centromere_scripts:v0.1"

      Int memSizeGB=16
      Int threadCount=8
      Int diskSizeGB=64
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        
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
