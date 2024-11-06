version 1.0

workflow runCentrolign {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Run centrolign for multiple samples, on a single chromosome"
    }
    call centrolign

    output {
        File centrolignGFA=centrolign.centrolignGFA
    }
}

task centrolign {
    input {
      File centrolignOptions=""
      File newickGuideTree
      File assemblyFasta
      String chrom

      String sampleID
      String dockerImage="miramastoras/centrolign:v0.2.1"

      Int memSizeGB=256
      Int threadCount=8
      Int diskSizeGB=256
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        centrolign ~{centrolignOptions} \
            -T ~{newickGuideTree} \
            ~{assemblyFasta} \
            > ~{sampleID}.~{chrom}.gfa
    >>>
    output {
        File centrolignGFA=glob("*.gfa")[0]
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
