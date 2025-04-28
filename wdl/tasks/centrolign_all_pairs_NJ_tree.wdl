version 1.0

workflow runCentrolignAllPairsNJ {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Use centrolign all pairs alignments to create a NJ tree"
    }
    call infer_tree

    output {
        File NJ_Tree=infer_tree.NWK
    }
}

task infer_tree {
    input {
      File pairwiseCigarsTarGz
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

        python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_cigar/ \
        > /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr12_r2_centrolign_all_pairs_nj_tree.nwk
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
