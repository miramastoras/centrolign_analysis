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

}
