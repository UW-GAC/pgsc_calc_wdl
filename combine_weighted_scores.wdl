version 1.0

workflow combine_weighted_scores {
    input {
        Array[File] scorefiles
        File weights
    }

    scatter (scorefile in scorefiles) {
        call weight_scores {
            input:
                scorefile = scorefile,
                weights = weights
        }
    }
}


task weight_scores {
    input {
        File scorefile
        File weights
        Int mem_gb = 64
    }

    Int disk_size = ceil(3*(size(scorefile, "GB") + size(weights, "GB"))) + 10

    command <<<
        R << RSCRIPT
        library(tidyverse)
        scores <- read_tsv("~{scorefile}")
        weights <- read_tsv("~{weights}")
        selected_scores <- scores %>%
            select(ID, effect_allele, any_of(weights$score)) %>%
            filter(!if_all(any_of(weights$score), ~ . == 0))
        for (i in 1:nrow(weights)) {
            if (weights$score[i] %in% names(selected_scores)) {
                selected_scores[[weights$score[i]]] <- selected_scores[[weights$score[i]]] * weights$weight[i]
            }
        }
        write_tsv(selected_scores, "weighted_scores.txt")
        RSCRIPT
    >>>

    output {
        File weighted_scores = "weighted_scores.txt"
    }

    runtime {
        docker: "rocker/tidyverse:4"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
    }
}
