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

    call combine_scores {
        input:
            scorefiles = weight_scores.weighted_scores
    }

    call format_scores {
        input:
            scorefile = combine_scores.combined_scores
    }

    output {
        File combined_scores = combine_scores.combined_scores
        File formatted_scores = format_scores.formatted_scores
    }
}


task weight_scores {
    input {
        File scorefile
        File weights
        Int mem_gb = 64
    }

    Int disk_size = ceil(10*(size(scorefile, "GB") + size(weights, "GB"))) + 10

    command <<<
        set -e
        mv ~{scorefile} scorefile.txt.gz
        gunzip scorefile.txt.gz
        wc -l scorefile.txt

        R --vanilla << RSCRIPT
        # More debugging
        options(datatable.verbose=TRUE)
        library(tidyverse)
        install.packages("R.utils", repos="https://cloud.r-project.org")
        weight_file <- "~{weights}"
        score_file <- "scorefile.txt"
        weights <- read_tsv(weight_file)
        score_vars_head <- read_tsv(score_file, n_max=10)
        pgs <- intersect(names(score_vars_head), weights[["score"]])
        cols <- c("ID", "effect_allele", pgs)
        scores <- data.table::fread(score_file, select=cols)
        print(dim(scores))
        # Debugging
        print(sapply(scores, function(x) sum(x != 0)))
        selected_scores <- scores %>%
            filter(!if_all(any_of(pgs), ~ . == 0))
        print(sapply(selected_scores, function(x) sum(x != 0)))
        rm(scores)
        for (i in 1:nrow(weights)) {
            if (weights[["score"]][i] %in% names(selected_scores)) {
                selected_scores[[weights[["score"]][i]]] <- selected_scores[[weights[["score"]][i]]] * weights[["weight"]][i]
            }
        }
        print(sapply(selected_scores, function(x) sum(x != 0)))
        write_tsv(selected_scores, "weighted_scores.txt")
        RSCRIPT
    >>>

    output {
        File weighted_scores = "weighted_scores.txt"
    }

    runtime {
        docker: "rocker/tidyverse:4.5.1"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
    }
}


task combine_scores {
    input {
        Array[File] scorefiles
        Int mem_gb = 64
    }

    Int disk_size = ceil(3*(size(scorefiles, "GB"))) + 10

    command <<<
        R << RSCRIPT
        library(tidyverse)
        weighted_scorefiles <- readLines("~{write_lines(scorefiles)}")
        all_scores <- read_tsv(weighted_scorefiles[1])
        for (i in 2:length(weighted_scorefiles)) {
            weighted_scores <- read_tsv(weighted_scorefiles[i])
            chk_alleles <- inner_join(all_scores[,1:2], weighted_scores[,1:2], by=c("ID"))
            stopifnot(all(chk_alleles[["effect_allele.x"]] == chk_alleles[["effect_allele.y"]]))
            all_scores <- full_join(all_scores, weighted_scores, by=c("ID", "effect_allele"))
        }
        all_scores[is.na(all_scores)] <- 0
        snps <- all_scores[,1:2]
        scores <- rowSums(all_scores[,3:ncol(all_scores)])
        new_scores <- bind_cols(snps, score=scores)
        write_tsv(new_scores, "combined_scores.txt")
        RSCRIPT
    >>>

    output {
        File combined_scores = "combined_scores.txt"
    }

    runtime {
        docker: "rocker/tidyverse:4.5.1"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
    }
}


task format_scores {
    input {
        File scorefile
        Int mem_gb = 64
    }

    Int disk_size = ceil(3*(size(scorefile, "GB"))) + 10

    command <<<
        R << RSCRIPT
        library(tidyverse)
        scores <- read_tsv("~{scorefile}")
        scores <- scores %>%
            separate_wider_delim(ID, delim=":", names=c("chr_name", "chr_position", "ref", "alt"), cols_remove=FALSE) %>%
            mutate(other_allele = ifelse(effect_allele == alt, ref, alt)) %>%
            rename(effect_weight = score) %>%
        select(ID, chr_name, chr_position, effect_allele, other_allele, effect_weight)
        write_tsv(scores, "formatted_scores.txt")
        RSCRIPT
    >>>

    output {
        File formatted_scores = "formatted_scores.txt"
    }

    runtime {
        docker: "rocker/tidyverse:4"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
    }
}
