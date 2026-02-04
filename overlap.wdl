version 1.0

import "calc_scores.wdl" as cs

workflow overlap {
    input {
        File scorefile
        File variants
        Int mem_gb = 64
    }

    call cs.compute_overlap {
        input:
            scorefile = scorefile,
            variants = variants,
            mem_gb = mem_gb
    }

    output {
        File overlap = compute_overlap.overlap
    }
}
