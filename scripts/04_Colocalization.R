############################################################
# COLOCALIZATION ANALYSIS
# AKR1C2, MAPK3, MAP3K11, FADS1, BID
#
# Exposure: eQTLGen whole-blood cis-eQTL
# Outcome : CRC GWAS GCST90018808
#
# Main analysis:
# All SNPs with exact rsID + GRCh37 position + allele-pair match
############################################################


############################
# 0. Clean environment
############################

rm(list = ls())
gc()


############################
# 1. Packages
############################

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

if (!requireNamespace("coloc", quietly = TRUE)) {
  install.packages("coloc")
}

library(data.table)
library(coloc)


############################
# 2. Working directory
############################

setwd(
  "C:/Users/92066/Desktop/单细胞孟德尔/MR/colocalization"
)

cat("Working directory:\n")
print(getwd())


############################
# 3. Read three input files
############################

eqtl <- fread(
  "AKR1C2_MAPK3_MAP3K11_FADS1_BID_full_cis_eQTL.txt.gz"
)

af <- fread(
  "five_genes_AF.txt.gz"
)

crc <- fread(
  "five_genes_CRC_GWAS.txt"
)


############################
# 4. Verify required columns
############################

required_eqtl <- c(
  "Pvalue",
  "SNP",
  "SNPChr",
  "SNPPos",
  "AssessedAllele",
  "OtherAllele",
  "Zscore",
  "GeneSymbol",
  "NrSamples"
)

required_af <- c(
  "SNP",
  "hg19_chr",
  "hg19_pos",
  "AlleleA",
  "AlleleB",
  "AlleleB_all"
)

required_crc <- c(
  "CHR",
  "POS",
  "SNP",
  "REF",
  "ALT",
  "BETA",
  "SE",
  "LP",
  "EAF"
)

stopifnot(
  all(required_eqtl %in% names(eqtl)),
  all(required_af %in% names(af)),
  all(required_crc %in% names(crc))
)

cat("\nInput dimensions:\n")
cat("eQTL:", nrow(eqtl), "rows\n")
cat("AF  :", nrow(af),   "rows\n")
cat("CRC :", nrow(crc),  "rows\n")


# Confirm that the five intended genes are present in the eQTL input
expected_genes <- c(
  "AKR1C2",
  "MAPK3",
  "MAP3K11",
  "FADS1",
  "BID"
)

missing_genes <- setdiff(
  expected_genes,
  unique(eqtl$GeneSymbol)
)

if (length(missing_genes) > 0) {
  stop(
    paste(
      "Missing gene(s) in eQTL input:",
      paste(missing_genes, collapse = ", ")
    )
  )
}

cat("\nGenes detected in eQTL input:\n")
print(sort(intersect(expected_genes, unique(eqtl$GeneSymbol))))


############################
# 5. Basic SNP checks
############################

cat("\nDuplicated SNP records before matching:\n")

cat(
  "eQTL:",
  sum(duplicated(eqtl$SNP)),
  "\n"
)

cat(
  "AF:",
  sum(duplicated(af$SNP)),
  "\n"
)

cat(
  "CRC:",
  sum(duplicated(crc$SNP)),
  "\n"
)

# eQTL and AF must be unique
stopifnot(
  sum(duplicated(eqtl$SNP)) == 0,
  sum(duplicated(af$SNP)) == 0
)

# CRC can contain duplicated rsIDs due to multiallelic sites.
# We DO NOT simply remove duplicates here.


############################
# 6. Merge eQTL with allele-frequency data
############################

eqtl2 <- merge(
  eqtl,
  af[, .(
    SNP,
    hg19_chr,
    hg19_pos,
    AlleleA,
    AlleleB,
    AlleleB_all
  )],
  by = "SNP",
  all.x = TRUE,
  sort = FALSE
)

cat(
  "\neQTL + AF:",
  nrow(eqtl2),
  "rows\n"
)

stopifnot(
  nrow(eqtl2) == nrow(eqtl)
)


############################
# 7. Check GRCh37 positions
############################

eqtl2[, coord_af_match :=
        as.character(SNPChr) ==
        as.character(hg19_chr) &
        SNPPos == hg19_pos
]

cat("\neQTL vs AF coordinate matching:\n")
print(table(eqtl2$coord_af_match, useNA = "ifany"))

stopifnot(
  all(eqtl2$coord_af_match)
)


############################
# 8. Determine eQTL allele orientation
############################

eqtl2[, allele_status := fifelse(
  
  AssessedAllele == AlleleB &
    OtherAllele == AlleleA,
  
  "Assessed=AlleleB",
  
  fifelse(
    
    AssessedAllele == AlleleA &
      OtherAllele == AlleleB,
    
    "Assessed=AlleleA",
    
    "Mismatch"
  )
)]


cat("\neQTL allele orientation:\n")
print(
  table(
    eqtl2$allele_status,
    useNA = "ifany"
  )
)

# No unidentified allele configuration allowed
stopifnot(
  sum(eqtl2$allele_status == "Mismatch") == 0
)


############################
# 9. eQTL effect allele frequency
############################

# eQTL Z-score is referenced to AssessedAllele.
# Convert AF file frequency to frequency of AssessedAllele.

eqtl2[, EAF_eqtl := fifelse(
  
  allele_status == "Assessed=AlleleB",
  
  AlleleB_all,
  
  1 - AlleleB_all
)]


eqtl2[, MAF_eqtl :=
        pmin(
          EAF_eqtl,
          1 - EAF_eqtl
        )
]


stopifnot(
  all(is.finite(eqtl2$EAF_eqtl)),
  all(eqtl2$EAF_eqtl > 0),
  all(eqtl2$EAF_eqtl < 1),
  all(eqtl2$MAF_eqtl > 0),
  all(eqtl2$MAF_eqtl <= 0.5)
)

cat("\neQTL MAF summary:\n")
print(summary(eqtl2$MAF_eqtl))


############################
# 10. Merge with CRC GWAS
############################

# allow.cartesian=TRUE is deliberate:
# GCST90018808 contains some multiallelic rsIDs.

dat <- merge(
  eqtl2,
  crc,
  by = "SNP",
  all = FALSE,
  allow.cartesian = TRUE,
  sort = FALSE
)

cat(
  "\nRows immediately after eQTL-CRC rsID merge:",
  nrow(dat),
  "\n"
)


############################
# 11. Exact coordinate matching
############################

dat[, coord_match :=
      as.character(SNPChr) ==
      as.character(CHR) &
      SNPPos == POS
]


############################
# 12. Exact allele-pair matching
############################

dat[, allele_match :=
      
      (
        AssessedAllele == ALT &
          OtherAllele == REF
      ) |
      
      (
        AssessedAllele == REF &
          OtherAllele == ALT
      )
]


############################
# 13. Final matched master dataset
############################

dat2 <- dat[
  coord_match == TRUE &
    allele_match == TRUE
]


cat("\nAfter exact position + allele matching:\n")
cat("Rows        :", nrow(dat2), "\n")
cat("Unique SNPs :", uniqueN(dat2$SNP), "\n")
cat(
  "Duplicates  :",
  sum(duplicated(dat2$SNP)),
  "\n"
)

# Critical check: one row per SNP
stopifnot(
  nrow(dat2) == uniqueN(dat2$SNP),
  sum(duplicated(dat2$SNP)) == 0
)


cat("\nSNPs per gene:\n")
print(
  dat2[
    ,
    .N,
    by = GeneSymbol
  ][order(GeneSymbol)]
)


############################
# 14. Harmonize CRC beta
############################

# GCST90018808 VCF:
# BETA is effect of ALT allele.
#
# Convert CRC beta to effect of
# the same AssessedAllele used by eQTLGen.

dat2[, beta_crc_h := fifelse(
  
  ALT == AssessedAllele,
  
  BETA,
  
  -BETA
)]


############################
# 15. Harmonize CRC EAF
############################

# CRC EAF = ALT allele frequency.
# Convert to AssessedAllele frequency.

dat2[, EAF_crc_h := fifelse(
  
  ALT == AssessedAllele,
  
  EAF,
  
  1 - EAF
)]


dat2[, MAF_crc :=
       pmin(
         EAF_crc_h,
         1 - EAF_crc_h
       )
]


############################
# 16. Allele-frequency discrepancy
############################

dat2[, AF_diff :=
       abs(
         EAF_eqtl -
           EAF_crc_h
       )
]

cat("\nAbsolute EAF difference:\n")
print(summary(dat2$AF_diff))

cat("\nNumber with |Delta EAF| > 0.20:\n")

print(
  dat2[
    AF_diff > 0.20,
    .N,
    by = GeneSymbol
  ][order(GeneSymbol)]
)


############################
# 17. Reconstruct eQTL beta and SE
############################

# eQTLGen standardized expression:
#
# beta =
# Z / sqrt[2*p*(1-p)*(N + Z^2)]
#
# SE =
# 1 / sqrt[2*p*(1-p)*(N + Z^2)]
#
# p = frequency of AssessedAllele
# N = SNP-specific NrSamples

dat2[, denom_eqtl := sqrt(
  
  2 *
    EAF_eqtl *
    (1 - EAF_eqtl) *
    (
      NrSamples +
        Zscore^2
    )
)]


dat2[, beta_eqtl :=
       Zscore /
       denom_eqtl
]

dat2[, se_eqtl :=
       1 /
       denom_eqtl
]


############################
# 18. Validate reconstructed eQTL effects
############################

z_reconstructed_error <-
  dat2$beta_eqtl /
  dat2$se_eqtl -
  dat2$Zscore


cat(
  "\nMaximum |beta/SE - Z|:",
  max(abs(z_reconstructed_error), na.rm = TRUE),
  "\n"
)


stopifnot(
  all(is.finite(dat2$beta_eqtl)),
  all(is.finite(dat2$se_eqtl)),
  all(dat2$se_eqtl > 0),
  max(
    abs(z_reconstructed_error),
    na.rm = TRUE
  ) < 1e-10
)


############################
# 19. CRC QC
############################

stopifnot(
  all(is.finite(dat2$beta_crc_h)),
  all(is.finite(dat2$SE)),
  all(dat2$SE > 0),
  all(is.finite(dat2$MAF_crc)),
  all(dat2$MAF_crc > 0),
  all(dat2$MAF_crc <= 0.5)
)


############################
# 20. CRC study parameters
############################

N_crc <- 470002
N_case <- 6581
N_control <- 463421

stopifnot(
  N_case + N_control == N_crc
)

s_crc <- N_case / N_crc


############################
# 21. General coloc function
############################

run_coloc <- function(
    gene,
    input_data,
    p1 = 1e-4,
    p2 = 1e-4,
    p12 = 1e-5
) {
  
  d <- copy(
    input_data[
      GeneSymbol == gene
    ]
  )
  
  
  # Final completeness/QC
  d <- d[
    is.finite(beta_eqtl) &
      is.finite(se_eqtl) &
      is.finite(beta_crc_h) &
      is.finite(SE) &
      is.finite(MAF_eqtl) &
      is.finite(MAF_crc) &
      se_eqtl > 0 &
      SE > 0 &
      MAF_eqtl > 0 &
      MAF_eqtl <= 0.5 &
      MAF_crc > 0 &
      MAF_crc <= 0.5
  ]
  
  
  if (nrow(d) < 100) {
    stop(
      paste(
        gene,
        "has fewer than 100 usable SNPs."
      )
    )
  }
  
  
  if (anyDuplicated(d$SNP) != 0) {
    stop(
      paste(
        "Duplicated SNPs remain for",
        gene
      )
    )
  }
  
  
  ############################
  # Trait 1: eQTL
  ############################
  
  D1 <- list(
    beta = d$beta_eqtl,
    varbeta = d$se_eqtl^2,
    snp = d$SNP,
    position = d$SNPPos,
    MAF = d$MAF_eqtl,
    type = "quant",
    sdY = 1
  )
  
  
  ############################
  # Trait 2: CRC
  ############################
  
  D2 <- list(
    beta = d$beta_crc_h,
    varbeta = d$SE^2,
    snp = d$SNP,
    position = d$POS,
    MAF = d$MAF_crc,
    type = "cc",
    N = N_crc,
    s = s_crc
  )
  
  
  # Formal coloc input checks
  coloc::check_dataset(D1)
  coloc::check_dataset(D2)
  
  
  # Confirm exactly the same SNPs/order
  stopifnot(
    identical(
      D1$snp,
      D2$snp
    )
  )
  
  
  result <- coloc::coloc.abf(
    dataset1 = D1,
    dataset2 = D2,
    p1 = p1,
    p2 = p2,
    p12 = p12
  )
  
  
  return(
    list(
      gene = gene,
      nsnp = nrow(d),
      result = result,
      D1 = D1,
      D2 = D2,
      data = d
    )
  )
}


############################
# 22. Summary extraction function
############################

extract_coloc_summary <- function(
    result_list,
    analysis_name
) {
  
  rbindlist(
    lapply(
      names(result_list),
      function(g) {
        
        s <-
          result_list[[g]]$result$summary
        
        
        data.table(
          Analysis = analysis_name,
          Gene = g,
          NSNP = as.integer(
            s["nsnps"]
          ),
          PP.H0 = as.numeric(
            s["PP.H0.abf"]
          ),
          PP.H1 = as.numeric(
            s["PP.H1.abf"]
          ),
          PP.H2 = as.numeric(
            s["PP.H2.abf"]
          ),
          PP.H3 = as.numeric(
            s["PP.H3.abf"]
          ),
          PP.H4 = as.numeric(
            s["PP.H4.abf"]
          )
        )
      }
    )
  )
}


############################
# 23. Run MAIN analysis
############################

genes <- c(
  "AKR1C2",
  "MAPK3",
  "MAP3K11",
  "FADS1",
  "BID"
)


main_results <- lapply(
  genes,
  function(g) {
    
    cat(
      "\nRunning MAIN analysis:",
      g,
      "\n"
    )
    
    run_coloc(
      gene = g,
      input_data = dat2
    )
  }
)

names(main_results) <- genes


############################
# 24. MAIN summary
############################

main_summary <-
  extract_coloc_summary(
    main_results,
    analysis_name = "Main"
  )


main_summary[, PP4_conditional :=
               PP.H4 /
               (
                 PP.H3 +
                   PP.H4
               )
]


main_summary[, Interpretation :=
               fifelse(
                 PP.H4 >= 0.80,
                 
                 "Strong evidence for colocalization",
                 
                 fifelse(
                   PP.H3 >= 0.80,
                   
                   "Evidence for distinct causal signals",
                   
                   "Inconclusive"
                 )
               )
]


cat(
  "\n==============================\n",
  "MAIN COLOCALIZATION RESULTS\n",
  "==============================\n"
)

print(main_summary)


############################
# 25. Export MAIN result
############################

fwrite(
  main_summary,
  "Colocalization_MAIN_5genes.csv"
)


saveRDS(
  main_results,
  "Colocalization_MAIN_5genes.rds"
)


############################
# 26. SNP-level posterior results
############################

main_snp_results <- rbindlist(
  lapply(
    genes,
    function(g) {
      
      x <- as.data.table(
        main_results[[g]]$result$results
      )
      
      x[, Gene := g]
      
      x
    }
  ),
  fill = TRUE
)


fwrite(
  main_snp_results,
  "Colocalization_MAIN_SNP_level.csv"
)


############################
# 27. Top shared-variant candidates
############################

main_top_snps <- rbindlist(
  lapply(
    genes,
    function(g) {
      
      x <- as.data.table(
        main_results[[g]]$result$results
      )
      
      x <- x[
        order(-SNP.PP.H4)
      ]
      
      n_top <- min(
        10,
        nrow(x)
      )
      
      y <- x[
        seq_len(n_top),
        .(
          SNP = snp,
          SNP.PP.H4
        )
      ]
      
      y[, Gene := g]
      
      y[, Position :=
          main_results[[g]]$data$SNPPos[
            match(
              SNP,
              main_results[[g]]$data$SNP
            )
          ]
      ]
      
      setcolorder(
        y,
        c(
          "Gene",
          "SNP",
          "Position",
          "SNP.PP.H4"
        )
      )
      
      y
    }
  )
)


fwrite(
  main_top_snps,
  "Colocalization_MAIN_top10_SNPs.csv"
)


############################
# 28. Save cleaned dataset
############################

fwrite(
  dat2,
  "Colocalization_cleaned_master_data.csv"
)

saveRDS(
  dat2,
  "Colocalization_cleaned_master_data.rds"
)

############################################################
# END OF MAIN ANALYSIS
############################################################
############################################################
# SENSITIVITY ANALYSIS
#
# Exclude SNPs with:
# |EAF_eQTL - EAF_CRC| > 0.20
############################################################


############################
# 1. Create sensitivity dataset
############################

dat_sensitivity <- copy(
  dat2[
    AF_diff <= 0.20
  ]
)


cat(
  "\nMain dataset SNPs:",
  nrow(dat2),
  "\n"
)

cat(
  "Sensitivity dataset SNPs:",
  nrow(dat_sensitivity),
  "\n"
)

cat(
  "Excluded SNPs:",
  nrow(dat2) -
    nrow(dat_sensitivity),
  "\n"
)


cat("\nSensitivity SNP counts by gene:\n")

print(
  dat_sensitivity[
    ,
    .N,
    by = GeneSymbol
  ][order(GeneSymbol)]
)


############################
# 2. Confirm no duplicate SNPs
############################

stopifnot(
  sum(
    duplicated(
      dat_sensitivity$SNP
    )
  ) == 0
)


############################
# 3. Run sensitivity coloc
############################

sensitivity_results <- lapply(
  genes,
  function(g) {
    
    cat(
      "\nRunning SENSITIVITY analysis:",
      g,
      "\n"
    )
    
    run_coloc(
      gene = g,
      input_data = dat_sensitivity
    )
  }
)


names(
  sensitivity_results
) <- genes


############################
# 4. Sensitivity summary
############################

sensitivity_summary <-
  extract_coloc_summary(
    sensitivity_results,
    analysis_name =
      "Sensitivity_AFdiff_le_0.20"
  )


sensitivity_summary[
  ,
  PP4_conditional :=
    PP.H4 /
    (
      PP.H3 +
        PP.H4
    )
]


sensitivity_summary[
  ,
  Interpretation :=
    fifelse(
      PP.H4 >= 0.80,
      
      "Strong evidence for colocalization",
      
      fifelse(
        PP.H3 >= 0.80,
        
        "Evidence for distinct causal signals",
        
        "Inconclusive"
      )
    )
]


cat(
  "\n==============================\n",
  "SENSITIVITY RESULTS\n",
  "==============================\n"
)

print(
  sensitivity_summary
)


############################
# 5. Export sensitivity result
############################

fwrite(
  sensitivity_summary,
  "Colocalization_SENSITIVITY_AFdiff020_5genes.csv"
)


saveRDS(
  sensitivity_results,
  "Colocalization_SENSITIVITY_AFdiff020_5genes.rds"
)


############################
# 6. Main vs sensitivity comparison
############################

comparison <- merge(
  main_summary[
    ,
    .(
      Gene,
      NSNP_main = NSNP,
      H0_main = PP.H0,
      H1_main = PP.H1,
      H2_main = PP.H2,
      H3_main = PP.H3,
      H4_main = PP.H4
    )
  ],
  
  sensitivity_summary[
    ,
    .(
      Gene,
      NSNP_sensitivity = NSNP,
      H0_sensitivity = PP.H0,
      H1_sensitivity = PP.H1,
      H2_sensitivity = PP.H2,
      H3_sensitivity = PP.H3,
      H4_sensitivity = PP.H4
    )
  ],
  
  by = "Gene"
)


comparison[
  ,
  Delta_H4 :=
    H4_sensitivity -
    H4_main
]


comparison[
  ,
  Delta_H3 :=
    H3_sensitivity -
    H3_main
]


cat(
  "\n==============================\n",
  "MAIN VS SENSITIVITY\n",
  "==============================\n"
)

print(
  comparison
)


############################
# 7. Export comparison
############################

fwrite(
  comparison,
  "Colocalization_MAIN_vs_SENSITIVITY.csv"
)


############################################################
# END OF SENSITIVITY ANALYSIS
############################################################