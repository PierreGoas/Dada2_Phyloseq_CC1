Dada2 tutorial
================

  - [Tutoriel Dada2](#tutoriel-dada2)
      - [Importation des librairies et préparation des
        données](#importation-des-librairies-et-préparation-des-données)
      - [Inspection des profils de qualité des
        reads](#inspection-des-profils-de-qualité-des-reads)
      - [Filtration et tri des
        donnnées](#filtration-et-tri-des-donnnées)
      - [Calcul des erreurs et visualisation de ces
        derniers](#calcul-des-erreurs-et-visualisation-de-ces-derniers)
      - [Application de Dada2 aux
        données](#application-de-dada2-aux-données)
      - [Alignement des séquences Forward et
        Reverse](#alignement-des-séquences-forward-et-reverse)
      - [Table d’observation](#table-dobservation)
      - [Elimination des séquences
        chimériques](#elimination-des-séquences-chimériques)
      - [Résumé des opérations effectuées
        précédement](#résumé-des-opérations-effectuées-précédement)
      - [Téléchargement des données
        Silva](#téléchargement-des-données-silva)
      - [Téléchargement des données Silva permmetant l’assignement
        d’espèces aux des séquences
        analysées](#téléchargement-des-données-silva-permmetant-lassignement-despèces-aux-des-séquences-analysées)
      - [Evaluation de la précisison](#evaluation-de-la-précisison)
  - [Premier pas sur Phyloseq (fin tuto
    Dada2)](#premier-pas-sur-phyloseq-fin-tuto-dada2)
      - [Importation des librairies et préparation des
        données](#importation-des-librairies-et-préparation-des-données-1)
      - [Visualisation de l’alpha
        diversité](#visualisation-de-lalpha-diversité)
      - [Ordination](#ordination)
      - [Représention des échantillons en “histogramme” (bar
        plot)](#représention-des-échantillons-en-histogramme-bar-plot)

# Tutoriel Dada2

## Importation des librairies et préparation des données

``` r
library("Rcpp")
library("dada2")
```

On importe les galeries nécessaires pour la réalisation de l’excercice

``` r
path <- "~/MiSeq_SOP"
list.files(path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

Chargement des données fastq.

``` r
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

On sépare les forward et les reverse. Les forward seront associés à fnFS
et les reverse à fnRs. Association des fnFs et fnRS entre eux
(correspondance entre les données)

## Inspection des profils de qualité des reads

``` r
plotQualityProfile(fnFs[1:4])
```

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
Affichage des graphes de qualité des données fnFs. En gris : heat map de
la fréquence de chaque score de qualité de chaque position de base En
vert : score de qualité moyen de chaque base En rouge : proportion
graduée des reads jusqu’à une certaine position (ici, elle est plate
car les reads faites sur Illumina sont de même longueur donc ce n’est
pas très utile pour ce cas présent) Résultat : bon. On va couper les
derniers nucléotides car on observe une chute sur la partie droite du
graphe. Cette coupure de nucléotide permettra de donner des résultats
plus fiable pour la suite. Ici, on va tronquer les nucléotides à partir
de la position 240 (on retire les 10 derniers nucléotides).

``` r
plotQualityProfile(fnRs[1:4])
```

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
Affichage du graphe de qualité des données fnRs Même codes couleurs que
fnFs Résultat : les reads reverse sont de moins bonne qualité en
particulier vers la fin. Ce n’est pas grave mais les nucléotides à
partir de la position 160 seront supprimés (soit 90 nucléotides)

## Filtration et tri des donnnées

``` r
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

Les données fnFs et fnRs vont être filtrées. Tout d’abord, elles vont
être placées sous un nouveau nom (ici, filtFs et filtRs)

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

Les données sont filtrés, ie que les derniers nucléotides vont être
retirés (fnFs : à partir de 240; fnRd : à partir de 160). Seulement 2
erreurs peuvent être acceptées lors de la filtration

## Calcul des erreurs et visualisation de ces derniers

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

Estimation du taux d’erreur de filtFs

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

Estimation du taux d’erreur de filtRs

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->
Représentation des fréquences d’erreurs estimées. Points gris : taux
d’erreurs observées pour chacun des scores de qualité consensus. Ligne
noire : taux d’erreurs estimé après que l’algorithme ait réuni toutes
les informations liée aux taux d’erreurs estimés Ligne rouge : taux
d’erreurs attendu selon le Q-score Corrélation entre les taux
d’erreurs estimés et observés, diminution quand le score de qualité
augmente -\> on peut poursuivre.

## Application de Dada2 aux données

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

Application de l’algorithme d’inférence aux données filtrées et ajustées
des reads forward

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

Application de l’algorithme d’inférence aux données filtrées et ajustées
des reads reverse

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

Dada2 a déduit 128 séquences variantes à partir des 1979 séquences
uniques.

## Alignement des séquences Forward et Reverse

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 6551 paired-reads (in 106 unique pairings) successfully merged out of 6907 (in 199 pairings) input.

    ## 5025 paired-reads (in 100 unique pairings) successfully merged out of 5188 (in 156 pairings) input.

    ## 4973 paired-reads (in 80 unique pairings) successfully merged out of 5268 (in 166 pairings) input.

    ## 2595 paired-reads (in 52 unique pairings) successfully merged out of 2756 (in 109 pairings) input.

    ## 2553 paired-reads (in 60 unique pairings) successfully merged out of 2785 (in 119 pairings) input.

    ## 3622 paired-reads (in 53 unique pairings) successfully merged out of 4103 (in 157 pairings) input.

    ## 6079 paired-reads (in 81 unique pairings) successfully merged out of 6515 (in 198 pairings) input.

    ## 3961 paired-reads (in 90 unique pairings) successfully merged out of 4384 (in 188 pairings) input.

    ## 14231 paired-reads (in 143 unique pairings) successfully merged out of 15358 (in 351 pairings) input.

    ## 10526 paired-reads (in 120 unique pairings) successfully merged out of 11166 (in 279 pairings) input.

    ## 11156 paired-reads (in 137 unique pairings) successfully merged out of 11799 (in 298 pairings) input.

    ## 4329 paired-reads (in 84 unique pairings) successfully merged out of 4788 (in 180 pairings) input.

    ## 17431 paired-reads (in 153 unique pairings) successfully merged out of 17812 (in 272 pairings) input.

    ## 5850 paired-reads (in 81 unique pairings) successfully merged out of 6095 (in 159 pairings) input.

    ## 3716 paired-reads (in 86 unique pairings) successfully merged out of 3894 (in 147 pairings) input.

    ## 6865 paired-reads (in 99 unique pairings) successfully merged out of 7193 (in 187 pairings) input.

    ## 4430 paired-reads (in 67 unique pairings) successfully merged out of 4605 (in 127 pairings) input.

    ## 4574 paired-reads (in 100 unique pairings) successfully merged out of 4736 (in 172 pairings) input.

    ## 6094 paired-reads (in 109 unique pairings) successfully merged out of 6314 (in 172 pairings) input.

    ## 4269 paired-reads (in 20 unique pairings) successfully merged out of 4281 (in 28 pairings) input.

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                       sequence
    ## 1 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACAGG
    ## 2 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAGG
    ## 3 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACCCCGGCCTGCCGTTGAAACTGGCGGCCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ## 4 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTTTTAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCCGGCCGTTGAAACTGGGGGCCTTGAGTGGGCGAGAAGAAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCCTTCCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCGAACAGG
    ## 5 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGACTCTCAAGTCAGCGGTCAAATCGCGGGGCTCAACCCCGTTCCGCCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAAGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTCATGCACGAAAGCGTGGGTATCGAACAGG
    ## 6 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1       579       1       1    148         0      0      1   TRUE
    ## 2       470       2       2    148         0      0      2   TRUE
    ## 3       449       3       4    148         0      0      1   TRUE
    ## 4       430       4       3    148         0      0      2   TRUE
    ## 5       345       5       6    148         0      0      1   TRUE
    ## 6       282       6       5    148         0      0      2   TRUE

Alignement des séquences pour former des contigues. L’alignement se fait
avec les séquences forward et reverse se correspondant sur au moins 12
bases.

## Table d’observation

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]  20 293

Construction d’ASV (table d’OTU). Ici, sur les 20 échantillons, 293 ASV
ont été crée

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  88 196   6   2

Visualisation de la distribution des longueurs de séquences. Ici, les
différentes séquences obtenues sont d’une longueur située dans la région
V4 de l’ARN16S

## Elimination des séquences chimériques

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 61 bimeras out of 293 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]  20 232

Dada2 a corrigé les erreurs de substitutions et d’index. Mais il reste
les séquences chimériques (séquence contenant de l’information de
plusiseurs autres séquences) à retirer. L’identification se fait par
construction de deux séquences (une gauche et une droite), provenant de
séquences dites parents retrouvées abondament dans notre échantillon.
Ici, 20232 séquences chimériques ont été identifiées

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.964263

Ces séquences chimériques représentent 3.5% de notre jeu de donnée

## Résumé des opérations effectuées précédement

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##        input filtered denoisedF denoisedR merged nonchim
    ## F3D0    7793     7113      6996      6978   6551    6539
    ## F3D1    5869     5299      5227      5239   5025    5014
    ## F3D141  5958     5463      5339      5351   4973    4850
    ## F3D142  3183     2914      2799      2833   2595    2521
    ## F3D143  3178     2941      2822      2868   2553    2519
    ## F3D144  4827     4312      4146      4224   3622    3483

Récapitulatif des reads qui ont traversé les différentes étapes de
filtration et de tri.

## Téléchargement des données Silva

``` bash
wget https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz
```

    ## --2020-12-04 14:26:19--  https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 137973851 (132M) [application/octet-stream]
    ## Saving to: ‘silva_nr99_v138_train_set.fa.gz.3’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 14.8M 9s
    ##     50K .......... .......... .......... .......... ..........  0% 15.3M 9s
    ##    100K .......... .......... .......... .......... ..........  0% 10.2M 10s
    ##    150K .......... .......... .......... .......... ..........  0% 26.7M 9s
    ##    200K .......... .......... .......... .......... ..........  0% 19.6M 8s
    ##    250K .......... .......... .......... .......... ..........  0% 47.2M 7s
    ##    300K .......... .......... .......... .......... ..........  0% 14.0M 8s
    ##    350K .......... .......... .......... .......... ..........  0% 97.8M 7s
    ##    400K .......... .......... .......... .......... ..........  0% 13.7M 7s
    ##    450K .......... .......... .......... .......... ..........  0% 46.8M 7s
    ##    500K .......... .......... .......... .......... ..........  0% 16.0M 7s
    ##    550K .......... .......... .......... .......... ..........  0% 15.2M 7s
    ##    600K .......... .......... .......... .......... ..........  0% 60.4M 7s
    ##    650K .......... .......... .......... .......... ..........  0% 18.9M 7s
    ##    700K .......... .......... .......... .......... ..........  0% 33.2M 7s
    ##    750K .......... .......... .......... .......... ..........  0% 13.9M 7s
    ##    800K .......... .......... .......... .......... ..........  0% 56.6M 6s
    ##    850K .......... .......... .......... .......... ..........  0% 21.4M 6s
    ##    900K .......... .......... .......... .......... ..........  0% 32.1M 6s
    ##    950K .......... .......... .......... .......... ..........  0% 7.16M 7s
    ##   1000K .......... .......... .......... .......... ..........  0% 46.9M 7s
    ##   1050K .......... .......... .......... .......... ..........  0% 24.0M 7s
    ##   1100K .......... .......... .......... .......... ..........  0% 35.2M 6s
    ##   1150K .......... .......... .......... .......... ..........  0% 26.5M 6s
    ##   1200K .......... .......... .......... .......... ..........  0% 21.7M 6s
    ##   1250K .......... .......... .......... .......... ..........  0% 21.9M 6s
    ##   1300K .......... .......... .......... .......... ..........  1% 63.6M 6s
    ##   1350K .......... .......... .......... .......... ..........  1% 67.3M 6s
    ##   1400K .......... .......... .......... .......... ..........  1% 21.3M 6s
    ##   1450K .......... .......... .......... .......... ..........  1% 37.0M 6s
    ##   1500K .......... .......... .......... .......... ..........  1% 16.4M 6s
    ##   1550K .......... .......... .......... .......... ..........  1% 72.9M 6s
    ##   1600K .......... .......... .......... .......... ..........  1% 64.4M 6s
    ##   1650K .......... .......... .......... .......... ..........  1% 16.0M 6s
    ##   1700K .......... .......... .......... .......... ..........  1% 59.9M 6s
    ##   1750K .......... .......... .......... .......... ..........  1% 71.6M 6s
    ##   1800K .......... .......... .......... .......... ..........  1% 14.1M 6s
    ##   1850K .......... .......... .......... .......... ..........  1% 95.7M 6s
    ##   1900K .......... .......... .......... .......... ..........  1% 12.4M 6s
    ##   1950K .......... .......... .......... .......... ..........  1%  101M 6s
    ##   2000K .......... .......... .......... .......... ..........  1% 80.0M 6s
    ##   2050K .......... .......... .......... .......... ..........  1% 13.1M 6s
    ##   2100K .......... .......... .......... .......... ..........  1% 67.5M 6s
    ##   2150K .......... .......... .......... .......... ..........  1%  101M 5s
    ##   2200K .......... .......... .......... .......... ..........  1% 15.8M 5s
    ##   2250K .......... .......... .......... .......... ..........  1% 78.3M 5s
    ##   2300K .......... .......... .......... .......... ..........  1% 15.1M 5s
    ##   2350K .......... .......... .......... .......... ..........  1% 62.2M 5s
    ##   2400K .......... .......... .......... .......... ..........  1% 29.3M 5s
    ##   2450K .......... .......... .......... .......... ..........  1% 36.5M 5s
    ##   2500K .......... .......... .......... .......... ..........  1% 68.3M 5s
    ##   2550K .......... .......... .......... .......... ..........  1% 32.7M 5s
    ##   2600K .......... .......... .......... .......... ..........  1% 22.1M 5s
    ##   2650K .......... .......... .......... .......... ..........  2% 74.0M 5s
    ##   2700K .......... .......... .......... .......... ..........  2% 20.2M 5s
    ##   2750K .......... .......... .......... .......... ..........  2% 68.9M 5s
    ##   2800K .......... .......... .......... .......... ..........  2% 62.7M 5s
    ##   2850K .......... .......... .......... .......... ..........  2% 23.2M 5s
    ##   2900K .......... .......... .......... .......... ..........  2% 40.8M 5s
    ##   2950K .......... .......... .......... .......... ..........  2% 71.9M 5s
    ##   3000K .......... .......... .......... .......... ..........  2% 17.4M 5s
    ##   3050K .......... .......... .......... .......... ..........  2% 57.9M 5s
    ##   3100K .......... .......... .......... .......... ..........  2% 18.1M 5s
    ##   3150K .......... .......... .......... .......... ..........  2% 92.7M 5s
    ##   3200K .......... .......... .......... .......... ..........  2% 43.3M 5s
    ##   3250K .......... .......... .......... .......... ..........  2% 17.7M 5s
    ##   3300K .......... .......... .......... .......... ..........  2% 61.4M 5s
    ##   3350K .......... .......... .......... .......... ..........  2% 17.5M 5s
    ##   3400K .......... .......... .......... .......... ..........  2% 66.1M 5s
    ##   3450K .......... .......... .......... .......... ..........  2% 87.1M 5s
    ##   3500K .......... .......... .......... .......... ..........  2% 16.0M 5s
    ##   3550K .......... .......... .......... .......... ..........  2% 57.1M 5s
    ##   3600K .......... .......... .......... .......... ..........  2% 18.3M 5s
    ##   3650K .......... .......... .......... .......... ..........  2% 27.7M 5s
    ##   3700K .......... .......... .......... .......... ..........  2% 78.2M 5s
    ##   3750K .......... .......... .......... .......... ..........  2% 10.5M 5s
    ##   3800K .......... .......... .......... .......... ..........  2% 77.1M 5s
    ##   3850K .......... .......... .......... .......... ..........  2% 93.2M 5s
    ##   3900K .......... .......... .......... .......... ..........  2% 13.2M 5s
    ##   3950K .......... .......... .......... .......... ..........  2% 69.6M 5s
    ##   4000K .......... .......... .......... .......... ..........  3% 77.0M 5s
    ##   4050K .......... .......... .......... .......... ..........  3% 22.7M 5s
    ##   4100K .......... .......... .......... .......... ..........  3% 53.3M 5s
    ##   4150K .......... .......... .......... .......... ..........  3% 67.8M 5s
    ##   4200K .......... .......... .......... .......... ..........  3% 20.5M 5s
    ##   4250K .......... .......... .......... .......... ..........  3% 59.6M 5s
    ##   4300K .......... .......... .......... .......... ..........  3% 23.3M 5s
    ##   4350K .......... .......... .......... .......... ..........  3% 81.8M 5s
    ##   4400K .......... .......... .......... .......... ..........  3% 28.7M 5s
    ##   4450K .......... .......... .......... .......... ..........  3% 26.4M 5s
    ##   4500K .......... .......... .......... .......... ..........  3% 8.84M 5s
    ##   4550K .......... .......... .......... .......... ..........  3%  101M 5s
    ##   4600K .......... .......... .......... .......... ..........  3% 83.6M 5s
    ##   4650K .......... .......... .......... .......... ..........  3% 16.9M 5s
    ##   4700K .......... .......... .......... .......... ..........  3% 62.8M 5s
    ##   4750K .......... .......... .......... .......... ..........  3% 57.9M 5s
    ##   4800K .......... .......... .......... .......... ..........  3% 13.1M 5s
    ##   4850K .......... .......... .......... .......... ..........  3% 41.4M 5s
    ##   4900K .......... .......... .......... .......... ..........  3% 77.5M 5s
    ##   4950K .......... .......... .......... .......... ..........  3% 88.3M 5s
    ##   5000K .......... .......... .......... .......... ..........  3% 26.2M 5s
    ##   5050K .......... .......... .......... .......... ..........  3% 28.5M 5s
    ##   5100K .......... .......... .......... .......... ..........  3% 53.0M 5s
    ##   5150K .......... .......... .......... .......... ..........  3%  103M 5s
    ##   5200K .......... .......... .......... .......... ..........  3% 25.3M 5s
    ##   5250K .......... .......... .......... .......... ..........  3% 69.4M 5s
    ##   5300K .......... .......... .......... .......... ..........  3% 59.6M 5s
    ##   5350K .......... .......... .......... .......... ..........  4% 82.7M 5s
    ##   5400K .......... .......... .......... .......... ..........  4% 24.0M 5s
    ##   5450K .......... .......... .......... .......... ..........  4% 66.0M 4s
    ##   5500K .......... .......... .......... .......... ..........  4% 72.3M 4s
    ##   5550K .......... .......... .......... .......... ..........  4% 74.9M 4s
    ##   5600K .......... .......... .......... .......... ..........  4% 24.3M 4s
    ##   5650K .......... .......... .......... .......... ..........  4% 68.3M 4s
    ##   5700K .......... .......... .......... .......... ..........  4% 62.5M 4s
    ##   5750K .......... .......... .......... .......... ..........  4% 95.7M 4s
    ##   5800K .......... .......... .......... .......... ..........  4% 20.6M 4s
    ##   5850K .......... .......... .......... .......... ..........  4% 21.9M 4s
    ##   5900K .......... .......... .......... .......... ..........  4% 87.1M 4s
    ##   5950K .......... .......... .......... .......... ..........  4% 4.26M 5s
    ##   6000K .......... .......... .......... .......... ..........  4% 71.6M 5s
    ##   6050K .......... .......... .......... .......... ..........  4% 39.5M 5s
    ##   6100K .......... .......... .......... .......... ..........  4% 74.5M 5s
    ##   6150K .......... .......... .......... .......... ..........  4% 33.7M 5s
    ##   6200K .......... .......... .......... .......... ..........  4% 8.41M 5s
    ##   6250K .......... .......... .......... .......... ..........  4% 63.0M 5s
    ##   6300K .......... .......... .......... .......... ..........  4% 70.4M 5s
    ##   6350K .......... .......... .......... .......... ..........  4% 69.8M 5s
    ##   6400K .......... .......... .......... .......... ..........  4% 84.2M 4s
    ##   6450K .......... .......... .......... .......... ..........  4% 26.4M 4s
    ##   6500K .......... .......... .......... .......... ..........  4% 57.4M 4s
    ##   6550K .......... .......... .......... .......... ..........  4% 29.8M 4s
    ##   6600K .......... .......... .......... .......... ..........  4% 67.5M 4s
    ##   6650K .......... .......... .......... .......... ..........  4% 79.5M 4s
    ##   6700K .......... .......... .......... .......... ..........  5% 32.3M 4s
    ##   6750K .......... .......... .......... .......... ..........  5% 52.8M 4s
    ##   6800K .......... .......... .......... .......... ..........  5% 67.9M 4s
    ##   6850K .......... .......... .......... .......... ..........  5% 54.0M 4s
    ##   6900K .......... .......... .......... .......... ..........  5% 64.6M 4s
    ##   6950K .......... .......... .......... .......... ..........  5% 58.0M 4s
    ##   7000K .......... .......... .......... .......... ..........  5% 29.3M 4s
    ##   7050K .......... .......... .......... .......... ..........  5% 79.2M 4s
    ##   7100K .......... .......... .......... .......... ..........  5% 69.1M 4s
    ##   7150K .......... .......... .......... .......... ..........  5% 5.88M 4s
    ##   7200K .......... .......... .......... .......... ..........  5% 34.6M 4s
    ##   7250K .......... .......... .......... .......... ..........  5% 37.9M 4s
    ##   7300K .......... .......... .......... .......... ..........  5% 39.4M 4s
    ##   7350K .......... .......... .......... .......... ..........  5% 58.0M 4s
    ##   7400K .......... .......... .......... .......... ..........  5% 44.5M 4s
    ##   7450K .......... .......... .......... .......... ..........  5% 46.1M 4s
    ##   7500K .......... .......... .......... .......... ..........  5% 35.4M 4s
    ##   7550K .......... .......... .......... .......... ..........  5% 49.9M 4s
    ##   7600K .......... .......... .......... .......... ..........  5% 36.1M 4s
    ##   7650K .......... .......... .......... .......... ..........  5% 46.6M 4s
    ##   7700K .......... .......... .......... .......... ..........  5% 36.5M 4s
    ##   7750K .......... .......... .......... .......... ..........  5% 7.52M 4s
    ##   7800K .......... .......... .......... .......... ..........  5% 48.7M 4s
    ##   7850K .......... .......... .......... .......... ..........  5% 57.3M 4s
    ##   7900K .......... .......... .......... .......... ..........  5% 51.5M 4s
    ##   7950K .......... .......... .......... .......... ..........  5% 49.9M 4s
    ##   8000K .......... .......... .......... .......... ..........  5% 43.3M 4s
    ##   8050K .......... .......... .......... .......... ..........  6% 44.7M 4s
    ##   8100K .......... .......... .......... .......... ..........  6% 42.1M 4s
    ##   8150K .......... .......... .......... .......... ..........  6% 43.3M 4s
    ##   8200K .......... .......... .......... .......... ..........  6% 43.2M 4s
    ##   8250K .......... .......... .......... .......... ..........  6% 51.9M 4s
    ##   8300K .......... .......... .......... .......... ..........  6% 45.2M 4s
    ##   8350K .......... .......... .......... .......... ..........  6% 60.0M 4s
    ##   8400K .......... .......... .......... .......... ..........  6% 39.9M 4s
    ##   8450K .......... .......... .......... .......... ..........  6% 55.8M 4s
    ##   8500K .......... .......... .......... .......... ..........  6% 49.5M 4s
    ##   8550K .......... .......... .......... .......... ..........  6% 57.0M 4s
    ##   8600K .......... .......... .......... .......... ..........  6% 39.5M 4s
    ##   8650K .......... .......... .......... .......... ..........  6% 42.4M 4s
    ##   8700K .......... .......... .......... .......... ..........  6% 38.6M 4s
    ##   8750K .......... .......... .......... .......... ..........  6% 45.9M 4s
    ##   8800K .......... .......... .......... .......... ..........  6% 36.6M 4s
    ##   8850K .......... .......... .......... .......... ..........  6% 51.8M 4s
    ##   8900K .......... .......... .......... .......... ..........  6% 37.6M 4s
    ##   8950K .......... .......... .......... .......... ..........  6% 46.2M 4s
    ##   9000K .......... .......... .......... .......... ..........  6% 41.6M 4s
    ##   9050K .......... .......... .......... .......... ..........  6% 59.6M 4s
    ##   9100K .......... .......... .......... .......... ..........  6% 46.5M 4s
    ##   9150K .......... .......... .......... .......... ..........  6% 65.5M 4s
    ##   9200K .......... .......... .......... .......... ..........  6% 48.7M 4s
    ##   9250K .......... .......... .......... .......... ..........  6% 50.8M 4s
    ##   9300K .......... .......... .......... .......... ..........  6% 44.4M 4s
    ##   9350K .......... .......... .......... .......... ..........  6% 48.8M 4s
    ##   9400K .......... .......... .......... .......... ..........  7% 48.2M 4s
    ##   9450K .......... .......... .......... .......... ..........  7% 59.1M 4s
    ##   9500K .......... .......... .......... .......... ..........  7% 48.7M 4s
    ##   9550K .......... .......... .......... .......... ..........  7% 55.3M 4s
    ##   9600K .......... .......... .......... .......... ..........  7% 51.9M 4s
    ##   9650K .......... .......... .......... .......... ..........  7% 58.4M 4s
    ##   9700K .......... .......... .......... .......... ..........  7% 52.2M 4s
    ##   9750K .......... .......... .......... .......... ..........  7% 64.8M 4s
    ##   9800K .......... .......... .......... .......... ..........  7% 53.3M 4s
    ##   9850K .......... .......... .......... .......... ..........  7% 47.9M 4s
    ##   9900K .......... .......... .......... .......... ..........  7% 49.3M 4s
    ##   9950K .......... .......... .......... .......... ..........  7% 56.1M 4s
    ##  10000K .......... .......... .......... .......... ..........  7% 52.5M 4s
    ##  10050K .......... .......... .......... .......... ..........  7% 61.1M 4s
    ##  10100K .......... .......... .......... .......... ..........  7% 48.5M 4s
    ##  10150K .......... .......... .......... .......... ..........  7% 47.5M 4s
    ##  10200K .......... .......... .......... .......... ..........  7% 58.3M 4s
    ##  10250K .......... .......... .......... .......... ..........  7% 81.8M 4s
    ##  10300K .......... .......... .......... .......... ..........  7% 68.7M 4s
    ##  10350K .......... .......... .......... .......... ..........  7% 69.5M 4s
    ##  10400K .......... .......... .......... .......... ..........  7% 59.2M 4s
    ##  10450K .......... .......... .......... .......... ..........  7% 59.4M 4s
    ##  10500K .......... .......... .......... .......... ..........  7% 68.0M 4s
    ##  10550K .......... .......... .......... .......... ..........  7% 83.1M 4s
    ##  10600K .......... .......... .......... .......... ..........  7% 63.4M 4s
    ##  10650K .......... .......... .......... .......... ..........  7% 70.2M 4s
    ##  10700K .......... .......... .......... .......... ..........  7% 61.1M 4s
    ##  10750K .......... .......... .......... .......... ..........  8% 58.2M 4s
    ##  10800K .......... .......... .......... .......... ..........  8% 53.4M 4s
    ##  10850K .......... .......... .......... .......... ..........  8% 52.8M 4s
    ##  10900K .......... .......... .......... .......... ..........  8% 63.5M 4s
    ##  10950K .......... .......... .......... .......... ..........  8% 66.4M 4s
    ##  11000K .......... .......... .......... .......... ..........  8% 71.6M 4s
    ##  11050K .......... .......... .......... .......... ..........  8% 91.7M 4s
    ##  11100K .......... .......... .......... .......... ..........  8% 62.9M 4s
    ##  11150K .......... .......... .......... .......... ..........  8% 54.5M 4s
    ##  11200K .......... .......... .......... .......... ..........  8% 67.4M 4s
    ##  11250K .......... .......... .......... .......... ..........  8% 84.5M 4s
    ##  11300K .......... .......... .......... .......... ..........  8% 67.7M 4s
    ##  11350K .......... .......... .......... .......... ..........  8% 80.6M 4s
    ##  11400K .......... .......... .......... .......... ..........  8% 65.7M 4s
    ##  11450K .......... .......... .......... .......... ..........  8% 47.8M 4s
    ##  11500K .......... .......... .......... .......... ..........  8% 63.4M 4s
    ##  11550K .......... .......... .......... .......... ..........  8% 65.8M 4s
    ##  11600K .......... .......... .......... .......... ..........  8% 59.0M 4s
    ##  11650K .......... .......... .......... .......... ..........  8% 67.2M 4s
    ##  11700K .......... .......... .......... .......... ..........  8% 60.6M 4s
    ##  11750K .......... .......... .......... .......... ..........  8% 59.1M 4s
    ##  11800K .......... .......... .......... .......... ..........  8% 68.3M 4s
    ##  11850K .......... .......... .......... .......... ..........  8% 70.0M 4s
    ##  11900K .......... .......... .......... .......... ..........  8% 49.2M 4s
    ##  11950K .......... .......... .......... .......... ..........  8% 67.2M 4s
    ##  12000K .......... .......... .......... .......... ..........  8% 40.7M 4s
    ##  12050K .......... .......... .......... .......... ..........  8% 67.7M 4s
    ##  12100K .......... .......... .......... .......... ..........  9% 73.1M 3s
    ##  12150K .......... .......... .......... .......... ..........  9% 83.4M 3s
    ##  12200K .......... .......... .......... .......... ..........  9% 60.4M 3s
    ##  12250K .......... .......... .......... .......... ..........  9% 47.8M 3s
    ##  12300K .......... .......... .......... .......... ..........  9% 34.1M 3s
    ##  12350K .......... .......... .......... .......... ..........  9% 36.9M 3s
    ##  12400K .......... .......... .......... .......... ..........  9% 38.4M 3s
    ##  12450K .......... .......... .......... .......... ..........  9% 19.3M 3s
    ##  12500K .......... .......... .......... .......... ..........  9% 39.8M 3s
    ##  12550K .......... .......... .......... .......... ..........  9% 39.7M 3s
    ##  12600K .......... .......... .......... .......... ..........  9% 27.1M 3s
    ##  12650K .......... .......... .......... .......... ..........  9% 44.9M 3s
    ##  12700K .......... .......... .......... .......... ..........  9% 43.9M 3s
    ##  12750K .......... .......... .......... .......... ..........  9% 35.2M 3s
    ##  12800K .......... .......... .......... .......... ..........  9% 33.6M 3s
    ##  12850K .......... .......... .......... .......... ..........  9% 38.9M 3s
    ##  12900K .......... .......... .......... .......... ..........  9% 32.3M 3s
    ##  12950K .......... .......... .......... .......... ..........  9% 36.8M 3s
    ##  13000K .......... .......... .......... .......... ..........  9% 39.1M 3s
    ##  13050K .......... .......... .......... .......... ..........  9% 41.0M 3s
    ##  13100K .......... .......... .......... .......... ..........  9% 35.7M 3s
    ##  13150K .......... .......... .......... .......... ..........  9% 51.6M 3s
    ##  13200K .......... .......... .......... .......... ..........  9% 47.6M 3s
    ##  13250K .......... .......... .......... .......... ..........  9% 47.4M 3s
    ##  13300K .......... .......... .......... .......... ..........  9% 35.9M 3s
    ##  13350K .......... .......... .......... .......... ..........  9% 45.0M 3s
    ##  13400K .......... .......... .......... .......... ..........  9% 39.3M 3s
    ##  13450K .......... .......... .......... .......... .......... 10% 46.2M 3s
    ##  13500K .......... .......... .......... .......... .......... 10% 39.8M 3s
    ##  13550K .......... .......... .......... .......... .......... 10% 52.2M 3s
    ##  13600K .......... .......... .......... .......... .......... 10% 47.3M 3s
    ##  13650K .......... .......... .......... .......... .......... 10% 52.3M 3s
    ##  13700K .......... .......... .......... .......... .......... 10% 55.9M 3s
    ##  13750K .......... .......... .......... .......... .......... 10% 51.0M 3s
    ##  13800K .......... .......... .......... .......... .......... 10% 52.9M 3s
    ##  13850K .......... .......... .......... .......... .......... 10% 47.9M 3s
    ##  13900K .......... .......... .......... .......... .......... 10% 41.4M 3s
    ##  13950K .......... .......... .......... .......... .......... 10% 38.3M 3s
    ##  14000K .......... .......... .......... .......... .......... 10% 41.4M 3s
    ##  14050K .......... .......... .......... .......... .......... 10% 45.9M 3s
    ##  14100K .......... .......... .......... .......... .......... 10% 49.7M 3s
    ##  14150K .......... .......... .......... .......... .......... 10% 52.3M 3s
    ##  14200K .......... .......... .......... .......... .......... 10% 53.8M 3s
    ##  14250K .......... .......... .......... .......... .......... 10% 56.9M 3s
    ##  14300K .......... .......... .......... .......... .......... 10% 59.4M 3s
    ##  14350K .......... .......... .......... .......... .......... 10% 67.5M 3s
    ##  14400K .......... .......... .......... .......... .......... 10% 55.4M 3s
    ##  14450K .......... .......... .......... .......... .......... 10% 49.8M 3s
    ##  14500K .......... .......... .......... .......... .......... 10% 58.8M 3s
    ##  14550K .......... .......... .......... .......... .......... 10% 56.8M 3s
    ##  14600K .......... .......... .......... .......... .......... 10% 47.3M 3s
    ##  14650K .......... .......... .......... .......... .......... 10% 51.4M 3s
    ##  14700K .......... .......... .......... .......... .......... 10% 55.5M 3s
    ##  14750K .......... .......... .......... .......... .......... 10% 58.6M 3s
    ##  14800K .......... .......... .......... .......... .......... 11% 58.6M 3s
    ##  14850K .......... .......... .......... .......... .......... 11% 67.6M 3s
    ##  14900K .......... .......... .......... .......... .......... 11% 51.2M 3s
    ##  14950K .......... .......... .......... .......... .......... 11% 57.2M 3s
    ##  15000K .......... .......... .......... .......... .......... 11% 54.7M 3s
    ##  15050K .......... .......... .......... .......... .......... 11% 55.6M 3s
    ##  15100K .......... .......... .......... .......... .......... 11% 50.9M 3s
    ##  15150K .......... .......... .......... .......... .......... 11% 50.9M 3s
    ##  15200K .......... .......... .......... .......... .......... 11% 43.9M 3s
    ##  15250K .......... .......... .......... .......... .......... 11% 46.1M 3s
    ##  15300K .......... .......... .......... .......... .......... 11% 49.0M 3s
    ##  15350K .......... .......... .......... .......... .......... 11% 62.6M 3s
    ##  15400K .......... .......... .......... .......... .......... 11% 52.0M 3s
    ##  15450K .......... .......... .......... .......... .......... 11% 58.8M 3s
    ##  15500K .......... .......... .......... .......... .......... 11% 55.7M 3s
    ##  15550K .......... .......... .......... .......... .......... 11% 53.2M 3s
    ##  15600K .......... .......... .......... .......... .......... 11% 51.7M 3s
    ##  15650K .......... .......... .......... .......... .......... 11% 53.8M 3s
    ##  15700K .......... .......... .......... .......... .......... 11% 51.1M 3s
    ##  15750K .......... .......... .......... .......... .......... 11% 48.9M 3s
    ##  15800K .......... .......... .......... .......... .......... 11% 47.1M 3s
    ##  15850K .......... .......... .......... .......... .......... 11% 75.0M 3s
    ##  15900K .......... .......... .......... .......... .......... 11% 69.3M 3s
    ##  15950K .......... .......... .......... .......... .......... 11% 52.1M 3s
    ##  16000K .......... .......... .......... .......... .......... 11% 58.0M 3s
    ##  16050K .......... .......... .......... .......... .......... 11% 56.3M 3s
    ##  16100K .......... .......... .......... .......... .......... 11% 58.6M 3s
    ##  16150K .......... .......... .......... .......... .......... 12% 67.1M 3s
    ##  16200K .......... .......... .......... .......... .......... 12% 62.9M 3s
    ##  16250K .......... .......... .......... .......... .......... 12% 52.2M 3s
    ##  16300K .......... .......... .......... .......... .......... 12% 53.5M 3s
    ##  16350K .......... .......... .......... .......... .......... 12% 62.3M 3s
    ##  16400K .......... .......... .......... .......... .......... 12% 46.5M 3s
    ##  16450K .......... .......... .......... .......... .......... 12% 51.7M 3s
    ##  16500K .......... .......... .......... .......... .......... 12% 51.0M 3s
    ##  16550K .......... .......... .......... .......... .......... 12% 65.5M 3s
    ##  16600K .......... .......... .......... .......... .......... 12% 73.2M 3s
    ##  16650K .......... .......... .......... .......... .......... 12% 74.5M 3s
    ##  16700K .......... .......... .......... .......... .......... 12% 58.2M 3s
    ##  16750K .......... .......... .......... .......... .......... 12% 52.3M 3s
    ##  16800K .......... .......... .......... .......... .......... 12% 53.6M 3s
    ##  16850K .......... .......... .......... .......... .......... 12% 70.9M 3s
    ##  16900K .......... .......... .......... .......... .......... 12% 59.7M 3s
    ##  16950K .......... .......... .......... .......... .......... 12% 63.1M 3s
    ##  17000K .......... .......... .......... .......... .......... 12% 63.3M 3s
    ##  17050K .......... .......... .......... .......... .......... 12% 62.9M 3s
    ##  17100K .......... .......... .......... .......... .......... 12% 59.0M 3s
    ##  17150K .......... .......... .......... .......... .......... 12% 62.0M 3s
    ##  17200K .......... .......... .......... .......... .......... 12% 56.6M 3s
    ##  17250K .......... .......... .......... .......... .......... 12% 65.6M 3s
    ##  17300K .......... .......... .......... .......... .......... 12% 73.2M 3s
    ##  17350K .......... .......... .......... .......... .......... 12% 64.1M 3s
    ##  17400K .......... .......... .......... .......... .......... 12% 60.9M 3s
    ##  17450K .......... .......... .......... .......... .......... 12% 61.3M 3s
    ##  17500K .......... .......... .......... .......... .......... 13% 57.8M 3s
    ##  17550K .......... .......... .......... .......... .......... 13% 53.8M 3s
    ##  17600K .......... .......... .......... .......... .......... 13% 60.0M 3s
    ##  17650K .......... .......... .......... .......... .......... 13% 65.4M 3s
    ##  17700K .......... .......... .......... .......... .......... 13% 63.5M 3s
    ##  17750K .......... .......... .......... .......... .......... 13% 70.0M 3s
    ##  17800K .......... .......... .......... .......... .......... 13% 66.2M 3s
    ##  17850K .......... .......... .......... .......... .......... 13% 61.1M 3s
    ##  17900K .......... .......... .......... .......... .......... 13% 64.1M 3s
    ##  17950K .......... .......... .......... .......... .......... 13% 73.6M 3s
    ##  18000K .......... .......... .......... .......... .......... 13% 65.1M 3s
    ##  18050K .......... .......... .......... .......... .......... 13% 71.5M 3s
    ##  18100K .......... .......... .......... .......... .......... 13% 66.1M 3s
    ##  18150K .......... .......... .......... .......... .......... 13% 90.5M 3s
    ##  18200K .......... .......... .......... .......... .......... 13% 68.8M 3s
    ##  18250K .......... .......... .......... .......... .......... 13% 67.4M 3s
    ##  18300K .......... .......... .......... .......... .......... 13% 61.9M 3s
    ##  18350K .......... .......... .......... .......... .......... 13% 70.5M 3s
    ##  18400K .......... .......... .......... .......... .......... 13% 59.6M 3s
    ##  18450K .......... .......... .......... .......... .......... 13% 63.3M 3s
    ##  18500K .......... .......... .......... .......... .......... 13% 61.3M 3s
    ##  18550K .......... .......... .......... .......... .......... 13% 68.9M 3s
    ##  18600K .......... .......... .......... .......... .......... 13% 71.7M 3s
    ##  18650K .......... .......... .......... .......... .......... 13% 65.4M 3s
    ##  18700K .......... .......... .......... .......... .......... 13% 60.0M 3s
    ##  18750K .......... .......... .......... .......... .......... 13% 66.3M 3s
    ##  18800K .......... .......... .......... .......... .......... 13% 56.0M 3s
    ##  18850K .......... .......... .......... .......... .......... 14% 58.9M 3s
    ##  18900K .......... .......... .......... .......... .......... 14% 69.6M 3s
    ##  18950K .......... .......... .......... .......... .......... 14% 78.1M 3s
    ##  19000K .......... .......... .......... .......... .......... 14% 57.0M 3s
    ##  19050K .......... .......... .......... .......... .......... 14% 82.9M 3s
    ##  19100K .......... .......... .......... .......... .......... 14% 70.8M 3s
    ##  19150K .......... .......... .......... .......... .......... 14% 73.6M 3s
    ##  19200K .......... .......... .......... .......... .......... 14% 89.1M 3s
    ##  19250K .......... .......... .......... .......... .......... 14% 92.7M 3s
    ##  19300K .......... .......... .......... .......... .......... 14% 69.4M 3s
    ##  19350K .......... .......... .......... .......... .......... 14% 79.7M 3s
    ##  19400K .......... .......... .......... .......... .......... 14% 72.0M 3s
    ##  19450K .......... .......... .......... .......... .......... 14% 68.3M 3s
    ##  19500K .......... .......... .......... .......... .......... 14% 69.0M 3s
    ##  19550K .......... .......... .......... .......... .......... 14% 79.5M 3s
    ##  19600K .......... .......... .......... .......... .......... 14% 73.3M 3s
    ##  19650K .......... .......... .......... .......... .......... 14% 81.1M 3s
    ##  19700K .......... .......... .......... .......... .......... 14% 71.3M 3s
    ##  19750K .......... .......... .......... .......... .......... 14% 67.3M 3s
    ##  19800K .......... .......... .......... .......... .......... 14% 68.5M 3s
    ##  19850K .......... .......... .......... .......... .......... 14% 74.5M 3s
    ##  19900K .......... .......... .......... .......... .......... 14% 62.7M 3s
    ##  19950K .......... .......... .......... .......... .......... 14% 73.3M 3s
    ##  20000K .......... .......... .......... .......... .......... 14% 72.5M 3s
    ##  20050K .......... .......... .......... .......... .......... 14% 68.6M 3s
    ##  20100K .......... .......... .......... .......... .......... 14% 72.8M 3s
    ##  20150K .......... .......... .......... .......... .......... 14% 71.0M 3s
    ##  20200K .......... .......... .......... .......... .......... 15% 73.5M 3s
    ##  20250K .......... .......... .......... .......... .......... 15% 67.7M 3s
    ##  20300K .......... .......... .......... .......... .......... 15% 64.7M 3s
    ##  20350K .......... .......... .......... .......... .......... 15% 68.3M 3s
    ##  20400K .......... .......... .......... .......... .......... 15% 67.0M 3s
    ##  20450K .......... .......... .......... .......... .......... 15% 68.4M 3s
    ##  20500K .......... .......... .......... .......... .......... 15% 63.1M 3s
    ##  20550K .......... .......... .......... .......... .......... 15% 72.4M 3s
    ##  20600K .......... .......... .......... .......... .......... 15% 68.6M 3s
    ##  20650K .......... .......... .......... .......... .......... 15% 74.5M 3s
    ##  20700K .......... .......... .......... .......... .......... 15% 78.3M 3s
    ##  20750K .......... .......... .......... .......... .......... 15% 70.9M 3s
    ##  20800K .......... .......... .......... .......... .......... 15% 69.3M 3s
    ##  20850K .......... .......... .......... .......... .......... 15% 71.7M 3s
    ##  20900K .......... .......... .......... .......... .......... 15% 69.1M 3s
    ##  20950K .......... .......... .......... .......... .......... 15% 66.8M 3s
    ##  21000K .......... .......... .......... .......... .......... 15% 78.1M 3s
    ##  21050K .......... .......... .......... .......... .......... 15% 69.2M 3s
    ##  21100K .......... .......... .......... .......... .......... 15% 75.1M 3s
    ##  21150K .......... .......... .......... .......... .......... 15% 75.8M 3s
    ##  21200K .......... .......... .......... .......... .......... 15% 61.0M 3s
    ##  21250K .......... .......... .......... .......... .......... 15% 79.7M 3s
    ##  21300K .......... .......... .......... .......... .......... 15% 75.2M 3s
    ##  21350K .......... .......... .......... .......... .......... 15% 88.7M 3s
    ##  21400K .......... .......... .......... .......... .......... 15% 74.7M 3s
    ##  21450K .......... .......... .......... .......... .......... 15% 71.5M 3s
    ##  21500K .......... .......... .......... .......... .......... 15% 71.5M 3s
    ##  21550K .......... .......... .......... .......... .......... 16% 93.3M 3s
    ##  21600K .......... .......... .......... .......... .......... 16% 68.0M 3s
    ##  21650K .......... .......... .......... .......... .......... 16% 76.9M 3s
    ##  21700K .......... .......... .......... .......... .......... 16% 73.1M 3s
    ##  21750K .......... .......... .......... .......... .......... 16% 66.1M 3s
    ##  21800K .......... .......... .......... .......... .......... 16% 65.1M 3s
    ##  21850K .......... .......... .......... .......... .......... 16% 78.5M 3s
    ##  21900K .......... .......... .......... .......... .......... 16% 80.0M 3s
    ##  21950K .......... .......... .......... .......... .......... 16% 81.3M 3s
    ##  22000K .......... .......... .......... .......... .......... 16% 86.2M 3s
    ##  22050K .......... .......... .......... .......... .......... 16% 82.1M 3s
    ##  22100K .......... .......... .......... .......... .......... 16% 62.9M 3s
    ##  22150K .......... .......... .......... .......... .......... 16% 68.8M 3s
    ##  22200K .......... .......... .......... .......... .......... 16% 76.7M 3s
    ##  22250K .......... .......... .......... .......... .......... 16% 94.2M 3s
    ##  22300K .......... .......... .......... .......... .......... 16% 54.7M 3s
    ##  22350K .......... .......... .......... .......... .......... 16% 72.1M 3s
    ##  22400K .......... .......... .......... .......... .......... 16% 28.9M 3s
    ##  22450K .......... .......... .......... .......... .......... 16% 44.5M 3s
    ##  22500K .......... .......... .......... .......... .......... 16% 35.8M 3s
    ##  22550K .......... .......... .......... .......... .......... 16% 38.2M 3s
    ##  22600K .......... .......... .......... .......... .......... 16% 36.6M 3s
    ##  22650K .......... .......... .......... .......... .......... 16% 36.7M 3s
    ##  22700K .......... .......... .......... .......... .......... 16% 37.9M 3s
    ##  22750K .......... .......... .......... .......... .......... 16% 36.8M 3s
    ##  22800K .......... .......... .......... .......... .......... 16% 39.8M 3s
    ##  22850K .......... .......... .......... .......... .......... 16% 42.7M 3s
    ##  22900K .......... .......... .......... .......... .......... 17% 51.1M 3s
    ##  22950K .......... .......... .......... .......... .......... 17% 58.4M 3s
    ##  23000K .......... .......... .......... .......... .......... 17% 43.2M 3s
    ##  23050K .......... .......... .......... .......... .......... 17% 43.8M 3s
    ##  23100K .......... .......... .......... .......... .......... 17% 48.6M 3s
    ##  23150K .......... .......... .......... .......... .......... 17% 52.1M 3s
    ##  23200K .......... .......... .......... .......... .......... 17% 44.9M 3s
    ##  23250K .......... .......... .......... .......... .......... 17% 48.5M 3s
    ##  23300K .......... .......... .......... .......... .......... 17% 39.1M 3s
    ##  23350K .......... .......... .......... .......... .......... 17% 40.1M 3s
    ##  23400K .......... .......... .......... .......... .......... 17% 46.6M 3s
    ##  23450K .......... .......... .......... .......... .......... 17% 58.1M 3s
    ##  23500K .......... .......... .......... .......... .......... 17% 57.4M 3s
    ##  23550K .......... .......... .......... .......... .......... 17% 51.5M 3s
    ##  23600K .......... .......... .......... .......... .......... 17% 50.3M 3s
    ##  23650K .......... .......... .......... .......... .......... 17% 40.5M 3s
    ##  23700K .......... .......... .......... .......... .......... 17% 50.8M 3s
    ##  23750K .......... .......... .......... .......... .......... 17% 51.9M 3s
    ##  23800K .......... .......... .......... .......... .......... 17% 47.1M 3s
    ##  23850K .......... .......... .......... .......... .......... 17% 58.3M 3s
    ##  23900K .......... .......... .......... .......... .......... 17% 49.9M 3s
    ##  23950K .......... .......... .......... .......... .......... 17% 57.3M 3s
    ##  24000K .......... .......... .......... .......... .......... 17% 53.0M 3s
    ##  24050K .......... .......... .......... .......... .......... 17% 50.7M 3s
    ##  24100K .......... .......... .......... .......... .......... 17% 47.0M 3s
    ##  24150K .......... .......... .......... .......... .......... 17% 51.3M 3s
    ##  24200K .......... .......... .......... .......... .......... 17% 49.2M 3s
    ##  24250K .......... .......... .......... .......... .......... 18% 57.3M 3s
    ##  24300K .......... .......... .......... .......... .......... 18% 52.5M 3s
    ##  24350K .......... .......... .......... .......... .......... 18% 56.9M 3s
    ##  24400K .......... .......... .......... .......... .......... 18% 52.5M 3s
    ##  24450K .......... .......... .......... .......... .......... 18% 52.1M 3s
    ##  24500K .......... .......... .......... .......... .......... 18% 48.1M 3s
    ##  24550K .......... .......... .......... .......... .......... 18% 70.0M 3s
    ##  24600K .......... .......... .......... .......... .......... 18% 48.6M 3s
    ##  24650K .......... .......... .......... .......... .......... 18% 59.5M 3s
    ##  24700K .......... .......... .......... .......... .......... 18% 52.3M 3s
    ##  24750K .......... .......... .......... .......... .......... 18% 64.4M 3s
    ##  24800K .......... .......... .......... .......... .......... 18% 54.8M 3s
    ##  24850K .......... .......... .......... .......... .......... 18% 68.3M 3s
    ##  24900K .......... .......... .......... .......... .......... 18% 53.4M 3s
    ##  24950K .......... .......... .......... .......... .......... 18% 50.3M 3s
    ##  25000K .......... .......... .......... .......... .......... 18% 48.5M 3s
    ##  25050K .......... .......... .......... .......... .......... 18% 55.8M 3s
    ##  25100K .......... .......... .......... .......... .......... 18% 50.7M 3s
    ##  25150K .......... .......... .......... .......... .......... 18% 62.4M 3s
    ##  25200K .......... .......... .......... .......... .......... 18% 55.8M 3s
    ##  25250K .......... .......... .......... .......... .......... 18% 59.8M 3s
    ##  25300K .......... .......... .......... .......... .......... 18% 61.4M 2s
    ##  25350K .......... .......... .......... .......... .......... 18% 52.8M 2s
    ##  25400K .......... .......... .......... .......... .......... 18% 53.5M 2s
    ##  25450K .......... .......... .......... .......... .......... 18% 74.0M 2s
    ##  25500K .......... .......... .......... .......... .......... 18% 61.9M 2s
    ##  25550K .......... .......... .......... .......... .......... 18% 65.8M 2s
    ##  25600K .......... .......... .......... .......... .......... 19% 51.5M 2s
    ##  25650K .......... .......... .......... .......... .......... 19% 50.5M 2s
    ##  25700K .......... .......... .......... .......... .......... 19% 56.1M 2s
    ##  25750K .......... .......... .......... .......... .......... 19% 76.4M 2s
    ##  25800K .......... .......... .......... .......... .......... 19% 54.0M 2s
    ##  25850K .......... .......... .......... .......... .......... 19% 61.6M 2s
    ##  25900K .......... .......... .......... .......... .......... 19% 56.6M 2s
    ##  25950K .......... .......... .......... .......... .......... 19% 51.6M 2s
    ##  26000K .......... .......... .......... .......... .......... 19% 54.3M 2s
    ##  26050K .......... .......... .......... .......... .......... 19% 59.5M 2s
    ##  26100K .......... .......... .......... .......... .......... 19% 55.8M 2s
    ##  26150K .......... .......... .......... .......... .......... 19% 61.5M 2s
    ##  26200K .......... .......... .......... .......... .......... 19% 60.2M 2s
    ##  26250K .......... .......... .......... .......... .......... 19% 78.4M 2s
    ##  26300K .......... .......... .......... .......... .......... 19% 63.7M 2s
    ##  26350K .......... .......... .......... .......... .......... 19% 81.6M 2s
    ##  26400K .......... .......... .......... .......... .......... 19% 55.0M 2s
    ##  26450K .......... .......... .......... .......... .......... 19% 84.9M 2s
    ##  26500K .......... .......... .......... .......... .......... 19% 53.3M 2s
    ##  26550K .......... .......... .......... .......... .......... 19% 67.8M 2s
    ##  26600K .......... .......... .......... .......... .......... 19% 63.1M 2s
    ##  26650K .......... .......... .......... .......... .......... 19% 75.4M 2s
    ##  26700K .......... .......... .......... .......... .......... 19% 58.1M 2s
    ##  26750K .......... .......... .......... .......... .......... 19% 84.2M 2s
    ##  26800K .......... .......... .......... .......... .......... 19% 56.6M 2s
    ##  26850K .......... .......... .......... .......... .......... 19% 74.6M 2s
    ##  26900K .......... .......... .......... .......... .......... 20% 70.7M 2s
    ##  26950K .......... .......... .......... .......... .......... 20% 82.5M 2s
    ##  27000K .......... .......... .......... .......... .......... 20% 74.5M 2s
    ##  27050K .......... .......... .......... .......... .......... 20% 64.4M 2s
    ##  27100K .......... .......... .......... .......... .......... 20% 68.6M 2s
    ##  27150K .......... .......... .......... .......... .......... 20% 75.3M 2s
    ##  27200K .......... .......... .......... .......... .......... 20% 58.3M 2s
    ##  27250K .......... .......... .......... .......... .......... 20% 73.8M 2s
    ##  27300K .......... .......... .......... .......... .......... 20% 57.6M 2s
    ##  27350K .......... .......... .......... .......... .......... 20% 69.1M 2s
    ##  27400K .......... .......... .......... .......... .......... 20% 61.5M 2s
    ##  27450K .......... .......... .......... .......... .......... 20% 63.9M 2s
    ##  27500K .......... .......... .......... .......... .......... 20% 57.3M 2s
    ##  27550K .......... .......... .......... .......... .......... 20% 71.4M 2s
    ##  27600K .......... .......... .......... .......... .......... 20% 67.6M 2s
    ##  27650K .......... .......... .......... .......... .......... 20% 74.0M 2s
    ##  27700K .......... .......... .......... .......... .......... 20% 80.3M 2s
    ##  27750K .......... .......... .......... .......... .......... 20% 73.9M 2s
    ##  27800K .......... .......... .......... .......... .......... 20% 80.2M 2s
    ##  27850K .......... .......... .......... .......... .......... 20% 87.8M 2s
    ##  27900K .......... .......... .......... .......... .......... 20% 36.0M 2s
    ##  27950K .......... .......... .......... .......... .......... 20% 62.5M 2s
    ##  28000K .......... .......... .......... .......... .......... 20% 71.2M 2s
    ##  28050K .......... .......... .......... .......... .......... 20% 97.9M 2s
    ##  28100K .......... .......... .......... .......... .......... 20% 62.8M 2s
    ##  28150K .......... .......... .......... .......... .......... 20% 86.0M 2s
    ##  28200K .......... .......... .......... .......... .......... 20% 92.2M 2s
    ##  28250K .......... .......... .......... .......... .......... 21% 59.8M 2s
    ##  28300K .......... .......... .......... .......... .......... 21% 39.3M 2s
    ##  28350K .......... .......... .......... .......... .......... 21% 70.6M 2s
    ##  28400K .......... .......... .......... .......... .......... 21% 60.9M 2s
    ##  28450K .......... .......... .......... .......... .......... 21% 70.5M 2s
    ##  28500K .......... .......... .......... .......... .......... 21% 63.5M 2s
    ##  28550K .......... .......... .......... .......... .......... 21% 59.6M 2s
    ##  28600K .......... .......... .......... .......... .......... 21% 45.6M 2s
    ##  28650K .......... .......... .......... .......... .......... 21% 65.7M 2s
    ##  28700K .......... .......... .......... .......... .......... 21% 58.5M 2s
    ##  28750K .......... .......... .......... .......... .......... 21% 71.9M 2s
    ##  28800K .......... .......... .......... .......... .......... 21% 71.9M 2s
    ##  28850K .......... .......... .......... .......... .......... 21% 25.0M 2s
    ##  28900K .......... .......... .......... .......... .......... 21% 49.1M 2s
    ##  28950K .......... .......... .......... .......... .......... 21% 66.8M 2s
    ##  29000K .......... .......... .......... .......... .......... 21% 84.2M 2s
    ##  29050K .......... .......... .......... .......... .......... 21% 86.5M 2s
    ##  29100K .......... .......... .......... .......... .......... 21% 86.6M 2s
    ##  29150K .......... .......... .......... .......... .......... 21% 29.6M 2s
    ##  29200K .......... .......... .......... .......... .......... 21% 63.3M 2s
    ##  29250K .......... .......... .......... .......... .......... 21% 76.2M 2s
    ##  29300K .......... .......... .......... .......... .......... 21% 66.1M 2s
    ##  29350K .......... .......... .......... .......... .......... 21% 75.2M 2s
    ##  29400K .......... .......... .......... .......... .......... 21% 78.9M 2s
    ##  29450K .......... .......... .......... .......... .......... 21% 61.1M 2s
    ##  29500K .......... .......... .......... .......... .......... 21% 58.3M 2s
    ##  29550K .......... .......... .......... .......... .......... 21% 62.7M 2s
    ##  29600K .......... .......... .......... .......... .......... 22% 55.7M 2s
    ##  29650K .......... .......... .......... .......... .......... 22% 67.2M 2s
    ##  29700K .......... .......... .......... .......... .......... 22% 63.7M 2s
    ##  29750K .......... .......... .......... .......... .......... 22% 73.3M 2s
    ##  29800K .......... .......... .......... .......... .......... 22% 72.6M 2s
    ##  29850K .......... .......... .......... .......... .......... 22% 53.5M 2s
    ##  29900K .......... .......... .......... .......... .......... 22% 69.1M 2s
    ##  29950K .......... .......... .......... .......... .......... 22% 64.9M 2s
    ##  30000K .......... .......... .......... .......... .......... 22% 68.2M 2s
    ##  30050K .......... .......... .......... .......... .......... 22% 56.7M 2s
    ##  30100K .......... .......... .......... .......... .......... 22% 54.7M 2s
    ##  30150K .......... .......... .......... .......... .......... 22% 78.7M 2s
    ##  30200K .......... .......... .......... .......... .......... 22% 52.1M 2s
    ##  30250K .......... .......... .......... .......... .......... 22% 71.7M 2s
    ##  30300K .......... .......... .......... .......... .......... 22% 77.6M 2s
    ##  30350K .......... .......... .......... .......... .......... 22% 49.6M 2s
    ##  30400K .......... .......... .......... .......... .......... 22% 63.4M 2s
    ##  30450K .......... .......... .......... .......... .......... 22% 73.5M 2s
    ##  30500K .......... .......... .......... .......... .......... 22% 52.1M 2s
    ##  30550K .......... .......... .......... .......... .......... 22% 91.5M 2s
    ##  30600K .......... .......... .......... .......... .......... 22% 60.5M 2s
    ##  30650K .......... .......... .......... .......... .......... 22% 36.4M 2s
    ##  30700K .......... .......... .......... .......... .......... 22% 66.3M 2s
    ##  30750K .......... .......... .......... .......... .......... 22% 60.1M 2s
    ##  30800K .......... .......... .......... .......... .......... 22% 64.2M 2s
    ##  30850K .......... .......... .......... .......... .......... 22% 81.5M 2s
    ##  30900K .......... .......... .......... .......... .......... 22% 71.5M 2s
    ##  30950K .......... .......... .......... .......... .......... 23% 66.1M 2s
    ##  31000K .......... .......... .......... .......... .......... 23% 68.8M 2s
    ##  31050K .......... .......... .......... .......... .......... 23% 54.0M 2s
    ##  31100K .......... .......... .......... .......... .......... 23% 60.9M 2s
    ##  31150K .......... .......... .......... .......... .......... 23% 81.5M 2s
    ##  31200K .......... .......... .......... .......... .......... 23% 66.7M 2s
    ##  31250K .......... .......... .......... .......... .......... 23% 19.5M 2s
    ##  31300K .......... .......... .......... .......... .......... 23% 59.9M 2s
    ##  31350K .......... .......... .......... .......... .......... 23% 71.7M 2s
    ##  31400K .......... .......... .......... .......... .......... 23% 62.9M 2s
    ##  31450K .......... .......... .......... .......... .......... 23% 75.8M 2s
    ##  31500K .......... .......... .......... .......... .......... 23% 72.5M 2s
    ##  31550K .......... .......... .......... .......... .......... 23% 69.7M 2s
    ##  31600K .......... .......... .......... .......... .......... 23% 33.8M 2s
    ##  31650K .......... .......... .......... .......... .......... 23% 53.0M 2s
    ##  31700K .......... .......... .......... .......... .......... 23% 53.1M 2s
    ##  31750K .......... .......... .......... .......... .......... 23% 70.0M 2s
    ##  31800K .......... .......... .......... .......... .......... 23% 78.5M 2s
    ##  31850K .......... .......... .......... .......... .......... 23% 72.3M 2s
    ##  31900K .......... .......... .......... .......... .......... 23% 62.5M 2s
    ##  31950K .......... .......... .......... .......... .......... 23% 73.3M 2s
    ##  32000K .......... .......... .......... .......... .......... 23% 60.4M 2s
    ##  32050K .......... .......... .......... .......... .......... 23% 53.0M 2s
    ##  32100K .......... .......... .......... .......... .......... 23% 61.6M 2s
    ##  32150K .......... .......... .......... .......... .......... 23% 73.6M 2s
    ##  32200K .......... .......... .......... .......... .......... 23% 85.6M 2s
    ##  32250K .......... .......... .......... .......... .......... 23% 24.5M 2s
    ##  32300K .......... .......... .......... .......... .......... 24% 63.6M 2s
    ##  32350K .......... .......... .......... .......... .......... 24% 68.2M 2s
    ##  32400K .......... .......... .......... .......... .......... 24% 87.1M 2s
    ##  32450K .......... .......... .......... .......... .......... 24%  100M 2s
    ##  32500K .......... .......... .......... .......... .......... 24% 76.6M 2s
    ##  32550K .......... .......... .......... .......... .......... 24% 34.4M 2s
    ##  32600K .......... .......... .......... .......... .......... 24% 44.1M 2s
    ##  32650K .......... .......... .......... .......... .......... 24% 63.4M 2s
    ##  32700K .......... .......... .......... .......... .......... 24% 69.8M 2s
    ##  32750K .......... .......... .......... .......... .......... 24% 81.2M 2s
    ##  32800K .......... .......... .......... .......... .......... 24% 65.8M 2s
    ##  32850K .......... .......... .......... .......... .......... 24% 13.9M 2s
    ##  32900K .......... .......... .......... .......... .......... 24% 55.6M 2s
    ##  32950K .......... .......... .......... .......... .......... 24% 64.2M 2s
    ##  33000K .......... .......... .......... .......... .......... 24% 70.6M 2s
    ##  33050K .......... .......... .......... .......... .......... 24% 68.9M 2s
    ##  33100K .......... .......... .......... .......... .......... 24% 78.2M 2s
    ##  33150K .......... .......... .......... .......... .......... 24% 79.1M 2s
    ##  33200K .......... .......... .......... .......... .......... 24% 27.2M 2s
    ##  33250K .......... .......... .......... .......... .......... 24% 26.5M 2s
    ##  33300K .......... .......... .......... .......... .......... 24% 64.1M 2s
    ##  33350K .......... .......... .......... .......... .......... 24% 65.0M 2s
    ##  33400K .......... .......... .......... .......... .......... 24% 71.5M 2s
    ##  33450K .......... .......... .......... .......... .......... 24% 61.3M 2s
    ##  33500K .......... .......... .......... .......... .......... 24% 57.3M 2s
    ##  33550K .......... .......... .......... .......... .......... 24% 82.0M 2s
    ##  33600K .......... .......... .......... .......... .......... 24% 74.2M 2s
    ##  33650K .......... .......... .......... .......... .......... 25% 71.1M 2s
    ##  33700K .......... .......... .......... .......... .......... 25% 68.2M 2s
    ##  33750K .......... .......... .......... .......... .......... 25% 91.4M 2s
    ##  33800K .......... .......... .......... .......... .......... 25% 88.3M 2s
    ##  33850K .......... .......... .......... .......... .......... 25% 83.5M 2s
    ##  33900K .......... .......... .......... .......... .......... 25% 72.5M 2s
    ##  33950K .......... .......... .......... .......... .......... 25% 33.0M 2s
    ##  34000K .......... .......... .......... .......... .......... 25% 84.2M 2s
    ##  34050K .......... .......... .......... .......... .......... 25% 65.3M 2s
    ##  34100K .......... .......... .......... .......... .......... 25% 69.5M 2s
    ##  34150K .......... .......... .......... .......... .......... 25% 78.8M 2s
    ##  34200K .......... .......... .......... .......... .......... 25% 83.1M 2s
    ##  34250K .......... .......... .......... .......... .......... 25% 78.3M 2s
    ##  34300K .......... .......... .......... .......... .......... 25% 84.4M 2s
    ##  34350K .......... .......... .......... .......... .......... 25% 44.7M 2s
    ##  34400K .......... .......... .......... .......... .......... 25% 65.3M 2s
    ##  34450K .......... .......... .......... .......... .......... 25% 46.8M 2s
    ##  34500K .......... .......... .......... .......... .......... 25% 81.1M 2s
    ##  34550K .......... .......... .......... .......... .......... 25% 83.3M 2s
    ##  34600K .......... .......... .......... .......... .......... 25% 74.0M 2s
    ##  34650K .......... .......... .......... .......... .......... 25% 77.1M 2s
    ##  34700K .......... .......... .......... .......... .......... 25% 21.9M 2s
    ##  34750K .......... .......... .......... .......... .......... 25% 65.0M 2s
    ##  34800K .......... .......... .......... .......... .......... 25% 31.1M 2s
    ##  34850K .......... .......... .......... .......... .......... 25% 52.2M 2s
    ##  34900K .......... .......... .......... .......... .......... 25% 61.9M 2s
    ##  34950K .......... .......... .......... .......... .......... 25% 50.8M 2s
    ##  35000K .......... .......... .......... .......... .......... 26% 65.6M 2s
    ##  35050K .......... .......... .......... .......... .......... 26% 82.4M 2s
    ##  35100K .......... .......... .......... .......... .......... 26% 65.6M 2s
    ##  35150K .......... .......... .......... .......... .......... 26% 61.8M 2s
    ##  35200K .......... .......... .......... .......... .......... 26% 66.9M 2s
    ##  35250K .......... .......... .......... .......... .......... 26% 64.4M 2s
    ##  35300K .......... .......... .......... .......... .......... 26% 78.4M 2s
    ##  35350K .......... .......... .......... .......... .......... 26% 96.6M 2s
    ##  35400K .......... .......... .......... .......... .......... 26% 23.0M 2s
    ##  35450K .......... .......... .......... .......... .......... 26% 59.3M 2s
    ##  35500K .......... .......... .......... .......... .......... 26% 75.9M 2s
    ##  35550K .......... .......... .......... .......... .......... 26% 62.6M 2s
    ##  35600K .......... .......... .......... .......... .......... 26% 65.7M 2s
    ##  35650K .......... .......... .......... .......... .......... 26% 92.5M 2s
    ##  35700K .......... .......... .......... .......... .......... 26% 81.5M 2s
    ##  35750K .......... .......... .......... .......... .......... 26% 67.4M 2s
    ##  35800K .......... .......... .......... .......... .......... 26% 68.9M 2s
    ##  35850K .......... .......... .......... .......... .......... 26% 76.0M 2s
    ##  35900K .......... .......... .......... .......... .......... 26% 74.8M 2s
    ##  35950K .......... .......... .......... .......... .......... 26% 78.3M 2s
    ##  36000K .......... .......... .......... .......... .......... 26% 81.6M 2s
    ##  36050K .......... .......... .......... .......... .......... 26% 66.8M 2s
    ##  36100K .......... .......... .......... .......... .......... 26% 69.1M 2s
    ##  36150K .......... .......... .......... .......... .......... 26% 81.3M 2s
    ##  36200K .......... .......... .......... .......... .......... 26% 77.5M 2s
    ##  36250K .......... .......... .......... .......... .......... 26% 82.9M 2s
    ##  36300K .......... .......... .......... .......... .......... 26% 75.6M 2s
    ##  36350K .......... .......... .......... .......... .......... 27% 78.2M 2s
    ##  36400K .......... .......... .......... .......... .......... 27% 73.7M 2s
    ##  36450K .......... .......... .......... .......... .......... 27% 67.0M 2s
    ##  36500K .......... .......... .......... .......... .......... 27% 69.4M 2s
    ##  36550K .......... .......... .......... .......... .......... 27% 61.5M 2s
    ##  36600K .......... .......... .......... .......... .......... 27% 63.2M 2s
    ##  36650K .......... .......... .......... .......... .......... 27% 71.7M 2s
    ##  36700K .......... .......... .......... .......... .......... 27% 66.7M 2s
    ##  36750K .......... .......... .......... .......... .......... 27% 74.5M 2s
    ##  36800K .......... .......... .......... .......... .......... 27% 82.5M 2s
    ##  36850K .......... .......... .......... .......... .......... 27% 89.6M 2s
    ##  36900K .......... .......... .......... .......... .......... 27% 72.3M 2s
    ##  36950K .......... .......... .......... .......... .......... 27% 73.8M 2s
    ##  37000K .......... .......... .......... .......... .......... 27% 68.9M 2s
    ##  37050K .......... .......... .......... .......... .......... 27% 81.9M 2s
    ##  37100K .......... .......... .......... .......... .......... 27% 78.6M 2s
    ##  37150K .......... .......... .......... .......... .......... 27% 85.7M 2s
    ##  37200K .......... .......... .......... .......... .......... 27% 74.0M 2s
    ##  37250K .......... .......... .......... .......... .......... 27% 15.9M 2s
    ##  37300K .......... .......... .......... .......... .......... 27% 71.1M 2s
    ##  37350K .......... .......... .......... .......... .......... 27% 75.0M 2s
    ##  37400K .......... .......... .......... .......... .......... 27% 60.6M 2s
    ##  37450K .......... .......... .......... .......... .......... 27% 82.3M 2s
    ##  37500K .......... .......... .......... .......... .......... 27% 58.5M 2s
    ##  37550K .......... .......... .......... .......... .......... 27% 72.7M 2s
    ##  37600K .......... .......... .......... .......... .......... 27% 63.4M 2s
    ##  37650K .......... .......... .......... .......... .......... 27% 89.8M 2s
    ##  37700K .......... .......... .......... .......... .......... 28% 70.8M 2s
    ##  37750K .......... .......... .......... .......... .......... 28% 69.2M 2s
    ##  37800K .......... .......... .......... .......... .......... 28% 68.4M 2s
    ##  37850K .......... .......... .......... .......... .......... 28% 71.2M 2s
    ##  37900K .......... .......... .......... .......... .......... 28% 67.0M 2s
    ##  37950K .......... .......... .......... .......... .......... 28% 55.8M 2s
    ##  38000K .......... .......... .......... .......... .......... 28% 65.4M 2s
    ##  38050K .......... .......... .......... .......... .......... 28% 72.0M 2s
    ##  38100K .......... .......... .......... .......... .......... 28% 56.1M 2s
    ##  38150K .......... .......... .......... .......... .......... 28% 78.5M 2s
    ##  38200K .......... .......... .......... .......... .......... 28% 65.0M 2s
    ##  38250K .......... .......... .......... .......... .......... 28% 85.0M 2s
    ##  38300K .......... .......... .......... .......... .......... 28% 72.2M 2s
    ##  38350K .......... .......... .......... .......... .......... 28% 66.2M 2s
    ##  38400K .......... .......... .......... .......... .......... 28% 63.4M 2s
    ##  38450K .......... .......... .......... .......... .......... 28% 83.4M 2s
    ##  38500K .......... .......... .......... .......... .......... 28% 66.0M 2s
    ##  38550K .......... .......... .......... .......... .......... 28% 86.2M 2s
    ##  38600K .......... .......... .......... .......... .......... 28% 76.4M 2s
    ##  38650K .......... .......... .......... .......... .......... 28% 19.4M 2s
    ##  38700K .......... .......... .......... .......... .......... 28% 13.3M 2s
    ##  38750K .......... .......... .......... .......... .......... 28% 71.4M 2s
    ##  38800K .......... .......... .......... .......... .......... 28% 65.4M 2s
    ##  38850K .......... .......... .......... .......... .......... 28% 80.6M 2s
    ##  38900K .......... .......... .......... .......... .......... 28% 71.4M 2s
    ##  38950K .......... .......... .......... .......... .......... 28% 79.7M 2s
    ##  39000K .......... .......... .......... .......... .......... 28% 63.0M 2s
    ##  39050K .......... .......... .......... .......... .......... 29% 70.3M 2s
    ##  39100K .......... .......... .......... .......... .......... 29% 56.2M 2s
    ##  39150K .......... .......... .......... .......... .......... 29% 46.0M 2s
    ##  39200K .......... .......... .......... .......... .......... 29% 84.5M 2s
    ##  39250K .......... .......... .......... .......... .......... 29% 68.5M 2s
    ##  39300K .......... .......... .......... .......... .......... 29% 61.4M 2s
    ##  39350K .......... .......... .......... .......... .......... 29% 95.3M 2s
    ##  39400K .......... .......... .......... .......... .......... 29% 66.3M 2s
    ##  39450K .......... .......... .......... .......... .......... 29% 76.3M 2s
    ##  39500K .......... .......... .......... .......... .......... 29% 27.3M 2s
    ##  39550K .......... .......... .......... .......... .......... 29% 58.7M 2s
    ##  39600K .......... .......... .......... .......... .......... 29% 66.7M 2s
    ##  39650K .......... .......... .......... .......... .......... 29% 67.8M 2s
    ##  39700K .......... .......... .......... .......... .......... 29% 75.8M 2s
    ##  39750K .......... .......... .......... .......... .......... 29% 79.5M 2s
    ##  39800K .......... .......... .......... .......... .......... 29% 86.9M 2s
    ##  39850K .......... .......... .......... .......... .......... 29% 83.6M 2s
    ##  39900K .......... .......... .......... .......... .......... 29% 73.2M 2s
    ##  39950K .......... .......... .......... .......... .......... 29% 23.9M 2s
    ##  40000K .......... .......... .......... .......... .......... 29% 64.1M 2s
    ##  40050K .......... .......... .......... .......... .......... 29% 81.5M 2s
    ##  40100K .......... .......... .......... .......... .......... 29% 75.0M 2s
    ##  40150K .......... .......... .......... .......... .......... 29% 87.3M 2s
    ##  40200K .......... .......... .......... .......... .......... 29% 86.8M 2s
    ##  40250K .......... .......... .......... .......... .......... 29% 23.2M 2s
    ##  40300K .......... .......... .......... .......... .......... 29% 63.6M 2s
    ##  40350K .......... .......... .......... .......... .......... 29% 62.4M 2s
    ##  40400K .......... .......... .......... .......... .......... 30% 77.1M 2s
    ##  40450K .......... .......... .......... .......... .......... 30% 82.1M 2s
    ##  40500K .......... .......... .......... .......... .......... 30% 68.4M 2s
    ##  40550K .......... .......... .......... .......... .......... 30% 65.4M 2s
    ##  40600K .......... .......... .......... .......... .......... 30% 65.6M 2s
    ##  40650K .......... .......... .......... .......... .......... 30% 68.5M 2s
    ##  40700K .......... .......... .......... .......... .......... 30% 75.9M 2s
    ##  40750K .......... .......... .......... .......... .......... 30% 87.7M 2s
    ##  40800K .......... .......... .......... .......... .......... 30% 65.8M 2s
    ##  40850K .......... .......... .......... .......... .......... 30% 57.9M 2s
    ##  40900K .......... .......... .......... .......... .......... 30% 41.3M 2s
    ##  40950K .......... .......... .......... .......... .......... 30% 70.3M 2s
    ##  41000K .......... .......... .......... .......... .......... 30% 70.3M 2s
    ##  41050K .......... .......... .......... .......... .......... 30% 80.0M 2s
    ##  41100K .......... .......... .......... .......... .......... 30% 66.1M 2s
    ##  41150K .......... .......... .......... .......... .......... 30% 61.4M 2s
    ##  41200K .......... .......... .......... .......... .......... 30% 58.8M 2s
    ##  41250K .......... .......... .......... .......... .......... 30% 70.2M 2s
    ##  41300K .......... .......... .......... .......... .......... 30% 56.1M 2s
    ##  41350K .......... .......... .......... .......... .......... 30% 54.7M 2s
    ##  41400K .......... .......... .......... .......... .......... 30% 61.2M 2s
    ##  41450K .......... .......... .......... .......... .......... 30% 72.7M 2s
    ##  41500K .......... .......... .......... .......... .......... 30% 78.4M 2s
    ##  41550K .......... .......... .......... .......... .......... 30% 70.0M 2s
    ##  41600K .......... .......... .......... .......... .......... 30% 75.0M 2s
    ##  41650K .......... .......... .......... .......... .......... 30% 82.0M 2s
    ##  41700K .......... .......... .......... .......... .......... 30% 69.2M 2s
    ##  41750K .......... .......... .......... .......... .......... 31% 71.6M 2s
    ##  41800K .......... .......... .......... .......... .......... 31% 68.0M 2s
    ##  41850K .......... .......... .......... .......... .......... 31% 91.1M 2s
    ##  41900K .......... .......... .......... .......... .......... 31% 57.5M 2s
    ##  41950K .......... .......... .......... .......... .......... 31% 66.2M 2s
    ##  42000K .......... .......... .......... .......... .......... 31% 75.6M 2s
    ##  42050K .......... .......... .......... .......... .......... 31% 74.3M 2s
    ##  42100K .......... .......... .......... .......... .......... 31% 70.9M 2s
    ##  42150K .......... .......... .......... .......... .......... 31% 74.3M 2s
    ##  42200K .......... .......... .......... .......... .......... 31% 78.2M 2s
    ##  42250K .......... .......... .......... .......... .......... 31% 92.1M 2s
    ##  42300K .......... .......... .......... .......... .......... 31% 71.5M 2s
    ##  42350K .......... .......... .......... .......... .......... 31% 67.2M 2s
    ##  42400K .......... .......... .......... .......... .......... 31% 84.8M 2s
    ##  42450K .......... .......... .......... .......... .......... 31% 31.6M 2s
    ##  42500K .......... .......... .......... .......... .......... 31% 82.8M 2s
    ##  42550K .......... .......... .......... .......... .......... 31% 97.1M 2s
    ##  42600K .......... .......... .......... .......... .......... 31% 65.4M 2s
    ##  42650K .......... .......... .......... .......... .......... 31% 74.1M 2s
    ##  42700K .......... .......... .......... .......... .......... 31% 65.4M 2s
    ##  42750K .......... .......... .......... .......... .......... 31% 75.1M 2s
    ##  42800K .......... .......... .......... .......... .......... 31% 57.4M 2s
    ##  42850K .......... .......... .......... .......... .......... 31% 60.8M 2s
    ##  42900K .......... .......... .......... .......... .......... 31% 69.7M 2s
    ##  42950K .......... .......... .......... .......... .......... 31% 45.5M 2s
    ##  43000K .......... .......... .......... .......... .......... 31% 72.9M 2s
    ##  43050K .......... .......... .......... .......... .......... 31% 79.9M 2s
    ##  43100K .......... .......... .......... .......... .......... 32% 68.3M 2s
    ##  43150K .......... .......... .......... .......... .......... 32% 83.8M 2s
    ##  43200K .......... .......... .......... .......... .......... 32% 70.8M 2s
    ##  43250K .......... .......... .......... .......... .......... 32% 20.3M 2s
    ##  43300K .......... .......... .......... .......... .......... 32% 66.5M 2s
    ##  43350K .......... .......... .......... .......... .......... 32%  101M 2s
    ##  43400K .......... .......... .......... .......... .......... 32% 86.5M 2s
    ##  43450K .......... .......... .......... .......... .......... 32% 90.2M 2s
    ##  43500K .......... .......... .......... .......... .......... 32% 81.4M 2s
    ##  43550K .......... .......... .......... .......... .......... 32% 87.4M 2s
    ##  43600K .......... .......... .......... .......... .......... 32% 34.2M 2s
    ##  43650K .......... .......... .......... .......... .......... 32% 72.2M 2s
    ##  43700K .......... .......... .......... .......... .......... 32% 31.2M 2s
    ##  43750K .......... .......... .......... .......... .......... 32% 72.3M 2s
    ##  43800K .......... .......... .......... .......... .......... 32% 39.7M 2s
    ##  43850K .......... .......... .......... .......... .......... 32% 75.7M 2s
    ##  43900K .......... .......... .......... .......... .......... 32% 63.1M 2s
    ##  43950K .......... .......... .......... .......... .......... 32% 98.2M 2s
    ##  44000K .......... .......... .......... .......... .......... 32% 79.0M 2s
    ##  44050K .......... .......... .......... .......... .......... 32% 13.6M 2s
    ##  44100K .......... .......... .......... .......... .......... 32% 72.4M 2s
    ##  44150K .......... .......... .......... .......... .......... 32% 34.8M 2s
    ##  44200K .......... .......... .......... .......... .......... 32% 57.7M 2s
    ##  44250K .......... .......... .......... .......... .......... 32% 78.7M 2s
    ##  44300K .......... .......... .......... .......... .......... 32% 63.3M 2s
    ##  44350K .......... .......... .......... .......... .......... 32% 91.3M 2s
    ##  44400K .......... .......... .......... .......... .......... 32% 69.9M 2s
    ##  44450K .......... .......... .......... .......... .......... 33% 83.2M 2s
    ##  44500K .......... .......... .......... .......... .......... 33% 60.4M 2s
    ##  44550K .......... .......... .......... .......... .......... 33% 92.9M 2s
    ##  44600K .......... .......... .......... .......... .......... 33% 73.0M 2s
    ##  44650K .......... .......... .......... .......... .......... 33% 78.9M 2s
    ##  44700K .......... .......... .......... .......... .......... 33% 63.9M 2s
    ##  44750K .......... .......... .......... .......... .......... 33% 65.3M 2s
    ##  44800K .......... .......... .......... .......... .......... 33% 60.7M 2s
    ##  44850K .......... .......... .......... .......... .......... 33% 72.4M 2s
    ##  44900K .......... .......... .......... .......... .......... 33% 72.5M 2s
    ##  44950K .......... .......... .......... .......... .......... 33% 70.2M 2s
    ##  45000K .......... .......... .......... .......... .......... 33% 79.1M 2s
    ##  45050K .......... .......... .......... .......... .......... 33% 68.4M 2s
    ##  45100K .......... .......... .......... .......... .......... 33% 60.7M 2s
    ##  45150K .......... .......... .......... .......... .......... 33% 79.2M 2s
    ##  45200K .......... .......... .......... .......... .......... 33% 68.1M 2s
    ##  45250K .......... .......... .......... .......... .......... 33% 92.7M 2s
    ##  45300K .......... .......... .......... .......... .......... 33% 64.5M 2s
    ##  45350K .......... .......... .......... .......... .......... 33% 64.4M 2s
    ##  45400K .......... .......... .......... .......... .......... 33% 68.6M 2s
    ##  45450K .......... .......... .......... .......... .......... 33% 67.3M 2s
    ##  45500K .......... .......... .......... .......... .......... 33% 62.3M 2s
    ##  45550K .......... .......... .......... .......... .......... 33% 75.6M 2s
    ##  45600K .......... .......... .......... .......... .......... 33% 67.7M 2s
    ##  45650K .......... .......... .......... .......... .......... 33% 78.7M 2s
    ##  45700K .......... .......... .......... .......... .......... 33% 12.1M 2s
    ##  45750K .......... .......... .......... .......... .......... 33% 26.5M 2s
    ##  45800K .......... .......... .......... .......... .......... 34% 56.0M 2s
    ##  45850K .......... .......... .......... .......... .......... 34% 72.0M 2s
    ##  45900K .......... .......... .......... .......... .......... 34% 74.8M 2s
    ##  45950K .......... .......... .......... .......... .......... 34% 79.5M 2s
    ##  46000K .......... .......... .......... .......... .......... 34% 18.2M 2s
    ##  46050K .......... .......... .......... .......... .......... 34% 40.9M 2s
    ##  46100K .......... .......... .......... .......... .......... 34% 38.5M 2s
    ##  46150K .......... .......... .......... .......... .......... 34% 67.8M 2s
    ##  46200K .......... .......... .......... .......... .......... 34% 53.1M 2s
    ##  46250K .......... .......... .......... .......... .......... 34% 69.9M 2s
    ##  46300K .......... .......... .......... .......... .......... 34% 68.2M 2s
    ##  46350K .......... .......... .......... .......... .......... 34% 83.2M 2s
    ##  46400K .......... .......... .......... .......... .......... 34% 85.2M 2s
    ##  46450K .......... .......... .......... .......... .......... 34% 24.4M 2s
    ##  46500K .......... .......... .......... .......... .......... 34% 64.7M 2s
    ##  46550K .......... .......... .......... .......... .......... 34% 65.2M 2s
    ##  46600K .......... .......... .......... .......... .......... 34% 66.5M 2s
    ##  46650K .......... .......... .......... .......... .......... 34% 76.3M 2s
    ##  46700K .......... .......... .......... .......... .......... 34% 85.8M 2s
    ##  46750K .......... .......... .......... .......... .......... 34%  103M 2s
    ##  46800K .......... .......... .......... .......... .......... 34% 70.0M 2s
    ##  46850K .......... .......... .......... .......... .......... 34% 93.6M 2s
    ##  46900K .......... .......... .......... .......... .......... 34% 47.3M 2s
    ##  46950K .......... .......... .......... .......... .......... 34% 66.6M 2s
    ##  47000K .......... .......... .......... .......... .......... 34% 59.7M 2s
    ##  47050K .......... .......... .......... .......... .......... 34% 71.1M 2s
    ##  47100K .......... .......... .......... .......... .......... 34% 64.5M 2s
    ##  47150K .......... .......... .......... .......... .......... 35% 65.5M 2s
    ##  47200K .......... .......... .......... .......... .......... 35% 79.2M 2s
    ##  47250K .......... .......... .......... .......... .......... 35% 43.0M 2s
    ##  47300K .......... .......... .......... .......... .......... 35% 60.0M 2s
    ##  47350K .......... .......... .......... .......... .......... 35% 82.5M 2s
    ##  47400K .......... .......... .......... .......... .......... 35% 63.9M 2s
    ##  47450K .......... .......... .......... .......... .......... 35% 70.4M 2s
    ##  47500K .......... .......... .......... .......... .......... 35% 68.0M 2s
    ##  47550K .......... .......... .......... .......... .......... 35% 74.7M 2s
    ##  47600K .......... .......... .......... .......... .......... 35% 67.0M 2s
    ##  47650K .......... .......... .......... .......... .......... 35% 80.8M 2s
    ##  47700K .......... .......... .......... .......... .......... 35% 74.3M 2s
    ##  47750K .......... .......... .......... .......... .......... 35% 80.7M 2s
    ##  47800K .......... .......... .......... .......... .......... 35% 58.6M 2s
    ##  47850K .......... .......... .......... .......... .......... 35% 75.8M 2s
    ##  47900K .......... .......... .......... .......... .......... 35% 67.0M 2s
    ##  47950K .......... .......... .......... .......... .......... 35% 72.8M 2s
    ##  48000K .......... .......... .......... .......... .......... 35% 66.3M 2s
    ##  48050K .......... .......... .......... .......... .......... 35% 87.1M 2s
    ##  48100K .......... .......... .......... .......... .......... 35% 84.4M 2s
    ##  48150K .......... .......... .......... .......... .......... 35% 75.4M 2s
    ##  48200K .......... .......... .......... .......... .......... 35% 69.4M 2s
    ##  48250K .......... .......... .......... .......... .......... 35% 74.1M 2s
    ##  48300K .......... .......... .......... .......... .......... 35% 65.4M 2s
    ##  48350K .......... .......... .......... .......... .......... 35% 73.2M 2s
    ##  48400K .......... .......... .......... .......... .......... 35% 67.4M 2s
    ##  48450K .......... .......... .......... .......... .......... 35% 62.8M 2s
    ##  48500K .......... .......... .......... .......... .......... 36% 64.7M 2s
    ##  48550K .......... .......... .......... .......... .......... 36% 79.8M 2s
    ##  48600K .......... .......... .......... .......... .......... 36% 71.3M 2s
    ##  48650K .......... .......... .......... .......... .......... 36% 70.1M 2s
    ##  48700K .......... .......... .......... .......... .......... 36% 62.8M 2s
    ##  48750K .......... .......... .......... .......... .......... 36% 62.6M 2s
    ##  48800K .......... .......... .......... .......... .......... 36% 77.7M 2s
    ##  48850K .......... .......... .......... .......... .......... 36% 75.6M 2s
    ##  48900K .......... .......... .......... .......... .......... 36% 67.1M 2s
    ##  48950K .......... .......... .......... .......... .......... 36% 74.4M 2s
    ##  49000K .......... .......... .......... .......... .......... 36% 61.3M 2s
    ##  49050K .......... .......... .......... .......... .......... 36% 77.2M 2s
    ##  49100K .......... .......... .......... .......... .......... 36% 71.0M 2s
    ##  49150K .......... .......... .......... .......... .......... 36% 67.0M 2s
    ##  49200K .......... .......... .......... .......... .......... 36% 79.8M 2s
    ##  49250K .......... .......... .......... .......... .......... 36% 81.4M 2s
    ##  49300K .......... .......... .......... .......... .......... 36% 72.6M 2s
    ##  49350K .......... .......... .......... .......... .......... 36% 77.8M 2s
    ##  49400K .......... .......... .......... .......... .......... 36% 90.2M 2s
    ##  49450K .......... .......... .......... .......... .......... 36% 74.0M 2s
    ##  49500K .......... .......... .......... .......... .......... 36% 80.7M 2s
    ##  49550K .......... .......... .......... .......... .......... 36% 82.0M 2s
    ##  49600K .......... .......... .......... .......... .......... 36% 75.0M 2s
    ##  49650K .......... .......... .......... .......... .......... 36% 76.6M 2s
    ##  49700K .......... .......... .......... .......... .......... 36% 63.3M 2s
    ##  49750K .......... .......... .......... .......... .......... 36% 56.4M 2s
    ##  49800K .......... .......... .......... .......... .......... 36% 67.7M 2s
    ##  49850K .......... .......... .......... .......... .......... 37% 81.1M 2s
    ##  49900K .......... .......... .......... .......... .......... 37% 69.6M 2s
    ##  49950K .......... .......... .......... .......... .......... 37% 83.9M 2s
    ##  50000K .......... .......... .......... .......... .......... 37% 71.4M 2s
    ##  50050K .......... .......... .......... .......... .......... 37% 69.6M 2s
    ##  50100K .......... .......... .......... .......... .......... 37% 69.0M 2s
    ##  50150K .......... .......... .......... .......... .......... 37% 80.6M 2s
    ##  50200K .......... .......... .......... .......... .......... 37% 68.7M 2s
    ##  50250K .......... .......... .......... .......... .......... 37% 83.3M 2s
    ##  50300K .......... .......... .......... .......... .......... 37% 75.5M 2s
    ##  50350K .......... .......... .......... .......... .......... 37% 70.3M 2s
    ##  50400K .......... .......... .......... .......... .......... 37% 75.9M 2s
    ##  50450K .......... .......... .......... .......... .......... 37% 74.3M 2s
    ##  50500K .......... .......... .......... .......... .......... 37% 58.5M 2s
    ##  50550K .......... .......... .......... .......... .......... 37% 76.8M 2s
    ##  50600K .......... .......... .......... .......... .......... 37% 76.1M 2s
    ##  50650K .......... .......... .......... .......... .......... 37% 64.1M 2s
    ##  50700K .......... .......... .......... .......... .......... 37% 69.2M 2s
    ##  50750K .......... .......... .......... .......... .......... 37% 74.5M 2s
    ##  50800K .......... .......... .......... .......... .......... 37% 77.5M 2s
    ##  50850K .......... .......... .......... .......... .......... 37% 77.7M 2s
    ##  50900K .......... .......... .......... .......... .......... 37% 82.1M 2s
    ##  50950K .......... .......... .......... .......... .......... 37% 46.1M 2s
    ##  51000K .......... .......... .......... .......... .......... 37% 88.5M 2s
    ##  51050K .......... .......... .......... .......... .......... 37% 89.1M 2s
    ##  51100K .......... .......... .......... .......... .......... 37% 64.0M 2s
    ##  51150K .......... .......... .......... .......... .......... 37% 76.4M 2s
    ##  51200K .......... .......... .......... .......... .......... 38% 69.9M 2s
    ##  51250K .......... .......... .......... .......... .......... 38% 63.4M 2s
    ##  51300K .......... .......... .......... .......... .......... 38% 64.5M 2s
    ##  51350K .......... .......... .......... .......... .......... 38% 74.6M 2s
    ##  51400K .......... .......... .......... .......... .......... 38% 73.3M 2s
    ##  51450K .......... .......... .......... .......... .......... 38% 83.7M 2s
    ##  51500K .......... .......... .......... .......... .......... 38% 59.5M 2s
    ##  51550K .......... .......... .......... .......... .......... 38% 63.3M 2s
    ##  51600K .......... .......... .......... .......... .......... 38% 70.8M 2s
    ##  51650K .......... .......... .......... .......... .......... 38% 79.9M 2s
    ##  51700K .......... .......... .......... .......... .......... 38% 67.6M 2s
    ##  51750K .......... .......... .......... .......... .......... 38% 78.1M 2s
    ##  51800K .......... .......... .......... .......... .......... 38% 70.6M 2s
    ##  51850K .......... .......... .......... .......... .......... 38% 54.9M 2s
    ##  51900K .......... .......... .......... .......... .......... 38% 67.5M 2s
    ##  51950K .......... .......... .......... .......... .......... 38% 75.1M 2s
    ##  52000K .......... .......... .......... .......... .......... 38% 59.1M 2s
    ##  52050K .......... .......... .......... .......... .......... 38% 59.1M 2s
    ##  52100K .......... .......... .......... .......... .......... 38% 74.4M 2s
    ##  52150K .......... .......... .......... .......... .......... 38% 90.5M 2s
    ##  52200K .......... .......... .......... .......... .......... 38% 77.4M 2s
    ##  52250K .......... .......... .......... .......... .......... 38% 76.8M 2s
    ##  52300K .......... .......... .......... .......... .......... 38% 79.7M 2s
    ##  52350K .......... .......... .......... .......... .......... 38% 69.6M 2s
    ##  52400K .......... .......... .......... .......... .......... 38% 87.0M 2s
    ##  52450K .......... .......... .......... .......... .......... 38% 67.7M 2s
    ##  52500K .......... .......... .......... .......... .......... 39% 77.2M 2s
    ##  52550K .......... .......... .......... .......... .......... 39%  102M 2s
    ##  52600K .......... .......... .......... .......... .......... 39% 58.6M 2s
    ##  52650K .......... .......... .......... .......... .......... 39% 62.1M 2s
    ##  52700K .......... .......... .......... .......... .......... 39% 44.6M 2s
    ##  52750K .......... .......... .......... .......... .......... 39% 70.2M 2s
    ##  52800K .......... .......... .......... .......... .......... 39% 64.1M 2s
    ##  52850K .......... .......... .......... .......... .......... 39% 70.7M 2s
    ##  52900K .......... .......... .......... .......... .......... 39% 72.7M 2s
    ##  52950K .......... .......... .......... .......... .......... 39% 56.6M 2s
    ##  53000K .......... .......... .......... .......... .......... 39% 80.3M 2s
    ##  53050K .......... .......... .......... .......... .......... 39% 84.8M 2s
    ##  53100K .......... .......... .......... .......... .......... 39% 57.5M 2s
    ##  53150K .......... .......... .......... .......... .......... 39% 71.5M 2s
    ##  53200K .......... .......... .......... .......... .......... 39% 69.7M 2s
    ##  53250K .......... .......... .......... .......... .......... 39% 55.5M 2s
    ##  53300K .......... .......... .......... .......... .......... 39% 70.0M 2s
    ##  53350K .......... .......... .......... .......... .......... 39% 69.3M 2s
    ##  53400K .......... .......... .......... .......... .......... 39% 64.7M 2s
    ##  53450K .......... .......... .......... .......... .......... 39% 77.6M 2s
    ##  53500K .......... .......... .......... .......... .......... 39% 41.9M 2s
    ##  53550K .......... .......... .......... .......... .......... 39% 60.1M 2s
    ##  53600K .......... .......... .......... .......... .......... 39% 69.6M 2s
    ##  53650K .......... .......... .......... .......... .......... 39% 79.2M 2s
    ##  53700K .......... .......... .......... .......... .......... 39% 72.7M 2s
    ##  53750K .......... .......... .......... .......... .......... 39% 75.8M 2s
    ##  53800K .......... .......... .......... .......... .......... 39% 68.9M 2s
    ##  53850K .......... .......... .......... .......... .......... 40% 72.0M 2s
    ##  53900K .......... .......... .......... .......... .......... 40% 79.3M 2s
    ##  53950K .......... .......... .......... .......... .......... 40% 78.7M 2s
    ##  54000K .......... .......... .......... .......... .......... 40% 4.96M 2s
    ##  54050K .......... .......... .......... .......... .......... 40% 61.0M 2s
    ##  54100K .......... .......... .......... .......... .......... 40% 57.6M 2s
    ##  54150K .......... .......... .......... .......... .......... 40% 60.0M 2s
    ##  54200K .......... .......... .......... .......... .......... 40% 80.0M 2s
    ##  54250K .......... .......... .......... .......... .......... 40% 75.7M 2s
    ##  54300K .......... .......... .......... .......... .......... 40% 78.0M 2s
    ##  54350K .......... .......... .......... .......... .......... 40% 79.3M 2s
    ##  54400K .......... .......... .......... .......... .......... 40% 87.2M 2s
    ##  54450K .......... .......... .......... .......... .......... 40% 80.8M 2s
    ##  54500K .......... .......... .......... .......... .......... 40% 77.7M 2s
    ##  54550K .......... .......... .......... .......... .......... 40% 38.6M 2s
    ##  54600K .......... .......... .......... .......... .......... 40% 57.8M 2s
    ##  54650K .......... .......... .......... .......... .......... 40% 52.0M 2s
    ##  54700K .......... .......... .......... .......... .......... 40% 64.1M 2s
    ##  54750K .......... .......... .......... .......... .......... 40% 64.5M 2s
    ##  54800K .......... .......... .......... .......... .......... 40% 68.5M 2s
    ##  54850K .......... .......... .......... .......... .......... 40% 61.3M 2s
    ##  54900K .......... .......... .......... .......... .......... 40% 70.9M 2s
    ##  54950K .......... .......... .......... .......... .......... 40% 62.6M 2s
    ##  55000K .......... .......... .......... .......... .......... 40% 68.9M 2s
    ##  55050K .......... .......... .......... .......... .......... 40% 71.7M 2s
    ##  55100K .......... .......... .......... .......... .......... 40% 78.0M 2s
    ##  55150K .......... .......... .......... .......... .......... 40% 62.6M 2s
    ##  55200K .......... .......... .......... .......... .......... 41% 65.6M 2s
    ##  55250K .......... .......... .......... .......... .......... 41% 41.0M 2s
    ##  55300K .......... .......... .......... .......... .......... 41% 70.1M 2s
    ##  55350K .......... .......... .......... .......... .......... 41% 86.9M 2s
    ##  55400K .......... .......... .......... .......... .......... 41% 81.9M 2s
    ##  55450K .......... .......... .......... .......... .......... 41% 72.6M 2s
    ##  55500K .......... .......... .......... .......... .......... 41% 73.5M 2s
    ##  55550K .......... .......... .......... .......... .......... 41% 67.6M 2s
    ##  55600K .......... .......... .......... .......... .......... 41% 51.5M 2s
    ##  55650K .......... .......... .......... .......... .......... 41% 63.2M 2s
    ##  55700K .......... .......... .......... .......... .......... 41% 62.1M 2s
    ##  55750K .......... .......... .......... .......... .......... 41% 69.3M 2s
    ##  55800K .......... .......... .......... .......... .......... 41% 61.2M 2s
    ##  55850K .......... .......... .......... .......... .......... 41% 58.4M 2s
    ##  55900K .......... .......... .......... .......... .......... 41% 72.5M 2s
    ##  55950K .......... .......... .......... .......... .......... 41% 67.1M 2s
    ##  56000K .......... .......... .......... .......... .......... 41% 69.5M 2s
    ##  56050K .......... .......... .......... .......... .......... 41% 64.5M 2s
    ##  56100K .......... .......... .......... .......... .......... 41% 71.7M 1s
    ##  56150K .......... .......... .......... .......... .......... 41% 72.4M 1s
    ##  56200K .......... .......... .......... .......... .......... 41% 67.4M 1s
    ##  56250K .......... .......... .......... .......... .......... 41% 71.8M 1s
    ##  56300K .......... .......... .......... .......... .......... 41% 88.5M 1s
    ##  56350K .......... .......... .......... .......... .......... 41% 65.1M 1s
    ##  56400K .......... .......... .......... .......... .......... 41% 65.5M 1s
    ##  56450K .......... .......... .......... .......... .......... 41% 71.4M 1s
    ##  56500K .......... .......... .......... .......... .......... 41% 63.2M 1s
    ##  56550K .......... .......... .......... .......... .......... 42% 70.5M 1s
    ##  56600K .......... .......... .......... .......... .......... 42% 70.4M 1s
    ##  56650K .......... .......... .......... .......... .......... 42% 64.0M 1s
    ##  56700K .......... .......... .......... .......... .......... 42% 63.5M 1s
    ##  56750K .......... .......... .......... .......... .......... 42% 78.8M 1s
    ##  56800K .......... .......... .......... .......... .......... 42% 61.7M 1s
    ##  56850K .......... .......... .......... .......... .......... 42% 71.7M 1s
    ##  56900K .......... .......... .......... .......... .......... 42% 88.9M 1s
    ##  56950K .......... .......... .......... .......... .......... 42% 73.9M 1s
    ##  57000K .......... .......... .......... .......... .......... 42% 78.5M 1s
    ##  57050K .......... .......... .......... .......... .......... 42% 83.7M 1s
    ##  57100K .......... .......... .......... .......... .......... 42% 46.4M 1s
    ##  57150K .......... .......... .......... .......... .......... 42% 67.4M 1s
    ##  57200K .......... .......... .......... .......... .......... 42% 60.0M 1s
    ##  57250K .......... .......... .......... .......... .......... 42% 66.4M 1s
    ##  57300K .......... .......... .......... .......... .......... 42% 71.3M 1s
    ##  57350K .......... .......... .......... .......... .......... 42% 69.9M 1s
    ##  57400K .......... .......... .......... .......... .......... 42% 76.8M 1s
    ##  57450K .......... .......... .......... .......... .......... 42% 63.9M 1s
    ##  57500K .......... .......... .......... .......... .......... 42% 70.2M 1s
    ##  57550K .......... .......... .......... .......... .......... 42% 67.5M 1s
    ##  57600K .......... .......... .......... .......... .......... 42% 61.5M 1s
    ##  57650K .......... .......... .......... .......... .......... 42% 72.7M 1s
    ##  57700K .......... .......... .......... .......... .......... 42% 71.4M 1s
    ##  57750K .......... .......... .......... .......... .......... 42% 69.2M 1s
    ##  57800K .......... .......... .......... .......... .......... 42% 63.7M 1s
    ##  57850K .......... .......... .......... .......... .......... 42% 70.3M 1s
    ##  57900K .......... .......... .......... .......... .......... 43% 59.3M 1s
    ##  57950K .......... .......... .......... .......... .......... 43% 71.1M 1s
    ##  58000K .......... .......... .......... .......... .......... 43% 64.7M 1s
    ##  58050K .......... .......... .......... .......... .......... 43% 65.3M 1s
    ##  58100K .......... .......... .......... .......... .......... 43% 76.3M 1s
    ##  58150K .......... .......... .......... .......... .......... 43% 69.5M 1s
    ##  58200K .......... .......... .......... .......... .......... 43% 62.8M 1s
    ##  58250K .......... .......... .......... .......... .......... 43% 73.7M 1s
    ##  58300K .......... .......... .......... .......... .......... 43% 73.6M 1s
    ##  58350K .......... .......... .......... .......... .......... 43% 86.6M 1s
    ##  58400K .......... .......... .......... .......... .......... 43% 70.0M 1s
    ##  58450K .......... .......... .......... .......... .......... 43% 90.2M 1s
    ##  58500K .......... .......... .......... .......... .......... 43% 83.8M 1s
    ##  58550K .......... .......... .......... .......... .......... 43% 5.22M 1s
    ##  58600K .......... .......... .......... .......... .......... 43% 70.3M 1s
    ##  58650K .......... .......... .......... .......... .......... 43% 85.6M 1s
    ##  58700K .......... .......... .......... .......... .......... 43% 67.3M 1s
    ##  58750K .......... .......... .......... .......... .......... 43% 83.6M 1s
    ##  58800K .......... .......... .......... .......... .......... 43% 74.5M 1s
    ##  58850K .......... .......... .......... .......... .......... 43% 69.2M 1s
    ##  58900K .......... .......... .......... .......... .......... 43% 68.5M 1s
    ##  58950K .......... .......... .......... .......... .......... 43% 64.0M 1s
    ##  59000K .......... .......... .......... .......... .......... 43% 66.4M 1s
    ##  59050K .......... .......... .......... .......... .......... 43% 75.7M 1s
    ##  59100K .......... .......... .......... .......... .......... 43% 59.7M 1s
    ##  59150K .......... .......... .......... .......... .......... 43% 76.0M 1s
    ##  59200K .......... .......... .......... .......... .......... 43% 68.0M 1s
    ##  59250K .......... .......... .......... .......... .......... 44% 79.0M 1s
    ##  59300K .......... .......... .......... .......... .......... 44% 80.1M 1s
    ##  59350K .......... .......... .......... .......... .......... 44%  103M 1s
    ##  59400K .......... .......... .......... .......... .......... 44% 77.3M 1s
    ##  59450K .......... .......... .......... .......... .......... 44% 84.8M 1s
    ##  59500K .......... .......... .......... .......... .......... 44% 81.4M 1s
    ##  59550K .......... .......... .......... .......... .......... 44% 86.0M 1s
    ##  59600K .......... .......... .......... .......... .......... 44% 71.1M 1s
    ##  59650K .......... .......... .......... .......... .......... 44% 90.1M 1s
    ##  59700K .......... .......... .......... .......... .......... 44% 68.8M 1s
    ##  59750K .......... .......... .......... .......... .......... 44% 66.7M 1s
    ##  59800K .......... .......... .......... .......... .......... 44% 54.3M 1s
    ##  59850K .......... .......... .......... .......... .......... 44% 72.1M 1s
    ##  59900K .......... .......... .......... .......... .......... 44% 71.0M 1s
    ##  59950K .......... .......... .......... .......... .......... 44% 98.7M 1s
    ##  60000K .......... .......... .......... .......... .......... 44% 82.3M 1s
    ##  60050K .......... .......... .......... .......... .......... 44% 74.6M 1s
    ##  60100K .......... .......... .......... .......... .......... 44% 87.4M 1s
    ##  60150K .......... .......... .......... .......... .......... 44% 64.9M 1s
    ##  60200K .......... .......... .......... .......... .......... 44% 44.3M 1s
    ##  60250K .......... .......... .......... .......... .......... 44% 85.3M 1s
    ##  60300K .......... .......... .......... .......... .......... 44% 87.7M 1s
    ##  60350K .......... .......... .......... .......... .......... 44% 38.8M 1s
    ##  60400K .......... .......... .......... .......... .......... 44% 75.2M 1s
    ##  60450K .......... .......... .......... .......... .......... 44% 97.9M 1s
    ##  60500K .......... .......... .......... .......... .......... 44% 72.1M 1s
    ##  60550K .......... .......... .......... .......... .......... 44% 63.0M 1s
    ##  60600K .......... .......... .......... .......... .......... 45% 73.4M 1s
    ##  60650K .......... .......... .......... .......... .......... 45% 69.8M 1s
    ##  60700K .......... .......... .......... .......... .......... 45% 70.6M 1s
    ##  60750K .......... .......... .......... .......... .......... 45% 65.6M 1s
    ##  60800K .......... .......... .......... .......... .......... 45% 85.0M 1s
    ##  60850K .......... .......... .......... .......... .......... 45% 87.2M 1s
    ##  60900K .......... .......... .......... .......... .......... 45% 76.1M 1s
    ##  60950K .......... .......... .......... .......... .......... 45% 79.1M 1s
    ##  61000K .......... .......... .......... .......... .......... 45% 46.0M 1s
    ##  61050K .......... .......... .......... .......... .......... 45% 64.1M 1s
    ##  61100K .......... .......... .......... .......... .......... 45% 51.9M 1s
    ##  61150K .......... .......... .......... .......... .......... 45% 58.8M 1s
    ##  61200K .......... .......... .......... .......... .......... 45% 90.2M 1s
    ##  61250K .......... .......... .......... .......... .......... 45% 48.8M 1s
    ##  61300K .......... .......... .......... .......... .......... 45% 61.7M 1s
    ##  61350K .......... .......... .......... .......... .......... 45% 77.8M 1s
    ##  61400K .......... .......... .......... .......... .......... 45% 70.4M 1s
    ##  61450K .......... .......... .......... .......... .......... 45% 76.3M 1s
    ##  61500K .......... .......... .......... .......... .......... 45% 60.9M 1s
    ##  61550K .......... .......... .......... .......... .......... 45% 80.3M 1s
    ##  61600K .......... .......... .......... .......... .......... 45% 30.9M 1s
    ##  61650K .......... .......... .......... .......... .......... 45% 84.2M 1s
    ##  61700K .......... .......... .......... .......... .......... 45% 44.0M 1s
    ##  61750K .......... .......... .......... .......... .......... 45% 68.9M 1s
    ##  61800K .......... .......... .......... .......... .......... 45% 59.1M 1s
    ##  61850K .......... .......... .......... .......... .......... 45% 84.6M 1s
    ##  61900K .......... .......... .......... .......... .......... 45% 71.2M 1s
    ##  61950K .......... .......... .......... .......... .......... 46% 83.3M 1s
    ##  62000K .......... .......... .......... .......... .......... 46% 70.0M 1s
    ##  62050K .......... .......... .......... .......... .......... 46% 50.5M 1s
    ##  62100K .......... .......... .......... .......... .......... 46% 68.0M 1s
    ##  62150K .......... .......... .......... .......... .......... 46% 66.9M 1s
    ##  62200K .......... .......... .......... .......... .......... 46% 66.1M 1s
    ##  62250K .......... .......... .......... .......... .......... 46% 79.6M 1s
    ##  62300K .......... .......... .......... .......... .......... 46% 49.8M 1s
    ##  62350K .......... .......... .......... .......... .......... 46% 75.6M 1s
    ##  62400K .......... .......... .......... .......... .......... 46% 72.5M 1s
    ##  62450K .......... .......... .......... .......... .......... 46% 75.2M 1s
    ##  62500K .......... .......... .......... .......... .......... 46% 71.1M 1s
    ##  62550K .......... .......... .......... .......... .......... 46% 44.2M 1s
    ##  62600K .......... .......... .......... .......... .......... 46% 39.8M 1s
    ##  62650K .......... .......... .......... .......... .......... 46% 66.5M 1s
    ##  62700K .......... .......... .......... .......... .......... 46% 57.6M 1s
    ##  62750K .......... .......... .......... .......... .......... 46% 73.6M 1s
    ##  62800K .......... .......... .......... .......... .......... 46% 66.8M 1s
    ##  62850K .......... .......... .......... .......... .......... 46% 81.9M 1s
    ##  62900K .......... .......... .......... .......... .......... 46% 33.4M 1s
    ##  62950K .......... .......... .......... .......... .......... 46% 46.3M 1s
    ##  63000K .......... .......... .......... .......... .......... 46% 69.6M 1s
    ##  63050K .......... .......... .......... .......... .......... 46% 15.6M 1s
    ##  63100K .......... .......... .......... .......... .......... 46% 8.64M 1s
    ##  63150K .......... .......... .......... .......... .......... 46% 88.4M 1s
    ##  63200K .......... .......... .......... .......... .......... 46% 46.0M 1s
    ##  63250K .......... .......... .......... .......... .......... 46% 56.2M 1s
    ##  63300K .......... .......... .......... .......... .......... 47% 53.4M 1s
    ##  63350K .......... .......... .......... .......... .......... 47% 70.8M 1s
    ##  63400K .......... .......... .......... .......... .......... 47% 66.9M 1s
    ##  63450K .......... .......... .......... .......... .......... 47% 71.4M 1s
    ##  63500K .......... .......... .......... .......... .......... 47% 66.9M 1s
    ##  63550K .......... .......... .......... .......... .......... 47% 69.4M 1s
    ##  63600K .......... .......... .......... .......... .......... 47% 83.9M 1s
    ##  63650K .......... .......... .......... .......... .......... 47%  102M 1s
    ##  63700K .......... .......... .......... .......... .......... 47% 35.7M 1s
    ##  63750K .......... .......... .......... .......... .......... 47% 89.7M 1s
    ##  63800K .......... .......... .......... .......... .......... 47% 30.7M 1s
    ##  63850K .......... .......... .......... .......... .......... 47% 66.2M 1s
    ##  63900K .......... .......... .......... .......... .......... 47% 50.8M 1s
    ##  63950K .......... .......... .......... .......... .......... 47% 72.3M 1s
    ##  64000K .......... .......... .......... .......... .......... 47% 65.4M 1s
    ##  64050K .......... .......... .......... .......... .......... 47% 69.9M 1s
    ##  64100K .......... .......... .......... .......... .......... 47% 62.6M 1s
    ##  64150K .......... .......... .......... .......... .......... 47% 68.2M 1s
    ##  64200K .......... .......... .......... .......... .......... 47% 64.0M 1s
    ##  64250K .......... .......... .......... .......... .......... 47% 64.9M 1s
    ##  64300K .......... .......... .......... .......... .......... 47% 55.6M 1s
    ##  64350K .......... .......... .......... .......... .......... 47% 74.1M 1s
    ##  64400K .......... .......... .......... .......... .......... 47% 57.8M 1s
    ##  64450K .......... .......... .......... .......... .......... 47% 68.5M 1s
    ##  64500K .......... .......... .......... .......... .......... 47% 63.0M 1s
    ##  64550K .......... .......... .......... .......... .......... 47% 68.8M 1s
    ##  64600K .......... .......... .......... .......... .......... 47% 76.0M 1s
    ##  64650K .......... .......... .......... .......... .......... 48% 73.0M 1s
    ##  64700K .......... .......... .......... .......... .......... 48% 75.6M 1s
    ##  64750K .......... .......... .......... .......... .......... 48% 67.0M 1s
    ##  64800K .......... .......... .......... .......... .......... 48% 76.2M 1s
    ##  64850K .......... .......... .......... .......... .......... 48% 79.3M 1s
    ##  64900K .......... .......... .......... .......... .......... 48% 71.6M 1s
    ##  64950K .......... .......... .......... .......... .......... 48% 7.96M 1s
    ##  65000K .......... .......... .......... .......... .......... 48% 72.7M 1s
    ##  65050K .......... .......... .......... .......... .......... 48% 76.1M 1s
    ##  65100K .......... .......... .......... .......... .......... 48% 76.7M 1s
    ##  65150K .......... .......... .......... .......... .......... 48% 96.3M 1s
    ##  65200K .......... .......... .......... .......... .......... 48% 82.9M 1s
    ##  65250K .......... .......... .......... .......... .......... 48% 64.5M 1s
    ##  65300K .......... .......... .......... .......... .......... 48% 71.8M 1s
    ##  65350K .......... .......... .......... .......... .......... 48% 61.1M 1s
    ##  65400K .......... .......... .......... .......... .......... 48% 70.5M 1s
    ##  65450K .......... .......... .......... .......... .......... 48% 81.2M 1s
    ##  65500K .......... .......... .......... .......... .......... 48% 75.5M 1s
    ##  65550K .......... .......... .......... .......... .......... 48% 81.3M 1s
    ##  65600K .......... .......... .......... .......... .......... 48% 87.7M 1s
    ##  65650K .......... .......... .......... .......... .......... 48% 79.1M 1s
    ##  65700K .......... .......... .......... .......... .......... 48% 67.4M 1s
    ##  65750K .......... .......... .......... .......... .......... 48% 77.1M 1s
    ##  65800K .......... .......... .......... .......... .......... 48% 61.5M 1s
    ##  65850K .......... .......... .......... .......... .......... 48% 64.1M 1s
    ##  65900K .......... .......... .......... .......... .......... 48% 68.7M 1s
    ##  65950K .......... .......... .......... .......... .......... 48% 69.2M 1s
    ##  66000K .......... .......... .......... .......... .......... 49% 58.2M 1s
    ##  66050K .......... .......... .......... .......... .......... 49% 75.6M 1s
    ##  66100K .......... .......... .......... .......... .......... 49% 72.0M 1s
    ##  66150K .......... .......... .......... .......... .......... 49% 67.3M 1s
    ##  66200K .......... .......... .......... .......... .......... 49% 66.8M 1s
    ##  66250K .......... .......... .......... .......... .......... 49% 77.9M 1s
    ##  66300K .......... .......... .......... .......... .......... 49% 68.8M 1s
    ##  66350K .......... .......... .......... .......... .......... 49% 70.4M 1s
    ##  66400K .......... .......... .......... .......... .......... 49% 63.8M 1s
    ##  66450K .......... .......... .......... .......... .......... 49% 72.6M 1s
    ##  66500K .......... .......... .......... .......... .......... 49% 67.1M 1s
    ##  66550K .......... .......... .......... .......... .......... 49% 78.8M 1s
    ##  66600K .......... .......... .......... .......... .......... 49% 69.5M 1s
    ##  66650K .......... .......... .......... .......... .......... 49% 81.8M 1s
    ##  66700K .......... .......... .......... .......... .......... 49% 81.5M 1s
    ##  66750K .......... .......... .......... .......... .......... 49% 70.5M 1s
    ##  66800K .......... .......... .......... .......... .......... 49% 70.9M 1s
    ##  66850K .......... .......... .......... .......... .......... 49% 15.0M 1s
    ##  66900K .......... .......... .......... .......... .......... 49% 32.2M 1s
    ##  66950K .......... .......... .......... .......... .......... 49% 41.4M 1s
    ##  67000K .......... .......... .......... .......... .......... 49% 61.1M 1s
    ##  67050K .......... .......... .......... .......... .......... 49% 80.3M 1s
    ##  67100K .......... .......... .......... .......... .......... 49% 87.9M 1s
    ##  67150K .......... .......... .......... .......... .......... 49% 6.60M 1s
    ##  67200K .......... .......... .......... .......... .......... 49% 38.1M 1s
    ##  67250K .......... .......... .......... .......... .......... 49% 48.4M 1s
    ##  67300K .......... .......... .......... .......... .......... 49% 54.0M 1s
    ##  67350K .......... .......... .......... .......... .......... 50% 71.7M 1s
    ##  67400K .......... .......... .......... .......... .......... 50% 70.2M 1s
    ##  67450K .......... .......... .......... .......... .......... 50% 30.0M 1s
    ##  67500K .......... .......... .......... .......... .......... 50% 62.1M 1s
    ##  67550K .......... .......... .......... .......... .......... 50% 66.7M 1s
    ##  67600K .......... .......... .......... .......... .......... 50% 64.9M 1s
    ##  67650K .......... .......... .......... .......... .......... 50% 75.9M 1s
    ##  67700K .......... .......... .......... .......... .......... 50% 75.7M 1s
    ##  67750K .......... .......... .......... .......... .......... 50% 16.5M 1s
    ##  67800K .......... .......... .......... .......... .......... 50% 59.5M 1s
    ##  67850K .......... .......... .......... .......... .......... 50% 92.3M 1s
    ##  67900K .......... .......... .......... .......... .......... 50% 83.5M 1s
    ##  67950K .......... .......... .......... .......... .......... 50%  102M 1s
    ##  68000K .......... .......... .......... .......... .......... 50% 68.5M 1s
    ##  68050K .......... .......... .......... .......... .......... 50% 24.1M 1s
    ##  68100K .......... .......... .......... .......... .......... 50% 86.4M 1s
    ##  68150K .......... .......... .......... .......... .......... 50% 47.1M 1s
    ##  68200K .......... .......... .......... .......... .......... 50% 59.1M 1s
    ##  68250K .......... .......... .......... .......... .......... 50% 78.6M 1s
    ##  68300K .......... .......... .......... .......... .......... 50% 8.73M 1s
    ##  68350K .......... .......... .......... .......... .......... 50% 63.7M 1s
    ##  68400K .......... .......... .......... .......... .......... 50% 59.9M 1s
    ##  68450K .......... .......... .......... .......... .......... 50% 69.5M 1s
    ##  68500K .......... .......... .......... .......... .......... 50% 77.1M 1s
    ##  68550K .......... .......... .......... .......... .......... 50%  101M 1s
    ##  68600K .......... .......... .......... .......... .......... 50% 13.9M 1s
    ##  68650K .......... .......... .......... .......... .......... 50% 38.5M 1s
    ##  68700K .......... .......... .......... .......... .......... 51% 59.2M 1s
    ##  68750K .......... .......... .......... .......... .......... 51% 66.9M 1s
    ##  68800K .......... .......... .......... .......... .......... 51% 70.6M 1s
    ##  68850K .......... .......... .......... .......... .......... 51% 17.5M 1s
    ##  68900K .......... .......... .......... .......... .......... 51% 76.9M 1s
    ##  68950K .......... .......... .......... .......... .......... 51% 43.1M 1s
    ##  69000K .......... .......... .......... .......... .......... 51% 66.3M 1s
    ##  69050K .......... .......... .......... .......... .......... 51%  105M 1s
    ##  69100K .......... .......... .......... .......... .......... 51% 72.3M 1s
    ##  69150K .......... .......... .......... .......... .......... 51% 16.3M 1s
    ##  69200K .......... .......... .......... .......... .......... 51% 64.2M 1s
    ##  69250K .......... .......... .......... .......... .......... 51% 62.8M 1s
    ##  69300K .......... .......... .......... .......... .......... 51% 42.9M 1s
    ##  69350K .......... .......... .......... .......... .......... 51% 4.97M 1s
    ##  69400K .......... .......... .......... .......... .......... 51% 73.5M 1s
    ##  69450K .......... .......... .......... .......... .......... 51% 42.4M 1s
    ##  69500K .......... .......... .......... .......... .......... 51% 59.9M 1s
    ##  69550K .......... .......... .......... .......... .......... 51% 85.5M 1s
    ##  69600K .......... .......... .......... .......... .......... 51% 64.8M 1s
    ##  69650K .......... .......... .......... .......... .......... 51% 85.7M 1s
    ##  69700K .......... .......... .......... .......... .......... 51% 17.3M 1s
    ##  69750K .......... .......... .......... .......... .......... 51% 54.9M 1s
    ##  69800K .......... .......... .......... .......... .......... 51% 82.4M 1s
    ##  69850K .......... .......... .......... .......... .......... 51% 96.8M 1s
    ##  69900K .......... .......... .......... .......... .......... 51% 26.4M 1s
    ##  69950K .......... .......... .......... .......... .......... 51% 63.9M 1s
    ##  70000K .......... .......... .......... .......... .......... 51% 62.9M 1s
    ##  70050K .......... .......... .......... .......... .......... 52% 84.3M 1s
    ##  70100K .......... .......... .......... .......... .......... 52% 43.2M 1s
    ##  70150K .......... .......... .......... .......... .......... 52% 67.1M 1s
    ##  70200K .......... .......... .......... .......... .......... 52% 32.7M 1s
    ##  70250K .......... .......... .......... .......... .......... 52% 95.3M 1s
    ##  70300K .......... .......... .......... .......... .......... 52% 4.24M 1s
    ##  70350K .......... .......... .......... .......... .......... 52% 86.4M 1s
    ##  70400K .......... .......... .......... .......... .......... 52% 77.9M 1s
    ##  70450K .......... .......... .......... .......... .......... 52% 84.9M 1s
    ##  70500K .......... .......... .......... .......... .......... 52% 85.1M 1s
    ##  70550K .......... .......... .......... .......... .......... 52% 86.1M 1s
    ##  70600K .......... .......... .......... .......... .......... 52% 12.6M 1s
    ##  70650K .......... .......... .......... .......... .......... 52% 52.1M 1s
    ##  70700K .......... .......... .......... .......... .......... 52% 64.1M 1s
    ##  70750K .......... .......... .......... .......... .......... 52% 86.2M 1s
    ##  70800K .......... .......... .......... .......... .......... 52% 12.8M 1s
    ##  70850K .......... .......... .......... .......... .......... 52% 76.4M 1s
    ##  70900K .......... .......... .......... .......... .......... 52% 83.8M 1s
    ##  70950K .......... .......... .......... .......... .......... 52% 83.4M 1s
    ##  71000K .......... .......... .......... .......... .......... 52% 89.2M 1s
    ##  71050K .......... .......... .......... .......... .......... 52% 81.4M 1s
    ##  71100K .......... .......... .......... .......... .......... 52% 25.3M 1s
    ##  71150K .......... .......... .......... .......... .......... 52% 68.8M 1s
    ##  71200K .......... .......... .......... .......... .......... 52% 30.6M 1s
    ##  71250K .......... .......... .......... .......... .......... 52% 65.9M 1s
    ##  71300K .......... .......... .......... .......... .......... 52% 79.9M 1s
    ##  71350K .......... .......... .......... .......... .......... 52% 96.7M 1s
    ##  71400K .......... .......... .......... .......... .......... 53% 55.5M 1s
    ##  71450K .......... .......... .......... .......... .......... 53% 73.6M 1s
    ##  71500K .......... .......... .......... .......... .......... 53% 13.8M 1s
    ##  71550K .......... .......... .......... .......... .......... 53% 53.6M 1s
    ##  71600K .......... .......... .......... .......... .......... 53% 61.7M 1s
    ##  71650K .......... .......... .......... .......... .......... 53% 93.8M 1s
    ##  71700K .......... .......... .......... .......... .......... 53% 22.3M 1s
    ##  71750K .......... .......... .......... .......... .......... 53% 10.4M 1s
    ##  71800K .......... .......... .......... .......... .......... 53% 68.9M 1s
    ##  71850K .......... .......... .......... .......... .......... 53% 72.1M 1s
    ##  71900K .......... .......... .......... .......... .......... 53% 72.8M 1s
    ##  71950K .......... .......... .......... .......... .......... 53% 83.0M 1s
    ##  72000K .......... .......... .......... .......... .......... 53% 77.9M 1s
    ##  72050K .......... .......... .......... .......... .......... 53% 92.2M 1s
    ##  72100K .......... .......... .......... .......... .......... 53% 15.9M 1s
    ##  72150K .......... .......... .......... .......... .......... 53% 33.9M 1s
    ##  72200K .......... .......... .......... .......... .......... 53% 63.6M 1s
    ##  72250K .......... .......... .......... .......... .......... 53% 62.2M 1s
    ##  72300K .......... .......... .......... .......... .......... 53% 69.5M 1s
    ##  72350K .......... .......... .......... .......... .......... 53% 81.1M 1s
    ##  72400K .......... .......... .......... .......... .......... 53% 76.4M 1s
    ##  72450K .......... .......... .......... .......... .......... 53% 19.0M 1s
    ##  72500K .......... .......... .......... .......... .......... 53% 75.3M 1s
    ##  72550K .......... .......... .......... .......... .......... 53% 42.4M 1s
    ##  72600K .......... .......... .......... .......... .......... 53% 75.3M 1s
    ##  72650K .......... .......... .......... .......... .......... 53% 71.5M 1s
    ##  72700K .......... .......... .......... .......... .......... 53% 77.0M 1s
    ##  72750K .......... .......... .......... .......... .......... 54% 90.3M 1s
    ##  72800K .......... .......... .......... .......... .......... 54% 24.8M 1s
    ##  72850K .......... .......... .......... .......... .......... 54% 56.7M 1s
    ##  72900K .......... .......... .......... .......... .......... 54% 70.0M 1s
    ##  72950K .......... .......... .......... .......... .......... 54% 60.7M 1s
    ##  73000K .......... .......... .......... .......... .......... 54% 7.46M 1s
    ##  73050K .......... .......... .......... .......... .......... 54% 84.2M 1s
    ##  73100K .......... .......... .......... .......... .......... 54% 47.6M 1s
    ##  73150K .......... .......... .......... .......... .......... 54% 45.2M 1s
    ##  73200K .......... .......... .......... .......... .......... 54% 70.7M 1s
    ##  73250K .......... .......... .......... .......... .......... 54% 72.9M 1s
    ##  73300K .......... .......... .......... .......... .......... 54% 62.8M 1s
    ##  73350K .......... .......... .......... .......... .......... 54% 80.9M 1s
    ##  73400K .......... .......... .......... .......... .......... 54% 69.4M 1s
    ##  73450K .......... .......... .......... .......... .......... 54% 58.3M 1s
    ##  73500K .......... .......... .......... .......... .......... 54% 34.8M 1s
    ##  73550K .......... .......... .......... .......... .......... 54% 42.3M 1s
    ##  73600K .......... .......... .......... .......... .......... 54% 41.7M 1s
    ##  73650K .......... .......... .......... .......... .......... 54% 43.5M 1s
    ##  73700K .......... .......... .......... .......... .......... 54% 38.8M 1s
    ##  73750K .......... .......... .......... .......... .......... 54% 35.6M 1s
    ##  73800K .......... .......... .......... .......... .......... 54% 42.4M 1s
    ##  73850K .......... .......... .......... .......... .......... 54% 41.7M 1s
    ##  73900K .......... .......... .......... .......... .......... 54% 35.3M 1s
    ##  73950K .......... .......... .......... .......... .......... 54% 50.3M 1s
    ##  74000K .......... .......... .......... .......... .......... 54% 43.4M 1s
    ##  74050K .......... .......... .......... .......... .......... 54% 30.9M 1s
    ##  74100K .......... .......... .......... .......... .......... 55% 47.0M 1s
    ##  74150K .......... .......... .......... .......... .......... 55% 49.1M 1s
    ##  74200K .......... .......... .......... .......... .......... 55% 43.4M 1s
    ##  74250K .......... .......... .......... .......... .......... 55% 55.5M 1s
    ##  74300K .......... .......... .......... .......... .......... 55% 48.4M 1s
    ##  74350K .......... .......... .......... .......... .......... 55% 45.0M 1s
    ##  74400K .......... .......... .......... .......... .......... 55% 52.6M 1s
    ##  74450K .......... .......... .......... .......... .......... 55% 50.6M 1s
    ##  74500K .......... .......... .......... .......... .......... 55% 11.6M 1s
    ##  74550K .......... .......... .......... .......... .......... 55% 39.4M 1s
    ##  74600K .......... .......... .......... .......... .......... 55% 36.2M 1s
    ##  74650K .......... .......... .......... .......... .......... 55% 41.9M 1s
    ##  74700K .......... .......... .......... .......... .......... 55% 41.6M 1s
    ##  74750K .......... .......... .......... .......... .......... 55% 40.6M 1s
    ##  74800K .......... .......... .......... .......... .......... 55% 43.3M 1s
    ##  74850K .......... .......... .......... .......... .......... 55% 52.2M 1s
    ##  74900K .......... .......... .......... .......... .......... 55% 46.4M 1s
    ##  74950K .......... .......... .......... .......... .......... 55% 41.8M 1s
    ##  75000K .......... .......... .......... .......... .......... 55% 51.9M 1s
    ##  75050K .......... .......... .......... .......... .......... 55% 54.6M 1s
    ##  75100K .......... .......... .......... .......... .......... 55% 49.1M 1s
    ##  75150K .......... .......... .......... .......... .......... 55% 57.2M 1s
    ##  75200K .......... .......... .......... .......... .......... 55% 56.4M 1s
    ##  75250K .......... .......... .......... .......... .......... 55% 11.3M 1s
    ##  75300K .......... .......... .......... .......... .......... 55% 42.7M 1s
    ##  75350K .......... .......... .......... .......... .......... 55% 47.7M 1s
    ##  75400K .......... .......... .......... .......... .......... 55% 60.0M 1s
    ##  75450K .......... .......... .......... .......... .......... 56% 68.4M 1s
    ##  75500K .......... .......... .......... .......... .......... 56% 17.4M 1s
    ##  75550K .......... .......... .......... .......... .......... 56% 32.6M 1s
    ##  75600K .......... .......... .......... .......... .......... 56% 46.7M 1s
    ##  75650K .......... .......... .......... .......... .......... 56% 63.7M 1s
    ##  75700K .......... .......... .......... .......... .......... 56% 51.7M 1s
    ##  75750K .......... .......... .......... .......... .......... 56% 49.2M 1s
    ##  75800K .......... .......... .......... .......... .......... 56% 50.3M 1s
    ##  75850K .......... .......... .......... .......... .......... 56% 46.9M 1s
    ##  75900K .......... .......... .......... .......... .......... 56% 51.8M 1s
    ##  75950K .......... .......... .......... .......... .......... 56% 20.2M 1s
    ##  76000K .......... .......... .......... .......... .......... 56% 40.4M 1s
    ##  76050K .......... .......... .......... .......... .......... 56% 55.9M 1s
    ##  76100K .......... .......... .......... .......... .......... 56% 47.3M 1s
    ##  76150K .......... .......... .......... .......... .......... 56% 53.2M 1s
    ##  76200K .......... .......... .......... .......... .......... 56% 52.9M 1s
    ##  76250K .......... .......... .......... .......... .......... 56% 65.5M 1s
    ##  76300K .......... .......... .......... .......... .......... 56% 45.7M 1s
    ##  76350K .......... .......... .......... .......... .......... 56% 74.5M 1s
    ##  76400K .......... .......... .......... .......... .......... 56% 51.9M 1s
    ##  76450K .......... .......... .......... .......... .......... 56% 44.4M 1s
    ##  76500K .......... .......... .......... .......... .......... 56% 46.4M 1s
    ##  76550K .......... .......... .......... .......... .......... 56% 52.0M 1s
    ##  76600K .......... .......... .......... .......... .......... 56% 71.4M 1s
    ##  76650K .......... .......... .......... .......... .......... 56% 62.5M 1s
    ##  76700K .......... .......... .......... .......... .......... 56% 65.8M 1s
    ##  76750K .......... .......... .......... .......... .......... 56% 63.0M 1s
    ##  76800K .......... .......... .......... .......... .......... 57% 56.6M 1s
    ##  76850K .......... .......... .......... .......... .......... 57% 6.58M 1s
    ##  76900K .......... .......... .......... .......... .......... 57% 60.7M 1s
    ##  76950K .......... .......... .......... .......... .......... 57% 41.8M 1s
    ##  77000K .......... .......... .......... .......... .......... 57% 55.7M 1s
    ##  77050K .......... .......... .......... .......... .......... 57% 83.5M 1s
    ##  77100K .......... .......... .......... .......... .......... 57% 59.8M 1s
    ##  77150K .......... .......... .......... .......... .......... 57% 63.6M 1s
    ##  77200K .......... .......... .......... .......... .......... 57% 64.9M 1s
    ##  77250K .......... .......... .......... .......... .......... 57% 51.8M 1s
    ##  77300K .......... .......... .......... .......... .......... 57% 59.2M 1s
    ##  77350K .......... .......... .......... .......... .......... 57% 56.7M 1s
    ##  77400K .......... .......... .......... .......... .......... 57% 58.4M 1s
    ##  77450K .......... .......... .......... .......... .......... 57% 70.5M 1s
    ##  77500K .......... .......... .......... .......... .......... 57% 38.7M 1s
    ##  77550K .......... .......... .......... .......... .......... 57% 70.1M 1s
    ##  77600K .......... .......... .......... .......... .......... 57% 48.5M 1s
    ##  77650K .......... .......... .......... .......... .......... 57% 69.4M 1s
    ##  77700K .......... .......... .......... .......... .......... 57% 6.71M 1s
    ##  77750K .......... .......... .......... .......... .......... 57% 37.9M 1s
    ##  77800K .......... .......... .......... .......... .......... 57% 50.3M 1s
    ##  77850K .......... .......... .......... .......... .......... 57% 57.4M 1s
    ##  77900K .......... .......... .......... .......... .......... 57% 54.0M 1s
    ##  77950K .......... .......... .......... .......... .......... 57% 58.5M 1s
    ##  78000K .......... .......... .......... .......... .......... 57% 12.7M 1s
    ##  78050K .......... .......... .......... .......... .......... 57% 25.4M 1s
    ##  78100K .......... .......... .......... .......... .......... 58% 65.1M 1s
    ##  78150K .......... .......... .......... .......... .......... 58% 59.9M 1s
    ##  78200K .......... .......... .......... .......... .......... 58% 45.2M 1s
    ##  78250K .......... .......... .......... .......... .......... 58% 56.7M 1s
    ##  78300K .......... .......... .......... .......... .......... 58% 68.1M 1s
    ##  78350K .......... .......... .......... .......... .......... 58% 45.9M 1s
    ##  78400K .......... .......... .......... .......... .......... 58% 35.5M 1s
    ##  78450K .......... .......... .......... .......... .......... 58% 36.1M 1s
    ##  78500K .......... .......... .......... .......... .......... 58% 36.2M 1s
    ##  78550K .......... .......... .......... .......... .......... 58% 39.7M 1s
    ##  78600K .......... .......... .......... .......... .......... 58% 39.9M 1s
    ##  78650K .......... .......... .......... .......... .......... 58% 36.3M 1s
    ##  78700K .......... .......... .......... .......... .......... 58% 38.5M 1s
    ##  78750K .......... .......... .......... .......... .......... 58% 41.9M 1s
    ##  78800K .......... .......... .......... .......... .......... 58% 42.7M 1s
    ##  78850K .......... .......... .......... .......... .......... 58% 41.9M 1s
    ##  78900K .......... .......... .......... .......... .......... 58% 35.5M 1s
    ##  78950K .......... .......... .......... .......... .......... 58% 49.6M 1s
    ##  79000K .......... .......... .......... .......... .......... 58% 36.6M 1s
    ##  79050K .......... .......... .......... .......... .......... 58% 43.9M 1s
    ##  79100K .......... .......... .......... .......... .......... 58% 40.8M 1s
    ##  79150K .......... .......... .......... .......... .......... 58% 41.1M 1s
    ##  79200K .......... .......... .......... .......... .......... 58% 42.1M 1s
    ##  79250K .......... .......... .......... .......... .......... 58% 44.0M 1s
    ##  79300K .......... .......... .......... .......... .......... 58% 33.3M 1s
    ##  79350K .......... .......... .......... .......... .......... 58% 41.8M 1s
    ##  79400K .......... .......... .......... .......... .......... 58% 40.4M 1s
    ##  79450K .......... .......... .......... .......... .......... 59% 41.6M 1s
    ##  79500K .......... .......... .......... .......... .......... 59% 38.4M 1s
    ##  79550K .......... .......... .......... .......... .......... 59% 42.2M 1s
    ##  79600K .......... .......... .......... .......... .......... 59% 39.1M 1s
    ##  79650K .......... .......... .......... .......... .......... 59% 41.5M 1s
    ##  79700K .......... .......... .......... .......... .......... 59% 41.8M 1s
    ##  79750K .......... .......... .......... .......... .......... 59% 38.8M 1s
    ##  79800K .......... .......... .......... .......... .......... 59% 42.1M 1s
    ##  79850K .......... .......... .......... .......... .......... 59% 55.8M 1s
    ##  79900K .......... .......... .......... .......... .......... 59% 43.3M 1s
    ##  79950K .......... .......... .......... .......... .......... 59% 51.8M 1s
    ##  80000K .......... .......... .......... .......... .......... 59% 44.8M 1s
    ##  80050K .......... .......... .......... .......... .......... 59% 46.2M 1s
    ##  80100K .......... .......... .......... .......... .......... 59% 37.2M 1s
    ##  80150K .......... .......... .......... .......... .......... 59% 51.1M 1s
    ##  80200K .......... .......... .......... .......... .......... 59% 48.5M 1s
    ##  80250K .......... .......... .......... .......... .......... 59% 53.1M 1s
    ##  80300K .......... .......... .......... .......... .......... 59% 41.6M 1s
    ##  80350K .......... .......... .......... .......... .......... 59% 37.2M 1s
    ##  80400K .......... .......... .......... .......... .......... 59% 42.0M 1s
    ##  80450K .......... .......... .......... .......... .......... 59% 49.4M 1s
    ##  80500K .......... .......... .......... .......... .......... 59% 47.2M 1s
    ##  80550K .......... .......... .......... .......... .......... 59% 48.7M 1s
    ##  80600K .......... .......... .......... .......... .......... 59% 44.6M 1s
    ##  80650K .......... .......... .......... .......... .......... 59% 42.4M 1s
    ##  80700K .......... .......... .......... .......... .......... 59% 46.3M 1s
    ##  80750K .......... .......... .......... .......... .......... 59% 68.5M 1s
    ##  80800K .......... .......... .......... .......... .......... 60% 62.3M 1s
    ##  80850K .......... .......... .......... .......... .......... 60% 59.4M 1s
    ##  80900K .......... .......... .......... .......... .......... 60% 50.5M 1s
    ##  80950K .......... .......... .......... .......... .......... 60% 62.8M 1s
    ##  81000K .......... .......... .......... .......... .......... 60% 49.2M 1s
    ##  81050K .......... .......... .......... .......... .......... 60% 60.0M 1s
    ##  81100K .......... .......... .......... .......... .......... 60% 48.6M 1s
    ##  81150K .......... .......... .......... .......... .......... 60% 48.6M 1s
    ##  81200K .......... .......... .......... .......... .......... 60% 51.3M 1s
    ##  81250K .......... .......... .......... .......... .......... 60% 47.3M 1s
    ##  81300K .......... .......... .......... .......... .......... 60% 45.2M 1s
    ##  81350K .......... .......... .......... .......... .......... 60% 46.6M 1s
    ##  81400K .......... .......... .......... .......... .......... 60% 54.0M 1s
    ##  81450K .......... .......... .......... .......... .......... 60% 53.2M 1s
    ##  81500K .......... .......... .......... .......... .......... 60% 55.7M 1s
    ##  81550K .......... .......... .......... .......... .......... 60% 60.8M 1s
    ##  81600K .......... .......... .......... .......... .......... 60% 43.8M 1s
    ##  81650K .......... .......... .......... .......... .......... 60% 53.4M 1s
    ##  81700K .......... .......... .......... .......... .......... 60% 72.1M 1s
    ##  81750K .......... .......... .......... .......... .......... 60% 61.4M 1s
    ##  81800K .......... .......... .......... .......... .......... 60% 49.1M 1s
    ##  81850K .......... .......... .......... .......... .......... 60% 63.6M 1s
    ##  81900K .......... .......... .......... .......... .......... 60% 45.3M 1s
    ##  81950K .......... .......... .......... .......... .......... 60% 63.6M 1s
    ##  82000K .......... .......... .......... .......... .......... 60% 68.0M 1s
    ##  82050K .......... .......... .......... .......... .......... 60% 62.9M 1s
    ##  82100K .......... .......... .......... .......... .......... 60% 57.5M 1s
    ##  82150K .......... .......... .......... .......... .......... 61% 54.0M 1s
    ##  82200K .......... .......... .......... .......... .......... 61% 54.5M 1s
    ##  82250K .......... .......... .......... .......... .......... 61% 70.3M 1s
    ##  82300K .......... .......... .......... .......... .......... 61% 75.6M 1s
    ##  82350K .......... .......... .......... .......... .......... 61% 78.3M 1s
    ##  82400K .......... .......... .......... .......... .......... 61% 72.5M 1s
    ##  82450K .......... .......... .......... .......... .......... 61% 72.0M 1s
    ##  82500K .......... .......... .......... .......... .......... 61% 69.5M 1s
    ##  82550K .......... .......... .......... .......... .......... 61% 67.1M 1s
    ##  82600K .......... .......... .......... .......... .......... 61% 74.0M 1s
    ##  82650K .......... .......... .......... .......... .......... 61% 82.6M 1s
    ##  82700K .......... .......... .......... .......... .......... 61% 68.2M 1s
    ##  82750K .......... .......... .......... .......... .......... 61% 66.4M 1s
    ##  82800K .......... .......... .......... .......... .......... 61% 64.5M 1s
    ##  82850K .......... .......... .......... .......... .......... 61% 66.7M 1s
    ##  82900K .......... .......... .......... .......... .......... 61% 61.7M 1s
    ##  82950K .......... .......... .......... .......... .......... 61% 80.5M 1s
    ##  83000K .......... .......... .......... .......... .......... 61% 45.6M 1s
    ##  83050K .......... .......... .......... .......... .......... 61% 58.3M 1s
    ##  83100K .......... .......... .......... .......... .......... 61% 41.4M 1s
    ##  83150K .......... .......... .......... .......... .......... 61% 71.3M 1s
    ##  83200K .......... .......... .......... .......... .......... 61% 66.4M 1s
    ##  83250K .......... .......... .......... .......... .......... 61% 67.2M 1s
    ##  83300K .......... .......... .......... .......... .......... 61% 55.0M 1s
    ##  83350K .......... .......... .......... .......... .......... 61% 64.3M 1s
    ##  83400K .......... .......... .......... .......... .......... 61% 65.3M 1s
    ##  83450K .......... .......... .......... .......... .......... 61% 69.1M 1s
    ##  83500K .......... .......... .......... .......... .......... 62% 69.4M 1s
    ##  83550K .......... .......... .......... .......... .......... 62% 78.0M 1s
    ##  83600K .......... .......... .......... .......... .......... 62% 64.7M 1s
    ##  83650K .......... .......... .......... .......... .......... 62% 52.8M 1s
    ##  83700K .......... .......... .......... .......... .......... 62% 66.0M 1s
    ##  83750K .......... .......... .......... .......... .......... 62% 69.4M 1s
    ##  83800K .......... .......... .......... .......... .......... 62% 66.3M 1s
    ##  83850K .......... .......... .......... .......... .......... 62% 59.3M 1s
    ##  83900K .......... .......... .......... .......... .......... 62% 67.1M 1s
    ##  83950K .......... .......... .......... .......... .......... 62% 57.8M 1s
    ##  84000K .......... .......... .......... .......... .......... 62% 69.3M 1s
    ##  84050K .......... .......... .......... .......... .......... 62% 78.8M 1s
    ##  84100K .......... .......... .......... .......... .......... 62% 72.9M 1s
    ##  84150K .......... .......... .......... .......... .......... 62% 85.7M 1s
    ##  84200K .......... .......... .......... .......... .......... 62% 44.3M 1s
    ##  84250K .......... .......... .......... .......... .......... 62% 73.5M 1s
    ##  84300K .......... .......... .......... .......... .......... 62% 82.3M 1s
    ##  84350K .......... .......... .......... .......... .......... 62% 82.2M 1s
    ##  84400K .......... .......... .......... .......... .......... 62% 65.9M 1s
    ##  84450K .......... .......... .......... .......... .......... 62% 79.8M 1s
    ##  84500K .......... .......... .......... .......... .......... 62% 64.1M 1s
    ##  84550K .......... .......... .......... .......... .......... 62% 65.8M 1s
    ##  84600K .......... .......... .......... .......... .......... 62% 71.9M 1s
    ##  84650K .......... .......... .......... .......... .......... 62% 20.8M 1s
    ##  84700K .......... .......... .......... .......... .......... 62% 73.7M 1s
    ##  84750K .......... .......... .......... .......... .......... 62% 64.4M 1s
    ##  84800K .......... .......... .......... .......... .......... 62% 67.3M 1s
    ##  84850K .......... .......... .......... .......... .......... 63% 68.3M 1s
    ##  84900K .......... .......... .......... .......... .......... 63% 79.9M 1s
    ##  84950K .......... .......... .......... .......... .......... 63% 31.6M 1s
    ##  85000K .......... .......... .......... .......... .......... 63% 30.6M 1s
    ##  85050K .......... .......... .......... .......... .......... 63% 49.3M 1s
    ##  85100K .......... .......... .......... .......... .......... 63% 81.3M 1s
    ##  85150K .......... .......... .......... .......... .......... 63% 88.6M 1s
    ##  85200K .......... .......... .......... .......... .......... 63% 88.0M 1s
    ##  85250K .......... .......... .......... .......... .......... 63% 91.8M 1s
    ##  85300K .......... .......... .......... .......... .......... 63% 64.5M 1s
    ##  85350K .......... .......... .......... .......... .......... 63% 66.8M 1s
    ##  85400K .......... .......... .......... .......... .......... 63% 49.9M 1s
    ##  85450K .......... .......... .......... .......... .......... 63% 52.3M 1s
    ##  85500K .......... .......... .......... .......... .......... 63% 56.1M 1s
    ##  85550K .......... .......... .......... .......... .......... 63% 65.2M 1s
    ##  85600K .......... .......... .......... .......... .......... 63% 52.6M 1s
    ##  85650K .......... .......... .......... .......... .......... 63% 63.4M 1s
    ##  85700K .......... .......... .......... .......... .......... 63% 56.1M 1s
    ##  85750K .......... .......... .......... .......... .......... 63% 59.8M 1s
    ##  85800K .......... .......... .......... .......... .......... 63% 61.7M 1s
    ##  85850K .......... .......... .......... .......... .......... 63% 77.7M 1s
    ##  85900K .......... .......... .......... .......... .......... 63% 60.6M 1s
    ##  85950K .......... .......... .......... .......... .......... 63% 67.9M 1s
    ##  86000K .......... .......... .......... .......... .......... 63% 86.7M 1s
    ##  86050K .......... .......... .......... .......... .......... 63% 58.9M 1s
    ##  86100K .......... .......... .......... .......... .......... 63% 41.7M 1s
    ##  86150K .......... .......... .......... .......... .......... 63% 50.5M 1s
    ##  86200K .......... .......... .......... .......... .......... 64% 68.8M 1s
    ##  86250K .......... .......... .......... .......... .......... 64% 40.7M 1s
    ##  86300K .......... .......... .......... .......... .......... 64% 62.6M 1s
    ##  86350K .......... .......... .......... .......... .......... 64% 78.6M 1s
    ##  86400K .......... .......... .......... .......... .......... 64% 72.7M 1s
    ##  86450K .......... .......... .......... .......... .......... 64% 17.0M 1s
    ##  86500K .......... .......... .......... .......... .......... 64% 62.0M 1s
    ##  86550K .......... .......... .......... .......... .......... 64% 54.8M 1s
    ##  86600K .......... .......... .......... .......... .......... 64% 6.67M 1s
    ##  86650K .......... .......... .......... .......... .......... 64% 69.8M 1s
    ##  86700K .......... .......... .......... .......... .......... 64% 66.9M 1s
    ##  86750K .......... .......... .......... .......... .......... 64% 64.4M 1s
    ##  86800K .......... .......... .......... .......... .......... 64% 63.9M 1s
    ##  86850K .......... .......... .......... .......... .......... 64% 69.1M 1s
    ##  86900K .......... .......... .......... .......... .......... 64% 82.3M 1s
    ##  86950K .......... .......... .......... .......... .......... 64% 78.4M 1s
    ##  87000K .......... .......... .......... .......... .......... 64% 70.3M 1s
    ##  87050K .......... .......... .......... .......... .......... 64% 78.8M 1s
    ##  87100K .......... .......... .......... .......... .......... 64% 71.1M 1s
    ##  87150K .......... .......... .......... .......... .......... 64% 70.1M 1s
    ##  87200K .......... .......... .......... .......... .......... 64% 76.9M 1s
    ##  87250K .......... .......... .......... .......... .......... 64% 81.8M 1s
    ##  87300K .......... .......... .......... .......... .......... 64% 66.5M 1s
    ##  87350K .......... .......... .......... .......... .......... 64% 64.8M 1s
    ##  87400K .......... .......... .......... .......... .......... 64% 71.7M 1s
    ##  87450K .......... .......... .......... .......... .......... 64% 69.5M 1s
    ##  87500K .......... .......... .......... .......... .......... 64% 74.4M 1s
    ##  87550K .......... .......... .......... .......... .......... 65% 19.1M 1s
    ##  87600K .......... .......... .......... .......... .......... 65% 53.8M 1s
    ##  87650K .......... .......... .......... .......... .......... 65% 51.5M 1s
    ##  87700K .......... .......... .......... .......... .......... 65% 51.6M 1s
    ##  87750K .......... .......... .......... .......... .......... 65% 63.4M 1s
    ##  87800K .......... .......... .......... .......... .......... 65% 31.9M 1s
    ##  87850K .......... .......... .......... .......... .......... 65% 67.3M 1s
    ##  87900K .......... .......... .......... .......... .......... 65% 45.4M 1s
    ##  87950K .......... .......... .......... .......... .......... 65% 64.2M 1s
    ##  88000K .......... .......... .......... .......... .......... 65% 6.81M 1s
    ##  88050K .......... .......... .......... .......... .......... 65% 42.7M 1s
    ##  88100K .......... .......... .......... .......... .......... 65% 39.5M 1s
    ##  88150K .......... .......... .......... .......... .......... 65% 74.8M 1s
    ##  88200K .......... .......... .......... .......... .......... 65% 41.8M 1s
    ##  88250K .......... .......... .......... .......... .......... 65% 11.0M 1s
    ##  88300K .......... .......... .......... .......... .......... 65% 67.5M 1s
    ##  88350K .......... .......... .......... .......... .......... 65% 85.8M 1s
    ##  88400K .......... .......... .......... .......... .......... 65% 72.7M 1s
    ##  88450K .......... .......... .......... .......... .......... 65% 63.7M 1s
    ##  88500K .......... .......... .......... .......... .......... 65% 42.1M 1s
    ##  88550K .......... .......... .......... .......... .......... 65% 74.5M 1s
    ##  88600K .......... .......... .......... .......... .......... 65% 17.3M 1s
    ##  88650K .......... .......... .......... .......... .......... 65% 74.9M 1s
    ##  88700K .......... .......... .......... .......... .......... 65% 68.7M 1s
    ##  88750K .......... .......... .......... .......... .......... 65% 67.6M 1s
    ##  88800K .......... .......... .......... .......... .......... 65% 80.3M 1s
    ##  88850K .......... .......... .......... .......... .......... 65% 18.0M 1s
    ##  88900K .......... .......... .......... .......... .......... 66% 64.3M 1s
    ##  88950K .......... .......... .......... .......... .......... 66% 80.2M 1s
    ##  89000K .......... .......... .......... .......... .......... 66% 79.6M 1s
    ##  89050K .......... .......... .......... .......... .......... 66% 15.4M 1s
    ##  89100K .......... .......... .......... .......... .......... 66% 67.5M 1s
    ##  89150K .......... .......... .......... .......... .......... 66% 56.1M 1s
    ##  89200K .......... .......... .......... .......... .......... 66% 53.7M 1s
    ##  89250K .......... .......... .......... .......... .......... 66% 74.7M 1s
    ##  89300K .......... .......... .......... .......... .......... 66% 57.8M 1s
    ##  89350K .......... .......... .......... .......... .......... 66% 44.1M 1s
    ##  89400K .......... .......... .......... .......... .......... 66% 54.0M 1s
    ##  89450K .......... .......... .......... .......... .......... 66% 33.5M 1s
    ##  89500K .......... .......... .......... .......... .......... 66% 61.9M 1s
    ##  89550K .......... .......... .......... .......... .......... 66% 73.2M 1s
    ##  89600K .......... .......... .......... .......... .......... 66% 71.4M 1s
    ##  89650K .......... .......... .......... .......... .......... 66% 32.6M 1s
    ##  89700K .......... .......... .......... .......... .......... 66% 54.6M 1s
    ##  89750K .......... .......... .......... .......... .......... 66% 76.2M 1s
    ##  89800K .......... .......... .......... .......... .......... 66% 35.1M 1s
    ##  89850K .......... .......... .......... .......... .......... 66% 17.4M 1s
    ##  89900K .......... .......... .......... .......... .......... 66% 46.4M 1s
    ##  89950K .......... .......... .......... .......... .......... 66% 55.2M 1s
    ##  90000K .......... .......... .......... .......... .......... 66% 48.2M 1s
    ##  90050K .......... .......... .......... .......... .......... 66% 74.8M 1s
    ##  90100K .......... .......... .......... .......... .......... 66% 39.5M 1s
    ##  90150K .......... .......... .......... .......... .......... 66% 20.0M 1s
    ##  90200K .......... .......... .......... .......... .......... 66% 51.7M 1s
    ##  90250K .......... .......... .......... .......... .......... 67% 60.3M 1s
    ##  90300K .......... .......... .......... .......... .......... 67% 88.2M 1s
    ##  90350K .......... .......... .......... .......... .......... 67% 23.2M 1s
    ##  90400K .......... .......... .......... .......... .......... 67% 65.5M 1s
    ##  90450K .......... .......... .......... .......... .......... 67% 68.4M 1s
    ##  90500K .......... .......... .......... .......... .......... 67% 83.8M 1s
    ##  90550K .......... .......... .......... .......... .......... 67% 9.16M 1s
    ##  90600K .......... .......... .......... .......... .......... 67% 60.8M 1s
    ##  90650K .......... .......... .......... .......... .......... 67% 65.4M 1s
    ##  90700K .......... .......... .......... .......... .......... 67% 78.1M 1s
    ##  90750K .......... .......... .......... .......... .......... 67% 42.6M 1s
    ##  90800K .......... .......... .......... .......... .......... 67% 66.2M 1s
    ##  90850K .......... .......... .......... .......... .......... 67% 56.0M 1s
    ##  90900K .......... .......... .......... .......... .......... 67% 70.5M 1s
    ##  90950K .......... .......... .......... .......... .......... 67% 70.0M 1s
    ##  91000K .......... .......... .......... .......... .......... 67% 25.9M 1s
    ##  91050K .......... .......... .......... .......... .......... 67% 73.6M 1s
    ##  91100K .......... .......... .......... .......... .......... 67% 76.4M 1s
    ##  91150K .......... .......... .......... .......... .......... 67% 30.0M 1s
    ##  91200K .......... .......... .......... .......... .......... 67% 57.8M 1s
    ##  91250K .......... .......... .......... .......... .......... 67% 82.7M 1s
    ##  91300K .......... .......... .......... .......... .......... 67% 52.1M 1s
    ##  91350K .......... .......... .......... .......... .......... 67% 8.74M 1s
    ##  91400K .......... .......... .......... .......... .......... 67% 34.4M 1s
    ##  91450K .......... .......... .......... .......... .......... 67% 55.0M 1s
    ##  91500K .......... .......... .......... .......... .......... 67% 74.2M 1s
    ##  91550K .......... .......... .......... .......... .......... 67% 61.3M 1s
    ##  91600K .......... .......... .......... .......... .......... 68% 43.1M 1s
    ##  91650K .......... .......... .......... .......... .......... 68% 69.6M 1s
    ##  91700K .......... .......... .......... .......... .......... 68% 61.8M 1s
    ##  91750K .......... .......... .......... .......... .......... 68% 45.0M 1s
    ##  91800K .......... .......... .......... .......... .......... 68% 17.0M 1s
    ##  91850K .......... .......... .......... .......... .......... 68% 87.7M 1s
    ##  91900K .......... .......... .......... .......... .......... 68% 85.0M 1s
    ##  91950K .......... .......... .......... .......... .......... 68% 81.9M 1s
    ##  92000K .......... .......... .......... .......... .......... 68% 56.1M 1s
    ##  92050K .......... .......... .......... .......... .......... 68% 40.3M 1s
    ##  92100K .......... .......... .......... .......... .......... 68% 73.4M 1s
    ##  92150K .......... .......... .......... .......... .......... 68% 45.1M 1s
    ##  92200K .......... .......... .......... .......... .......... 68% 48.4M 1s
    ##  92250K .......... .......... .......... .......... .......... 68% 45.0M 1s
    ##  92300K .......... .......... .......... .......... .......... 68% 69.9M 1s
    ##  92350K .......... .......... .......... .......... .......... 68% 34.2M 1s
    ##  92400K .......... .......... .......... .......... .......... 68% 43.2M 1s
    ##  92450K .......... .......... .......... .......... .......... 68% 10.7M 1s
    ##  92500K .......... .......... .......... .......... .......... 68% 55.4M 1s
    ##  92550K .......... .......... .......... .......... .......... 68% 64.5M 1s
    ##  92600K .......... .......... .......... .......... .......... 68% 85.7M 1s
    ##  92650K .......... .......... .......... .......... .......... 68% 26.4M 1s
    ##  92700K .......... .......... .......... .......... .......... 68% 11.4M 1s
    ##  92750K .......... .......... .......... .......... .......... 68% 30.8M 1s
    ##  92800K .......... .......... .......... .......... .......... 68% 74.8M 1s
    ##  92850K .......... .......... .......... .......... .......... 68% 22.1M 1s
    ##  92900K .......... .......... .......... .......... .......... 68% 59.1M 1s
    ##  92950K .......... .......... .......... .......... .......... 69% 55.3M 1s
    ##  93000K .......... .......... .......... .......... .......... 69% 51.7M 1s
    ##  93050K .......... .......... .......... .......... .......... 69% 67.6M 1s
    ##  93100K .......... .......... .......... .......... .......... 69% 61.5M 1s
    ##  93150K .......... .......... .......... .......... .......... 69% 47.3M 1s
    ##  93200K .......... .......... .......... .......... .......... 69% 49.0M 1s
    ##  93250K .......... .......... .......... .......... .......... 69% 46.5M 1s
    ##  93300K .......... .......... .......... .......... .......... 69% 69.7M 1s
    ##  93350K .......... .......... .......... .......... .......... 69% 22.4M 1s
    ##  93400K .......... .......... .......... .......... .......... 69% 48.4M 1s
    ##  93450K .......... .......... .......... .......... .......... 69% 59.9M 1s
    ##  93500K .......... .......... .......... .......... .......... 69% 59.6M 1s
    ##  93550K .......... .......... .......... .......... .......... 69% 39.0M 1s
    ##  93600K .......... .......... .......... .......... .......... 69% 50.7M 1s
    ##  93650K .......... .......... .......... .......... .......... 69% 45.7M 1s
    ##  93700K .......... .......... .......... .......... .......... 69% 50.9M 1s
    ##  93750K .......... .......... .......... .......... .......... 69% 54.2M 1s
    ##  93800K .......... .......... .......... .......... .......... 69% 17.0M 1s
    ##  93850K .......... .......... .......... .......... .......... 69% 64.5M 1s
    ##  93900K .......... .......... .......... .......... .......... 69% 64.5M 1s
    ##  93950K .......... .......... .......... .......... .......... 69% 24.2M 1s
    ##  94000K .......... .......... .......... .......... .......... 69% 36.7M 1s
    ##  94050K .......... .......... .......... .......... .......... 69% 52.8M 1s
    ##  94100K .......... .......... .......... .......... .......... 69% 74.7M 1s
    ##  94150K .......... .......... .......... .......... .......... 69% 73.8M 1s
    ##  94200K .......... .......... .......... .......... .......... 69% 66.2M 1s
    ##  94250K .......... .......... .......... .......... .......... 69% 54.1M 1s
    ##  94300K .......... .......... .......... .......... .......... 70% 59.6M 1s
    ##  94350K .......... .......... .......... .......... .......... 70% 24.2M 1s
    ##  94400K .......... .......... .......... .......... .......... 70% 65.8M 1s
    ##  94450K .......... .......... .......... .......... .......... 70% 69.6M 1s
    ##  94500K .......... .......... .......... .......... .......... 70% 74.8M 1s
    ##  94550K .......... .......... .......... .......... .......... 70% 51.5M 1s
    ##  94600K .......... .......... .......... .......... .......... 70% 18.6M 1s
    ##  94650K .......... .......... .......... .......... .......... 70% 57.0M 1s
    ##  94700K .......... .......... .......... .......... .......... 70% 74.4M 1s
    ##  94750K .......... .......... .......... .......... .......... 70% 71.2M 1s
    ##  94800K .......... .......... .......... .......... .......... 70% 22.0M 1s
    ##  94850K .......... .......... .......... .......... .......... 70% 56.7M 1s
    ##  94900K .......... .......... .......... .......... .......... 70% 51.9M 1s
    ##  94950K .......... .......... .......... .......... .......... 70% 62.0M 1s
    ##  95000K .......... .......... .......... .......... .......... 70% 66.8M 1s
    ##  95050K .......... .......... .......... .......... .......... 70% 40.4M 1s
    ##  95100K .......... .......... .......... .......... .......... 70% 47.5M 1s
    ##  95150K .......... .......... .......... .......... .......... 70% 82.6M 1s
    ##  95200K .......... .......... .......... .......... .......... 70% 53.8M 1s
    ##  95250K .......... .......... .......... .......... .......... 70% 19.4M 1s
    ##  95300K .......... .......... .......... .......... .......... 70% 71.4M 1s
    ##  95350K .......... .......... .......... .......... .......... 70% 64.3M 1s
    ##  95400K .......... .......... .......... .......... .......... 70% 34.4M 1s
    ##  95450K .......... .......... .......... .......... .......... 70% 43.4M 1s
    ##  95500K .......... .......... .......... .......... .......... 70% 39.9M 1s
    ##  95550K .......... .......... .......... .......... .......... 70% 62.9M 1s
    ##  95600K .......... .......... .......... .......... .......... 70% 61.6M 1s
    ##  95650K .......... .......... .......... .......... .......... 71% 72.3M 1s
    ##  95700K .......... .......... .......... .......... .......... 71% 57.5M 1s
    ##  95750K .......... .......... .......... .......... .......... 71% 49.4M 1s
    ##  95800K .......... .......... .......... .......... .......... 71% 47.7M 1s
    ##  95850K .......... .......... .......... .......... .......... 71% 59.3M 1s
    ##  95900K .......... .......... .......... .......... .......... 71% 35.8M 1s
    ##  95950K .......... .......... .......... .......... .......... 71% 74.7M 1s
    ##  96000K .......... .......... .......... .......... .......... 71% 56.1M 1s
    ##  96050K .......... .......... .......... .......... .......... 71% 50.9M 1s
    ##  96100K .......... .......... .......... .......... .......... 71% 31.2M 1s
    ##  96150K .......... .......... .......... .......... .......... 71% 57.8M 1s
    ##  96200K .......... .......... .......... .......... .......... 71% 58.1M 1s
    ##  96250K .......... .......... .......... .......... .......... 71% 75.4M 1s
    ##  96300K .......... .......... .......... .......... .......... 71% 48.1M 1s
    ##  96350K .......... .......... .......... .......... .......... 71% 65.1M 1s
    ##  96400K .......... .......... .......... .......... .......... 71% 60.1M 1s
    ##  96450K .......... .......... .......... .......... .......... 71% 29.8M 1s
    ##  96500K .......... .......... .......... .......... .......... 71% 52.0M 1s
    ##  96550K .......... .......... .......... .......... .......... 71% 65.1M 1s
    ##  96600K .......... .......... .......... .......... .......... 71% 63.7M 1s
    ##  96650K .......... .......... .......... .......... .......... 71% 26.2M 1s
    ##  96700K .......... .......... .......... .......... .......... 71% 51.0M 1s
    ##  96750K .......... .......... .......... .......... .......... 71% 69.6M 1s
    ##  96800K .......... .......... .......... .......... .......... 71% 51.4M 1s
    ##  96850K .......... .......... .......... .......... .......... 71% 78.4M 1s
    ##  96900K .......... .......... .......... .......... .......... 71% 19.9M 1s
    ##  96950K .......... .......... .......... .......... .......... 71% 69.3M 1s
    ##  97000K .......... .......... .......... .......... .......... 72% 70.1M 1s
    ##  97050K .......... .......... .......... .......... .......... 72% 79.3M 1s
    ##  97100K .......... .......... .......... .......... .......... 72% 60.8M 1s
    ##  97150K .......... .......... .......... .......... .......... 72% 24.1M 1s
    ##  97200K .......... .......... .......... .......... .......... 72% 60.1M 1s
    ##  97250K .......... .......... .......... .......... .......... 72% 77.5M 1s
    ##  97300K .......... .......... .......... .......... .......... 72% 78.0M 1s
    ##  97350K .......... .......... .......... .......... .......... 72% 22.4M 1s
    ##  97400K .......... .......... .......... .......... .......... 72% 48.0M 1s
    ##  97450K .......... .......... .......... .......... .......... 72% 52.9M 1s
    ##  97500K .......... .......... .......... .......... .......... 72% 61.5M 1s
    ##  97550K .......... .......... .......... .......... .......... 72% 68.3M 1s
    ##  97600K .......... .......... .......... .......... .......... 72% 65.3M 1s
    ##  97650K .......... .......... .......... .......... .......... 72% 81.9M 1s
    ##  97700K .......... .......... .......... .......... .......... 72% 45.9M 1s
    ##  97750K .......... .......... .......... .......... .......... 72% 63.4M 1s
    ##  97800K .......... .......... .......... .......... .......... 72% 63.5M 1s
    ##  97850K .......... .......... .......... .......... .......... 72% 67.0M 1s
    ##  97900K .......... .......... .......... .......... .......... 72% 33.5M 1s
    ##  97950K .......... .......... .......... .......... .......... 72% 45.1M 1s
    ##  98000K .......... .......... .......... .......... .......... 72% 31.7M 1s
    ##  98050K .......... .......... .......... .......... .......... 72% 34.5M 1s
    ##  98100K .......... .......... .......... .......... .......... 72% 38.7M 1s
    ##  98150K .......... .......... .......... .......... .......... 72% 49.6M 1s
    ##  98200K .......... .......... .......... .......... .......... 72% 35.4M 1s
    ##  98250K .......... .......... .......... .......... .......... 72% 50.7M 1s
    ##  98300K .......... .......... .......... .......... .......... 72% 14.1M 1s
    ##  98350K .......... .......... .......... .......... .......... 73% 37.2M 1s
    ##  98400K .......... .......... .......... .......... .......... 73% 41.0M 1s
    ##  98450K .......... .......... .......... .......... .......... 73% 49.9M 1s
    ##  98500K .......... .......... .......... .......... .......... 73% 41.0M 1s
    ##  98550K .......... .......... .......... .......... .......... 73% 45.2M 1s
    ##  98600K .......... .......... .......... .......... .......... 73% 38.4M 1s
    ##  98650K .......... .......... .......... .......... .......... 73% 41.8M 1s
    ##  98700K .......... .......... .......... .......... .......... 73% 41.9M 1s
    ##  98750K .......... .......... .......... .......... .......... 73% 52.9M 1s
    ##  98800K .......... .......... .......... .......... .......... 73% 44.9M 1s
    ##  98850K .......... .......... .......... .......... .......... 73% 10.6M 1s
    ##  98900K .......... .......... .......... .......... .......... 73% 34.0M 1s
    ##  98950K .......... .......... .......... .......... .......... 73% 43.1M 1s
    ##  99000K .......... .......... .......... .......... .......... 73% 42.8M 1s
    ##  99050K .......... .......... .......... .......... .......... 73% 47.6M 1s
    ##  99100K .......... .......... .......... .......... .......... 73% 44.4M 1s
    ##  99150K .......... .......... .......... .......... .......... 73% 54.2M 1s
    ##  99200K .......... .......... .......... .......... .......... 73% 44.9M 1s
    ##  99250K .......... .......... .......... .......... .......... 73% 47.2M 1s
    ##  99300K .......... .......... .......... .......... .......... 73% 40.1M 1s
    ##  99350K .......... .......... .......... .......... .......... 73% 37.4M 1s
    ##  99400K .......... .......... .......... .......... .......... 73% 43.0M 1s
    ##  99450K .......... .......... .......... .......... .......... 73% 53.6M 1s
    ##  99500K .......... .......... .......... .......... .......... 73% 46.9M 1s
    ##  99550K .......... .......... .......... .......... .......... 73% 54.9M 1s
    ##  99600K .......... .......... .......... .......... .......... 73% 11.4M 1s
    ##  99650K .......... .......... .......... .......... .......... 73% 41.8M 1s
    ##  99700K .......... .......... .......... .......... .......... 74% 44.2M 1s
    ##  99750K .......... .......... .......... .......... .......... 74% 57.7M 1s
    ##  99800K .......... .......... .......... .......... .......... 74% 54.4M 1s
    ##  99850K .......... .......... .......... .......... .......... 74% 55.7M 1s
    ##  99900K .......... .......... .......... .......... .......... 74% 37.4M 1s
    ##  99950K .......... .......... .......... .......... .......... 74% 45.0M 1s
    ## 100000K .......... .......... .......... .......... .......... 74% 43.6M 1s
    ## 100050K .......... .......... .......... .......... .......... 74% 56.9M 1s
    ## 100100K .......... .......... .......... .......... .......... 74% 53.0M 1s
    ## 100150K .......... .......... .......... .......... .......... 74% 50.0M 1s
    ## 100200K .......... .......... .......... .......... .......... 74% 47.6M 1s
    ## 100250K .......... .......... .......... .......... .......... 74% 44.4M 1s
    ## 100300K .......... .......... .......... .......... .......... 74% 47.6M 1s
    ## 100350K .......... .......... .......... .......... .......... 74% 69.5M 1s
    ## 100400K .......... .......... .......... .......... .......... 74% 49.9M 1s
    ## 100450K .......... .......... .......... .......... .......... 74% 44.7M 1s
    ## 100500K .......... .......... .......... .......... .......... 74% 46.1M 1s
    ## 100550K .......... .......... .......... .......... .......... 74% 46.2M 1s
    ## 100600K .......... .......... .......... .......... .......... 74% 43.6M 1s
    ## 100650K .......... .......... .......... .......... .......... 74% 65.1M 1s
    ## 100700K .......... .......... .......... .......... .......... 74% 61.2M 1s
    ## 100750K .......... .......... .......... .......... .......... 74% 54.7M 1s
    ## 100800K .......... .......... .......... .......... .......... 74% 40.4M 1s
    ## 100850K .......... .......... .......... .......... .......... 74% 49.0M 1s
    ## 100900K .......... .......... .......... .......... .......... 74% 50.7M 1s
    ## 100950K .......... .......... .......... .......... .......... 74% 78.8M 1s
    ## 101000K .......... .......... .......... .......... .......... 74% 42.4M 1s
    ## 101050K .......... .......... .......... .......... .......... 75% 43.1M 1s
    ## 101100K .......... .......... .......... .......... .......... 75% 53.8M 1s
    ## 101150K .......... .......... .......... .......... .......... 75% 68.2M 1s
    ## 101200K .......... .......... .......... .......... .......... 75% 64.4M 1s
    ## 101250K .......... .......... .......... .......... .......... 75% 72.4M 1s
    ## 101300K .......... .......... .......... .......... .......... 75% 54.4M 1s
    ## 101350K .......... .......... .......... .......... .......... 75% 45.6M 1s
    ## 101400K .......... .......... .......... .......... .......... 75% 63.2M 1s
    ## 101450K .......... .......... .......... .......... .......... 75% 70.8M 1s
    ## 101500K .......... .......... .......... .......... .......... 75% 53.7M 1s
    ## 101550K .......... .......... .......... .......... .......... 75% 67.5M 1s
    ## 101600K .......... .......... .......... .......... .......... 75% 69.8M 1s
    ## 101650K .......... .......... .......... .......... .......... 75% 67.0M 1s
    ## 101700K .......... .......... .......... .......... .......... 75% 51.8M 1s
    ## 101750K .......... .......... .......... .......... .......... 75% 54.7M 1s
    ## 101800K .......... .......... .......... .......... .......... 75% 39.7M 1s
    ## 101850K .......... .......... .......... .......... .......... 75% 60.7M 1s
    ## 101900K .......... .......... .......... .......... .......... 75% 55.8M 1s
    ## 101950K .......... .......... .......... .......... .......... 75% 62.6M 1s
    ## 102000K .......... .......... .......... .......... .......... 75% 62.1M 1s
    ## 102050K .......... .......... .......... .......... .......... 75% 65.1M 1s
    ## 102100K .......... .......... .......... .......... .......... 75% 58.5M 1s
    ## 102150K .......... .......... .......... .......... .......... 75% 75.4M 1s
    ## 102200K .......... .......... .......... .......... .......... 75% 63.5M 1s
    ## 102250K .......... .......... .......... .......... .......... 75% 73.2M 1s
    ## 102300K .......... .......... .......... .......... .......... 75% 57.3M 1s
    ## 102350K .......... .......... .......... .......... .......... 75% 54.0M 1s
    ## 102400K .......... .......... .......... .......... .......... 76% 50.5M 1s
    ## 102450K .......... .......... .......... .......... .......... 76% 66.3M 1s
    ## 102500K .......... .......... .......... .......... .......... 76% 46.4M 1s
    ## 102550K .......... .......... .......... .......... .......... 76% 63.5M 1s
    ## 102600K .......... .......... .......... .......... .......... 76% 70.3M 1s
    ## 102650K .......... .......... .......... .......... .......... 76% 84.2M 1s
    ## 102700K .......... .......... .......... .......... .......... 76% 59.2M 1s
    ## 102750K .......... .......... .......... .......... .......... 76% 65.4M 1s
    ## 102800K .......... .......... .......... .......... .......... 76% 61.0M 1s
    ## 102850K .......... .......... .......... .......... .......... 76% 61.6M 1s
    ## 102900K .......... .......... .......... .......... .......... 76% 53.9M 1s
    ## 102950K .......... .......... .......... .......... .......... 76% 75.1M 1s
    ## 103000K .......... .......... .......... .......... .......... 76% 42.4M 1s
    ## 103050K .......... .......... .......... .......... .......... 76% 72.0M 1s
    ## 103100K .......... .......... .......... .......... .......... 76% 63.5M 1s
    ## 103150K .......... .......... .......... .......... .......... 76% 60.3M 1s
    ## 103200K .......... .......... .......... .......... .......... 76% 59.9M 1s
    ## 103250K .......... .......... .......... .......... .......... 76% 58.6M 1s
    ## 103300K .......... .......... .......... .......... .......... 76% 63.9M 1s
    ## 103350K .......... .......... .......... .......... .......... 76% 63.1M 1s
    ## 103400K .......... .......... .......... .......... .......... 76% 57.0M 1s
    ## 103450K .......... .......... .......... .......... .......... 76% 70.4M 1s
    ## 103500K .......... .......... .......... .......... .......... 76% 57.6M 1s
    ## 103550K .......... .......... .......... .......... .......... 76% 76.2M 1s
    ## 103600K .......... .......... .......... .......... .......... 76% 64.6M 1s
    ## 103650K .......... .......... .......... .......... .......... 76% 71.9M 1s
    ## 103700K .......... .......... .......... .......... .......... 77% 57.4M 1s
    ## 103750K .......... .......... .......... .......... .......... 77% 55.3M 1s
    ## 103800K .......... .......... .......... .......... .......... 77% 68.3M 1s
    ## 103850K .......... .......... .......... .......... .......... 77% 41.9M 1s
    ## 103900K .......... .......... .......... .......... .......... 77% 57.9M 1s
    ## 103950K .......... .......... .......... .......... .......... 77% 88.6M 1s
    ## 104000K .......... .......... .......... .......... .......... 77% 33.0M 1s
    ## 104050K .......... .......... .......... .......... .......... 77% 82.5M 1s
    ## 104100K .......... .......... .......... .......... .......... 77% 71.9M 1s
    ## 104150K .......... .......... .......... .......... .......... 77% 89.5M 1s
    ## 104200K .......... .......... .......... .......... .......... 77% 22.9M 1s
    ## 104250K .......... .......... .......... .......... .......... 77% 77.0M 1s
    ## 104300K .......... .......... .......... .......... .......... 77% 59.2M 1s
    ## 104350K .......... .......... .......... .......... .......... 77% 62.2M 1s
    ## 104400K .......... .......... .......... .......... .......... 77% 65.6M 1s
    ## 104450K .......... .......... .......... .......... .......... 77% 88.2M 1s
    ## 104500K .......... .......... .......... .......... .......... 77% 36.5M 1s
    ## 104550K .......... .......... .......... .......... .......... 77% 45.9M 1s
    ## 104600K .......... .......... .......... .......... .......... 77% 65.3M 1s
    ## 104650K .......... .......... .......... .......... .......... 77% 75.6M 1s
    ## 104700K .......... .......... .......... .......... .......... 77% 70.0M 1s
    ## 104750K .......... .......... .......... .......... .......... 77% 80.0M 1s
    ## 104800K .......... .......... .......... .......... .......... 77% 30.4M 1s
    ## 104850K .......... .......... .......... .......... .......... 77% 59.5M 1s
    ## 104900K .......... .......... .......... .......... .......... 77% 59.6M 1s
    ## 104950K .......... .......... .......... .......... .......... 77% 48.2M 1s
    ## 105000K .......... .......... .......... .......... .......... 77% 44.0M 1s
    ## 105050K .......... .......... .......... .......... .......... 78% 68.5M 1s
    ## 105100K .......... .......... .......... .......... .......... 78% 39.7M 1s
    ## 105150K .......... .......... .......... .......... .......... 78% 65.8M 1s
    ## 105200K .......... .......... .......... .......... .......... 78% 68.5M 1s
    ## 105250K .......... .......... .......... .......... .......... 78% 71.1M 1s
    ## 105300K .......... .......... .......... .......... .......... 78% 61.4M 1s
    ## 105350K .......... .......... .......... .......... .......... 78% 56.3M 1s
    ## 105400K .......... .......... .......... .......... .......... 78% 34.9M 1s
    ## 105450K .......... .......... .......... .......... .......... 78% 90.2M 1s
    ## 105500K .......... .......... .......... .......... .......... 78% 82.1M 1s
    ## 105550K .......... .......... .......... .......... .......... 78% 63.1M 1s
    ## 105600K .......... .......... .......... .......... .......... 78% 70.0M 1s
    ## 105650K .......... .......... .......... .......... .......... 78% 35.1M 1s
    ## 105700K .......... .......... .......... .......... .......... 78% 76.8M 1s
    ## 105750K .......... .......... .......... .......... .......... 78% 56.0M 1s
    ## 105800K .......... .......... .......... .......... .......... 78% 58.5M 1s
    ## 105850K .......... .......... .......... .......... .......... 78% 61.7M 1s
    ## 105900K .......... .......... .......... .......... .......... 78% 68.8M 1s
    ## 105950K .......... .......... .......... .......... .......... 78% 22.9M 1s
    ## 106000K .......... .......... .......... .......... .......... 78% 62.2M 1s
    ## 106050K .......... .......... .......... .......... .......... 78% 73.0M 1s
    ## 106100K .......... .......... .......... .......... .......... 78% 73.3M 1s
    ## 106150K .......... .......... .......... .......... .......... 78% 65.2M 1s
    ## 106200K .......... .......... .......... .......... .......... 78% 64.5M 1s
    ## 106250K .......... .......... .......... .......... .......... 78% 34.3M 1s
    ## 106300K .......... .......... .......... .......... .......... 78% 64.1M 1s
    ## 106350K .......... .......... .......... .......... .......... 78% 70.4M 1s
    ## 106400K .......... .......... .......... .......... .......... 79% 58.0M 1s
    ## 106450K .......... .......... .......... .......... .......... 79% 71.6M 1s
    ## 106500K .......... .......... .......... .......... .......... 79% 88.7M 1s
    ## 106550K .......... .......... .......... .......... .......... 79% 38.3M 1s
    ## 106600K .......... .......... .......... .......... .......... 79% 61.2M 1s
    ## 106650K .......... .......... .......... .......... .......... 79% 57.3M 1s
    ## 106700K .......... .......... .......... .......... .......... 79% 68.5M 1s
    ## 106750K .......... .......... .......... .......... .......... 79% 64.8M 1s
    ## 106800K .......... .......... .......... .......... .......... 79% 81.1M 1s
    ## 106850K .......... .......... .......... .......... .......... 79% 80.5M 1s
    ## 106900K .......... .......... .......... .......... .......... 79% 59.2M 1s
    ## 106950K .......... .......... .......... .......... .......... 79% 83.9M 1s
    ## 107000K .......... .......... .......... .......... .......... 79% 72.9M 1s
    ## 107050K .......... .......... .......... .......... .......... 79% 67.0M 1s
    ## 107100K .......... .......... .......... .......... .......... 79% 54.0M 1s
    ## 107150K .......... .......... .......... .......... .......... 79% 75.6M 1s
    ## 107200K .......... .......... .......... .......... .......... 79% 48.3M 1s
    ## 107250K .......... .......... .......... .......... .......... 79% 74.6M 1s
    ## 107300K .......... .......... .......... .......... .......... 79% 66.8M 1s
    ## 107350K .......... .......... .......... .......... .......... 79% 65.1M 1s
    ## 107400K .......... .......... .......... .......... .......... 79% 39.1M 1s
    ## 107450K .......... .......... .......... .......... .......... 79% 59.8M 1s
    ## 107500K .......... .......... .......... .......... .......... 79% 71.0M 1s
    ## 107550K .......... .......... .......... .......... .......... 79% 14.5M 1s
    ## 107600K .......... .......... .......... .......... .......... 79% 60.7M 1s
    ## 107650K .......... .......... .......... .......... .......... 79% 72.7M 1s
    ## 107700K .......... .......... .......... .......... .......... 79% 70.4M 1s
    ## 107750K .......... .......... .......... .......... .......... 80% 84.3M 1s
    ## 107800K .......... .......... .......... .......... .......... 80% 32.5M 1s
    ## 107850K .......... .......... .......... .......... .......... 80% 81.3M 1s
    ## 107900K .......... .......... .......... .......... .......... 80% 57.1M 1s
    ## 107950K .......... .......... .......... .......... .......... 80% 58.0M 1s
    ## 108000K .......... .......... .......... .......... .......... 80% 59.4M 1s
    ## 108050K .......... .......... .......... .......... .......... 80% 70.3M 1s
    ## 108100K .......... .......... .......... .......... .......... 80% 26.1M 1s
    ## 108150K .......... .......... .......... .......... .......... 80% 32.9M 1s
    ## 108200K .......... .......... .......... .......... .......... 80% 58.1M 1s
    ## 108250K .......... .......... .......... .......... .......... 80% 70.5M 1s
    ## 108300K .......... .......... .......... .......... .......... 80% 77.2M 1s
    ## 108350K .......... .......... .......... .......... .......... 80% 80.1M 1s
    ## 108400K .......... .......... .......... .......... .......... 80% 61.9M 1s
    ## 108450K .......... .......... .......... .......... .......... 80% 76.7M 1s
    ## 108500K .......... .......... .......... .......... .......... 80% 56.8M 1s
    ## 108550K .......... .......... .......... .......... .......... 80% 45.8M 1s
    ## 108600K .......... .......... .......... .......... .......... 80% 59.3M 1s
    ## 108650K .......... .......... .......... .......... .......... 80% 75.2M 1s
    ## 108700K .......... .......... .......... .......... .......... 80% 26.9M 1s
    ## 108750K .......... .......... .......... .......... .......... 80% 39.5M 1s
    ## 108800K .......... .......... .......... .......... .......... 80% 68.3M 1s
    ## 108850K .......... .......... .......... .......... .......... 80% 58.6M 1s
    ## 108900K .......... .......... .......... .......... .......... 80% 52.9M 1s
    ## 108950K .......... .......... .......... .......... .......... 80% 70.6M 1s
    ## 109000K .......... .......... .......... .......... .......... 80% 74.7M 1s
    ## 109050K .......... .......... .......... .......... .......... 80% 81.1M 1s
    ## 109100K .......... .......... .......... .......... .......... 81% 73.2M 1s
    ## 109150K .......... .......... .......... .......... .......... 81% 64.9M 1s
    ## 109200K .......... .......... .......... .......... .......... 81% 65.0M 1s
    ## 109250K .......... .......... .......... .......... .......... 81% 82.0M 1s
    ## 109300K .......... .......... .......... .......... .......... 81% 68.0M 1s
    ## 109350K .......... .......... .......... .......... .......... 81% 89.2M 1s
    ## 109400K .......... .......... .......... .......... .......... 81% 26.8M 1s
    ## 109450K .......... .......... .......... .......... .......... 81% 67.0M 1s
    ## 109500K .......... .......... .......... .......... .......... 81% 62.2M 1s
    ## 109550K .......... .......... .......... .......... .......... 81% 72.1M 1s
    ## 109600K .......... .......... .......... .......... .......... 81% 69.3M 0s
    ## 109650K .......... .......... .......... .......... .......... 81% 44.7M 0s
    ## 109700K .......... .......... .......... .......... .......... 81% 53.7M 0s
    ## 109750K .......... .......... .......... .......... .......... 81% 66.0M 0s
    ## 109800K .......... .......... .......... .......... .......... 81% 60.3M 0s
    ## 109850K .......... .......... .......... .......... .......... 81% 84.9M 0s
    ## 109900K .......... .......... .......... .......... .......... 81% 62.5M 0s
    ## 109950K .......... .......... .......... .......... .......... 81% 80.0M 0s
    ## 110000K .......... .......... .......... .......... .......... 81% 68.6M 0s
    ## 110050K .......... .......... .......... .......... .......... 81% 60.5M 0s
    ## 110100K .......... .......... .......... .......... .......... 81% 61.2M 0s
    ## 110150K .......... .......... .......... .......... .......... 81% 59.4M 0s
    ## 110200K .......... .......... .......... .......... .......... 81% 73.9M 0s
    ## 110250K .......... .......... .......... .......... .......... 81% 83.1M 0s
    ## 110300K .......... .......... .......... .......... .......... 81% 58.6M 0s
    ## 110350K .......... .......... .......... .......... .......... 81% 58.3M 0s
    ## 110400K .......... .......... .......... .......... .......... 81% 54.3M 0s
    ## 110450K .......... .......... .......... .......... .......... 82% 67.4M 0s
    ## 110500K .......... .......... .......... .......... .......... 82% 60.8M 0s
    ## 110550K .......... .......... .......... .......... .......... 82% 38.5M 0s
    ## 110600K .......... .......... .......... .......... .......... 82% 53.9M 0s
    ## 110650K .......... .......... .......... .......... .......... 82% 17.6M 0s
    ## 110700K .......... .......... .......... .......... .......... 82% 19.3M 0s
    ## 110750K .......... .......... .......... .......... .......... 82% 84.1M 0s
    ## 110800K .......... .......... .......... .......... .......... 82% 67.4M 0s
    ## 110850K .......... .......... .......... .......... .......... 82% 88.9M 0s
    ## 110900K .......... .......... .......... .......... .......... 82% 21.6M 0s
    ## 110950K .......... .......... .......... .......... .......... 82% 56.8M 0s
    ## 111000K .......... .......... .......... .......... .......... 82% 60.5M 0s
    ## 111050K .......... .......... .......... .......... .......... 82% 81.7M 0s
    ## 111100K .......... .......... .......... .......... .......... 82% 65.6M 0s
    ## 111150K .......... .......... .......... .......... .......... 82% 77.4M 0s
    ## 111200K .......... .......... .......... .......... .......... 82% 67.3M 0s
    ## 111250K .......... .......... .......... .......... .......... 82% 69.2M 0s
    ## 111300K .......... .......... .......... .......... .......... 82% 61.9M 0s
    ## 111350K .......... .......... .......... .......... .......... 82% 77.9M 0s
    ## 111400K .......... .......... .......... .......... .......... 82% 54.9M 0s
    ## 111450K .......... .......... .......... .......... .......... 82% 76.0M 0s
    ## 111500K .......... .......... .......... .......... .......... 82% 85.0M 0s
    ## 111550K .......... .......... .......... .......... .......... 82% 4.50M 0s
    ## 111600K .......... .......... .......... .......... .......... 82% 64.2M 0s
    ## 111650K .......... .......... .......... .......... .......... 82% 10.7M 0s
    ## 111700K .......... .......... .......... .......... .......... 82% 64.6M 0s
    ## 111750K .......... .......... .......... .......... .......... 82% 83.6M 0s
    ## 111800K .......... .......... .......... .......... .......... 83% 72.0M 0s
    ## 111850K .......... .......... .......... .......... .......... 83% 68.4M 0s
    ## 111900K .......... .......... .......... .......... .......... 83% 53.4M 0s
    ## 111950K .......... .......... .......... .......... .......... 83% 84.7M 0s
    ## 112000K .......... .......... .......... .......... .......... 83% 31.0M 0s
    ## 112050K .......... .......... .......... .......... .......... 83% 98.4M 0s
    ## 112100K .......... .......... .......... .......... .......... 83% 70.9M 0s
    ## 112150K .......... .......... .......... .......... .......... 83% 77.7M 0s
    ## 112200K .......... .......... .......... .......... .......... 83% 36.9M 0s
    ## 112250K .......... .......... .......... .......... .......... 83% 21.3M 0s
    ## 112300K .......... .......... .......... .......... .......... 83% 75.4M 0s
    ## 112350K .......... .......... .......... .......... .......... 83% 92.6M 0s
    ## 112400K .......... .......... .......... .......... .......... 83% 75.6M 0s
    ## 112450K .......... .......... .......... .......... .......... 83% 85.1M 0s
    ## 112500K .......... .......... .......... .......... .......... 83% 22.9M 0s
    ## 112550K .......... .......... .......... .......... .......... 83% 82.2M 0s
    ## 112600K .......... .......... .......... .......... .......... 83% 26.2M 0s
    ## 112650K .......... .......... .......... .......... .......... 83% 70.5M 0s
    ## 112700K .......... .......... .......... .......... .......... 83% 46.7M 0s
    ## 112750K .......... .......... .......... .......... .......... 83% 77.1M 0s
    ## 112800K .......... .......... .......... .......... .......... 83% 83.1M 0s
    ## 112850K .......... .......... .......... .......... .......... 83% 98.9M 0s
    ## 112900K .......... .......... .......... .......... .......... 83% 21.5M 0s
    ## 112950K .......... .......... .......... .......... .......... 83% 99.4M 0s
    ## 113000K .......... .......... .......... .......... .......... 83% 64.3M 0s
    ## 113050K .......... .......... .......... .......... .......... 83% 88.4M 0s
    ## 113100K .......... .......... .......... .......... .......... 83% 82.0M 0s
    ## 113150K .......... .......... .......... .......... .......... 84% 99.5M 0s
    ## 113200K .......... .......... .......... .......... .......... 84% 43.3M 0s
    ## 113250K .......... .......... .......... .......... .......... 84% 79.0M 0s
    ## 113300K .......... .......... .......... .......... .......... 84% 36.4M 0s
    ## 113350K .......... .......... .......... .......... .......... 84% 82.0M 0s
    ## 113400K .......... .......... .......... .......... .......... 84% 48.1M 0s
    ## 113450K .......... .......... .......... .......... .......... 84% 61.3M 0s
    ## 113500K .......... .......... .......... .......... .......... 84% 57.1M 0s
    ## 113550K .......... .......... .......... .......... .......... 84% 61.6M 0s
    ## 113600K .......... .......... .......... .......... .......... 84% 8.75M 0s
    ## 113650K .......... .......... .......... .......... .......... 84% 71.2M 0s
    ## 113700K .......... .......... .......... .......... .......... 84% 70.1M 0s
    ## 113750K .......... .......... .......... .......... .......... 84% 78.0M 0s
    ## 113800K .......... .......... .......... .......... .......... 84% 70.5M 0s
    ## 113850K .......... .......... .......... .......... .......... 84% 78.4M 0s
    ## 113900K .......... .......... .......... .......... .......... 84% 66.3M 0s
    ## 113950K .......... .......... .......... .......... .......... 84% 66.5M 0s
    ## 114000K .......... .......... .......... .......... .......... 84% 51.6M 0s
    ## 114050K .......... .......... .......... .......... .......... 84% 73.7M 0s
    ## 114100K .......... .......... .......... .......... .......... 84% 65.5M 0s
    ## 114150K .......... .......... .......... .......... .......... 84% 69.2M 0s
    ## 114200K .......... .......... .......... .......... .......... 84% 61.5M 0s
    ## 114250K .......... .......... .......... .......... .......... 84% 54.7M 0s
    ## 114300K .......... .......... .......... .......... .......... 84% 63.3M 0s
    ## 114350K .......... .......... .......... .......... .......... 84% 69.3M 0s
    ## 114400K .......... .......... .......... .......... .......... 84% 57.8M 0s
    ## 114450K .......... .......... .......... .......... .......... 84% 74.6M 0s
    ## 114500K .......... .......... .......... .......... .......... 85% 68.3M 0s
    ## 114550K .......... .......... .......... .......... .......... 85% 72.9M 0s
    ## 114600K .......... .......... .......... .......... .......... 85% 74.6M 0s
    ## 114650K .......... .......... .......... .......... .......... 85% 64.7M 0s
    ## 114700K .......... .......... .......... .......... .......... 85% 65.6M 0s
    ## 114750K .......... .......... .......... .......... .......... 85% 87.7M 0s
    ## 114800K .......... .......... .......... .......... .......... 85% 69.0M 0s
    ## 114850K .......... .......... .......... .......... .......... 85% 66.5M 0s
    ## 114900K .......... .......... .......... .......... .......... 85% 41.2M 0s
    ## 114950K .......... .......... .......... .......... .......... 85% 61.5M 0s
    ## 115000K .......... .......... .......... .......... .......... 85% 58.0M 0s
    ## 115050K .......... .......... .......... .......... .......... 85% 66.2M 0s
    ## 115100K .......... .......... .......... .......... .......... 85% 64.1M 0s
    ## 115150K .......... .......... .......... .......... .......... 85% 61.2M 0s
    ## 115200K .......... .......... .......... .......... .......... 85% 78.8M 0s
    ## 115250K .......... .......... .......... .......... .......... 85% 52.8M 0s
    ## 115300K .......... .......... .......... .......... .......... 85% 33.4M 0s
    ## 115350K .......... .......... .......... .......... .......... 85% 86.9M 0s
    ## 115400K .......... .......... .......... .......... .......... 85% 69.2M 0s
    ## 115450K .......... .......... .......... .......... .......... 85% 72.6M 0s
    ## 115500K .......... .......... .......... .......... .......... 85% 19.5M 0s
    ## 115550K .......... .......... .......... .......... .......... 85% 37.4M 0s
    ## 115600K .......... .......... .......... .......... .......... 85% 58.3M 0s
    ## 115650K .......... .......... .......... .......... .......... 85% 71.1M 0s
    ## 115700K .......... .......... .......... .......... .......... 85% 70.5M 0s
    ## 115750K .......... .......... .......... .......... .......... 85% 56.1M 0s
    ## 115800K .......... .......... .......... .......... .......... 85% 39.8M 0s
    ## 115850K .......... .......... .......... .......... .......... 86% 58.9M 0s
    ## 115900K .......... .......... .......... .......... .......... 86% 62.3M 0s
    ## 115950K .......... .......... .......... .......... .......... 86% 72.5M 0s
    ## 116000K .......... .......... .......... .......... .......... 86% 67.7M 0s
    ## 116050K .......... .......... .......... .......... .......... 86% 76.5M 0s
    ## 116100K .......... .......... .......... .......... .......... 86% 66.8M 0s
    ## 116150K .......... .......... .......... .......... .......... 86% 77.1M 0s
    ## 116200K .......... .......... .......... .......... .......... 86% 65.3M 0s
    ## 116250K .......... .......... .......... .......... .......... 86% 78.2M 0s
    ## 116300K .......... .......... .......... .......... .......... 86% 44.9M 0s
    ## 116350K .......... .......... .......... .......... .......... 86% 64.9M 0s
    ## 116400K .......... .......... .......... .......... .......... 86% 63.9M 0s
    ## 116450K .......... .......... .......... .......... .......... 86% 77.8M 0s
    ## 116500K .......... .......... .......... .......... .......... 86% 58.6M 0s
    ## 116550K .......... .......... .......... .......... .......... 86% 65.1M 0s
    ## 116600K .......... .......... .......... .......... .......... 86% 67.5M 0s
    ## 116650K .......... .......... .......... .......... .......... 86% 70.1M 0s
    ## 116700K .......... .......... .......... .......... .......... 86% 74.0M 0s
    ## 116750K .......... .......... .......... .......... .......... 86% 98.7M 0s
    ## 116800K .......... .......... .......... .......... .......... 86% 27.8M 0s
    ## 116850K .......... .......... .......... .......... .......... 86% 95.8M 0s
    ## 116900K .......... .......... .......... .......... .......... 86% 78.9M 0s
    ## 116950K .......... .......... .......... .......... .......... 86% 35.3M 0s
    ## 117000K .......... .......... .......... .......... .......... 86% 63.0M 0s
    ## 117050K .......... .......... .......... .......... .......... 86% 82.9M 0s
    ## 117100K .......... .......... .......... .......... .......... 86% 57.3M 0s
    ## 117150K .......... .......... .......... .......... .......... 86% 84.9M 0s
    ## 117200K .......... .......... .......... .......... .......... 87% 40.7M 0s
    ## 117250K .......... .......... .......... .......... .......... 87% 61.4M 0s
    ## 117300K .......... .......... .......... .......... .......... 87% 61.5M 0s
    ## 117350K .......... .......... .......... .......... .......... 87% 71.1M 0s
    ## 117400K .......... .......... .......... .......... .......... 87% 68.5M 0s
    ## 117450K .......... .......... .......... .......... .......... 87% 68.6M 0s
    ## 117500K .......... .......... .......... .......... .......... 87% 56.7M 0s
    ## 117550K .......... .......... .......... .......... .......... 87% 53.3M 0s
    ## 117600K .......... .......... .......... .......... .......... 87% 68.7M 0s
    ## 117650K .......... .......... .......... .......... .......... 87% 35.0M 0s
    ## 117700K .......... .......... .......... .......... .......... 87% 63.2M 0s
    ## 117750K .......... .......... .......... .......... .......... 87% 68.6M 0s
    ## 117800K .......... .......... .......... .......... .......... 87% 60.7M 0s
    ## 117850K .......... .......... .......... .......... .......... 87% 66.3M 0s
    ## 117900K .......... .......... .......... .......... .......... 87% 73.8M 0s
    ## 117950K .......... .......... .......... .......... .......... 87% 96.7M 0s
    ## 118000K .......... .......... .......... .......... .......... 87% 56.4M 0s
    ## 118050K .......... .......... .......... .......... .......... 87% 80.2M 0s
    ## 118100K .......... .......... .......... .......... .......... 87% 61.0M 0s
    ## 118150K .......... .......... .......... .......... .......... 87% 60.1M 0s
    ## 118200K .......... .......... .......... .......... .......... 87% 76.7M 0s
    ## 118250K .......... .......... .......... .......... .......... 87% 80.9M 0s
    ## 118300K .......... .......... .......... .......... .......... 87% 62.5M 0s
    ## 118350K .......... .......... .......... .......... .......... 87% 72.1M 0s
    ## 118400K .......... .......... .......... .......... .......... 87% 39.9M 0s
    ## 118450K .......... .......... .......... .......... .......... 87% 56.1M 0s
    ## 118500K .......... .......... .......... .......... .......... 87% 58.6M 0s
    ## 118550K .......... .......... .......... .......... .......... 88% 82.3M 0s
    ## 118600K .......... .......... .......... .......... .......... 88% 67.5M 0s
    ## 118650K .......... .......... .......... .......... .......... 88% 62.0M 0s
    ## 118700K .......... .......... .......... .......... .......... 88% 27.5M 0s
    ## 118750K .......... .......... .......... .......... .......... 88% 61.5M 0s
    ## 118800K .......... .......... .......... .......... .......... 88% 80.4M 0s
    ## 118850K .......... .......... .......... .......... .......... 88% 72.1M 0s
    ## 118900K .......... .......... .......... .......... .......... 88% 78.3M 0s
    ## 118950K .......... .......... .......... .......... .......... 88% 62.5M 0s
    ## 119000K .......... .......... .......... .......... .......... 88% 69.1M 0s
    ## 119050K .......... .......... .......... .......... .......... 88% 77.3M 0s
    ## 119100K .......... .......... .......... .......... .......... 88% 30.8M 0s
    ## 119150K .......... .......... .......... .......... .......... 88% 76.2M 0s
    ## 119200K .......... .......... .......... .......... .......... 88% 78.4M 0s
    ## 119250K .......... .......... .......... .......... .......... 88% 60.4M 0s
    ## 119300K .......... .......... .......... .......... .......... 88% 66.8M 0s
    ## 119350K .......... .......... .......... .......... .......... 88% 30.7M 0s
    ## 119400K .......... .......... .......... .......... .......... 88% 36.6M 0s
    ## 119450K .......... .......... .......... .......... .......... 88% 80.5M 0s
    ## 119500K .......... .......... .......... .......... .......... 88% 87.3M 0s
    ## 119550K .......... .......... .......... .......... .......... 88% 79.4M 0s
    ## 119600K .......... .......... .......... .......... .......... 88% 97.7M 0s
    ## 119650K .......... .......... .......... .......... .......... 88% 84.9M 0s
    ## 119700K .......... .......... .......... .......... .......... 88% 63.1M 0s
    ## 119750K .......... .......... .......... .......... .......... 88% 18.2M 0s
    ## 119800K .......... .......... .......... .......... .......... 88% 83.8M 0s
    ## 119850K .......... .......... .......... .......... .......... 88%  102M 0s
    ## 119900K .......... .......... .......... .......... .......... 89% 80.9M 0s
    ## 119950K .......... .......... .......... .......... .......... 89% 87.1M 0s
    ## 120000K .......... .......... .......... .......... .......... 89% 30.2M 0s
    ## 120050K .......... .......... .......... .......... .......... 89% 77.7M 0s
    ## 120100K .......... .......... .......... .......... .......... 89% 63.3M 0s
    ## 120150K .......... .......... .......... .......... .......... 89% 64.0M 0s
    ## 120200K .......... .......... .......... .......... .......... 89% 74.5M 0s
    ## 120250K .......... .......... .......... .......... .......... 89% 13.7M 0s
    ## 120300K .......... .......... .......... .......... .......... 89% 61.5M 0s
    ## 120350K .......... .......... .......... .......... .......... 89% 51.6M 0s
    ## 120400K .......... .......... .......... .......... .......... 89% 63.5M 0s
    ## 120450K .......... .......... .......... .......... .......... 89% 80.2M 0s
    ## 120500K .......... .......... .......... .......... .......... 89% 63.1M 0s
    ## 120550K .......... .......... .......... .......... .......... 89% 74.7M 0s
    ## 120600K .......... .......... .......... .......... .......... 89% 68.6M 0s
    ## 120650K .......... .......... .......... .......... .......... 89% 82.4M 0s
    ## 120700K .......... .......... .......... .......... .......... 89% 67.9M 0s
    ## 120750K .......... .......... .......... .......... .......... 89% 73.5M 0s
    ## 120800K .......... .......... .......... .......... .......... 89% 70.7M 0s
    ## 120850K .......... .......... .......... .......... .......... 89% 73.5M 0s
    ## 120900K .......... .......... .......... .......... .......... 89% 79.3M 0s
    ## 120950K .......... .......... .......... .......... .......... 89% 63.0M 0s
    ## 121000K .......... .......... .......... .......... .......... 89% 62.3M 0s
    ## 121050K .......... .......... .......... .......... .......... 89% 26.6M 0s
    ## 121100K .......... .......... .......... .......... .......... 89% 54.9M 0s
    ## 121150K .......... .......... .......... .......... .......... 89% 71.7M 0s
    ## 121200K .......... .......... .......... .......... .......... 89% 32.3M 0s
    ## 121250K .......... .......... .......... .......... .......... 90% 31.4M 0s
    ## 121300K .......... .......... .......... .......... .......... 90% 57.0M 0s
    ## 121350K .......... .......... .......... .......... .......... 90% 81.9M 0s
    ## 121400K .......... .......... .......... .......... .......... 90% 71.2M 0s
    ## 121450K .......... .......... .......... .......... .......... 90% 66.0M 0s
    ## 121500K .......... .......... .......... .......... .......... 90% 78.4M 0s
    ## 121550K .......... .......... .......... .......... .......... 90% 66.3M 0s
    ## 121600K .......... .......... .......... .......... .......... 90% 44.2M 0s
    ## 121650K .......... .......... .......... .......... .......... 90% 81.0M 0s
    ## 121700K .......... .......... .......... .......... .......... 90% 43.3M 0s
    ## 121750K .......... .......... .......... .......... .......... 90% 49.8M 0s
    ## 121800K .......... .......... .......... .......... .......... 90% 82.5M 0s
    ## 121850K .......... .......... .......... .......... .......... 90% 81.2M 0s
    ## 121900K .......... .......... .......... .......... .......... 90% 80.4M 0s
    ## 121950K .......... .......... .......... .......... .......... 90% 73.2M 0s
    ## 122000K .......... .......... .......... .......... .......... 90% 72.5M 0s
    ## 122050K .......... .......... .......... .......... .......... 90% 39.2M 0s
    ## 122100K .......... .......... .......... .......... .......... 90% 32.9M 0s
    ## 122150K .......... .......... .......... .......... .......... 90% 73.2M 0s
    ## 122200K .......... .......... .......... .......... .......... 90% 86.0M 0s
    ## 122250K .......... .......... .......... .......... .......... 90%  103M 0s
    ## 122300K .......... .......... .......... .......... .......... 90% 77.2M 0s
    ## 122350K .......... .......... .......... .......... .......... 90% 90.9M 0s
    ## 122400K .......... .......... .......... .......... .......... 90% 87.7M 0s
    ## 122450K .......... .......... .......... .......... .......... 90% 75.4M 0s
    ## 122500K .......... .......... .......... .......... .......... 90% 75.7M 0s
    ## 122550K .......... .......... .......... .......... .......... 90% 98.3M 0s
    ## 122600K .......... .......... .......... .......... .......... 91% 81.3M 0s
    ## 122650K .......... .......... .......... .......... .......... 91% 73.1M 0s
    ## 122700K .......... .......... .......... .......... .......... 91% 37.3M 0s
    ## 122750K .......... .......... .......... .......... .......... 91% 68.1M 0s
    ## 122800K .......... .......... .......... .......... .......... 91% 74.2M 0s
    ## 122850K .......... .......... .......... .......... .......... 91% 94.5M 0s
    ## 122900K .......... .......... .......... .......... .......... 91% 73.0M 0s
    ## 122950K .......... .......... .......... .......... .......... 91% 65.2M 0s
    ## 123000K .......... .......... .......... .......... .......... 91% 64.2M 0s
    ## 123050K .......... .......... .......... .......... .......... 91% 37.6M 0s
    ## 123100K .......... .......... .......... .......... .......... 91% 76.5M 0s
    ## 123150K .......... .......... .......... .......... .......... 91% 96.8M 0s
    ## 123200K .......... .......... .......... .......... .......... 91% 65.1M 0s
    ## 123250K .......... .......... .......... .......... .......... 91% 62.4M 0s
    ## 123300K .......... .......... .......... .......... .......... 91% 63.8M 0s
    ## 123350K .......... .......... .......... .......... .......... 91% 57.5M 0s
    ## 123400K .......... .......... .......... .......... .......... 91% 62.9M 0s
    ## 123450K .......... .......... .......... .......... .......... 91% 63.7M 0s
    ## 123500K .......... .......... .......... .......... .......... 91% 68.9M 0s
    ## 123550K .......... .......... .......... .......... .......... 91% 66.5M 0s
    ## 123600K .......... .......... .......... .......... .......... 91% 83.3M 0s
    ## 123650K .......... .......... .......... .......... .......... 91% 72.7M 0s
    ## 123700K .......... .......... .......... .......... .......... 91% 61.9M 0s
    ## 123750K .......... .......... .......... .......... .......... 91% 60.8M 0s
    ## 123800K .......... .......... .......... .......... .......... 91% 67.8M 0s
    ## 123850K .......... .......... .......... .......... .......... 91% 75.0M 0s
    ## 123900K .......... .......... .......... .......... .......... 91% 63.1M 0s
    ## 123950K .......... .......... .......... .......... .......... 92% 72.8M 0s
    ## 124000K .......... .......... .......... .......... .......... 92% 90.0M 0s
    ## 124050K .......... .......... .......... .......... .......... 92% 94.5M 0s
    ## 124100K .......... .......... .......... .......... .......... 92% 18.4M 0s
    ## 124150K .......... .......... .......... .......... .......... 92% 75.0M 0s
    ## 124200K .......... .......... .......... .......... .......... 92% 65.9M 0s
    ## 124250K .......... .......... .......... .......... .......... 92% 81.6M 0s
    ## 124300K .......... .......... .......... .......... .......... 92% 80.3M 0s
    ## 124350K .......... .......... .......... .......... .......... 92% 96.4M 0s
    ## 124400K .......... .......... .......... .......... .......... 92% 84.1M 0s
    ## 124450K .......... .......... .......... .......... .......... 92% 79.6M 0s
    ## 124500K .......... .......... .......... .......... .......... 92% 34.0M 0s
    ## 124550K .......... .......... .......... .......... .......... 92% 73.4M 0s
    ## 124600K .......... .......... .......... .......... .......... 92% 33.9M 0s
    ## 124650K .......... .......... .......... .......... .......... 92% 54.4M 0s
    ## 124700K .......... .......... .......... .......... .......... 92% 75.1M 0s
    ## 124750K .......... .......... .......... .......... .......... 92% 84.1M 0s
    ## 124800K .......... .......... .......... .......... .......... 92% 65.9M 0s
    ## 124850K .......... .......... .......... .......... .......... 92% 79.3M 0s
    ## 124900K .......... .......... .......... .......... .......... 92% 21.4M 0s
    ## 124950K .......... .......... .......... .......... .......... 92% 39.6M 0s
    ## 125000K .......... .......... .......... .......... .......... 92% 88.3M 0s
    ## 125050K .......... .......... .......... .......... .......... 92%  101M 0s
    ## 125100K .......... .......... .......... .......... .......... 92% 79.9M 0s
    ## 125150K .......... .......... .......... .......... .......... 92% 73.9M 0s
    ## 125200K .......... .......... .......... .......... .......... 92% 24.7M 0s
    ## 125250K .......... .......... .......... .......... .......... 92% 41.8M 0s
    ## 125300K .......... .......... .......... .......... .......... 93% 67.0M 0s
    ## 125350K .......... .......... .......... .......... .......... 93% 75.5M 0s
    ## 125400K .......... .......... .......... .......... .......... 93% 80.3M 0s
    ## 125450K .......... .......... .......... .......... .......... 93% 70.3M 0s
    ## 125500K .......... .......... .......... .......... .......... 93% 73.5M 0s
    ## 125550K .......... .......... .......... .......... .......... 93% 75.1M 0s
    ## 125600K .......... .......... .......... .......... .......... 93% 77.1M 0s
    ## 125650K .......... .......... .......... .......... .......... 93% 23.9M 0s
    ## 125700K .......... .......... .......... .......... .......... 93% 51.1M 0s
    ## 125750K .......... .......... .......... .......... .......... 93% 69.8M 0s
    ## 125800K .......... .......... .......... .......... .......... 93% 77.4M 0s
    ## 125850K .......... .......... .......... .......... .......... 93% 93.0M 0s
    ## 125900K .......... .......... .......... .......... .......... 93% 82.2M 0s
    ## 125950K .......... .......... .......... .......... .......... 93% 84.7M 0s
    ## 126000K .......... .......... .......... .......... .......... 93% 67.0M 0s
    ## 126050K .......... .......... .......... .......... .......... 93% 59.6M 0s
    ## 126100K .......... .......... .......... .......... .......... 93% 70.8M 0s
    ## 126150K .......... .......... .......... .......... .......... 93% 71.5M 0s
    ## 126200K .......... .......... .......... .......... .......... 93% 42.2M 0s
    ## 126250K .......... .......... .......... .......... .......... 93% 82.2M 0s
    ## 126300K .......... .......... .......... .......... .......... 93% 86.4M 0s
    ## 126350K .......... .......... .......... .......... .......... 93% 16.9M 0s
    ## 126400K .......... .......... .......... .......... .......... 93% 85.2M 0s
    ## 126450K .......... .......... .......... .......... .......... 93% 88.9M 0s
    ## 126500K .......... .......... .......... .......... .......... 93% 80.5M 0s
    ## 126550K .......... .......... .......... .......... .......... 93% 89.8M 0s
    ## 126600K .......... .......... .......... .......... .......... 93% 84.0M 0s
    ## 126650K .......... .......... .......... .......... .......... 94% 20.7M 0s
    ## 126700K .......... .......... .......... .......... .......... 94% 60.3M 0s
    ## 126750K .......... .......... .......... .......... .......... 94% 67.5M 0s
    ## 126800K .......... .......... .......... .......... .......... 94% 75.5M 0s
    ## 126850K .......... .......... .......... .......... .......... 94% 80.5M 0s
    ## 126900K .......... .......... .......... .......... .......... 94% 73.0M 0s
    ## 126950K .......... .......... .......... .......... .......... 94% 5.55M 0s
    ## 127000K .......... .......... .......... .......... .......... 94% 40.5M 0s
    ## 127050K .......... .......... .......... .......... .......... 94% 64.9M 0s
    ## 127100K .......... .......... .......... .......... .......... 94% 39.0M 0s
    ## 127150K .......... .......... .......... .......... .......... 94% 61.3M 0s
    ## 127200K .......... .......... .......... .......... .......... 94% 77.0M 0s
    ## 127250K .......... .......... .......... .......... .......... 94% 79.9M 0s
    ## 127300K .......... .......... .......... .......... .......... 94% 73.9M 0s
    ## 127350K .......... .......... .......... .......... .......... 94% 76.4M 0s
    ## 127400K .......... .......... .......... .......... .......... 94% 37.7M 0s
    ## 127450K .......... .......... .......... .......... .......... 94% 39.1M 0s
    ## 127500K .......... .......... .......... .......... .......... 94% 55.4M 0s
    ## 127550K .......... .......... .......... .......... .......... 94% 52.4M 0s
    ## 127600K .......... .......... .......... .......... .......... 94% 65.3M 0s
    ## 127650K .......... .......... .......... .......... .......... 94% 80.7M 0s
    ## 127700K .......... .......... .......... .......... .......... 94% 72.5M 0s
    ## 127750K .......... .......... .......... .......... .......... 94% 83.6M 0s
    ## 127800K .......... .......... .......... .......... .......... 94% 55.7M 0s
    ## 127850K .......... .......... .......... .......... .......... 94% 68.0M 0s
    ## 127900K .......... .......... .......... .......... .......... 94% 73.3M 0s
    ## 127950K .......... .......... .......... .......... .......... 94% 80.4M 0s
    ## 128000K .......... .......... .......... .......... .......... 95% 62.1M 0s
    ## 128050K .......... .......... .......... .......... .......... 95% 71.4M 0s
    ## 128100K .......... .......... .......... .......... .......... 95% 69.2M 0s
    ## 128150K .......... .......... .......... .......... .......... 95% 68.3M 0s
    ## 128200K .......... .......... .......... .......... .......... 95% 63.9M 0s
    ## 128250K .......... .......... .......... .......... .......... 95% 70.3M 0s
    ## 128300K .......... .......... .......... .......... .......... 95% 75.8M 0s
    ## 128350K .......... .......... .......... .......... .......... 95% 73.2M 0s
    ## 128400K .......... .......... .......... .......... .......... 95% 71.9M 0s
    ## 128450K .......... .......... .......... .......... .......... 95% 70.3M 0s
    ## 128500K .......... .......... .......... .......... .......... 95% 67.4M 0s
    ## 128550K .......... .......... .......... .......... .......... 95% 74.2M 0s
    ## 128600K .......... .......... .......... .......... .......... 95% 53.2M 0s
    ## 128650K .......... .......... .......... .......... .......... 95% 55.9M 0s
    ## 128700K .......... .......... .......... .......... .......... 95% 79.2M 0s
    ## 128750K .......... .......... .......... .......... .......... 95% 58.5M 0s
    ## 128800K .......... .......... .......... .......... .......... 95% 90.8M 0s
    ## 128850K .......... .......... .......... .......... .......... 95% 86.5M 0s
    ## 128900K .......... .......... .......... .......... .......... 95% 14.2M 0s
    ## 128950K .......... .......... .......... .......... .......... 95% 85.3M 0s
    ## 129000K .......... .......... .......... .......... .......... 95% 82.9M 0s
    ## 129050K .......... .......... .......... .......... .......... 95% 82.7M 0s
    ## 129100K .......... .......... .......... .......... .......... 95% 87.2M 0s
    ## 129150K .......... .......... .......... .......... .......... 95% 86.4M 0s
    ## 129200K .......... .......... .......... .......... .......... 95% 56.0M 0s
    ## 129250K .......... .......... .......... .......... .......... 95% 53.6M 0s
    ## 129300K .......... .......... .......... .......... .......... 95% 39.1M 0s
    ## 129350K .......... .......... .......... .......... .......... 96% 51.2M 0s
    ## 129400K .......... .......... .......... .......... .......... 96% 72.5M 0s
    ## 129450K .......... .......... .......... .......... .......... 96% 83.7M 0s
    ## 129500K .......... .......... .......... .......... .......... 96% 76.4M 0s
    ## 129550K .......... .......... .......... .......... .......... 96% 72.6M 0s
    ## 129600K .......... .......... .......... .......... .......... 96% 68.5M 0s
    ## 129650K .......... .......... .......... .......... .......... 96% 46.3M 0s
    ## 129700K .......... .......... .......... .......... .......... 96% 30.8M 0s
    ## 129750K .......... .......... .......... .......... .......... 96% 78.9M 0s
    ## 129800K .......... .......... .......... .......... .......... 96% 69.8M 0s
    ## 129850K .......... .......... .......... .......... .......... 96% 92.7M 0s
    ## 129900K .......... .......... .......... .......... .......... 96% 86.1M 0s
    ## 129950K .......... .......... .......... .......... .......... 96% 31.2M 0s
    ## 130000K .......... .......... .......... .......... .......... 96% 59.3M 0s
    ## 130050K .......... .......... .......... .......... .......... 96% 32.5M 0s
    ## 130100K .......... .......... .......... .......... .......... 96% 82.2M 0s
    ## 130150K .......... .......... .......... .......... .......... 96% 91.2M 0s
    ## 130200K .......... .......... .......... .......... .......... 96% 85.7M 0s
    ## 130250K .......... .......... .......... .......... .......... 96% 85.1M 0s
    ## 130300K .......... .......... .......... .......... .......... 96% 69.0M 0s
    ## 130350K .......... .......... .......... .......... .......... 96% 87.3M 0s
    ## 130400K .......... .......... .......... .......... .......... 96% 64.5M 0s
    ## 130450K .......... .......... .......... .......... .......... 96% 40.6M 0s
    ## 130500K .......... .......... .......... .......... .......... 96% 55.0M 0s
    ## 130550K .......... .......... .......... .......... .......... 96% 71.9M 0s
    ## 130600K .......... .......... .......... .......... .......... 96% 56.8M 0s
    ## 130650K .......... .......... .......... .......... .......... 97% 75.2M 0s
    ## 130700K .......... .......... .......... .......... .......... 97% 69.9M 0s
    ## 130750K .......... .......... .......... .......... .......... 97% 50.4M 0s
    ## 130800K .......... .......... .......... .......... .......... 97% 68.0M 0s
    ## 130850K .......... .......... .......... .......... .......... 97% 60.3M 0s
    ## 130900K .......... .......... .......... .......... .......... 97% 65.5M 0s
    ## 130950K .......... .......... .......... .......... .......... 97% 93.0M 0s
    ## 131000K .......... .......... .......... .......... .......... 97% 91.1M 0s
    ## 131050K .......... .......... .......... .......... .......... 97% 81.1M 0s
    ## 131100K .......... .......... .......... .......... .......... 97% 77.2M 0s
    ## 131150K .......... .......... .......... .......... .......... 97% 69.3M 0s
    ## 131200K .......... .......... .......... .......... .......... 97% 69.8M 0s
    ## 131250K .......... .......... .......... .......... .......... 97% 76.1M 0s
    ## 131300K .......... .......... .......... .......... .......... 97% 85.0M 0s
    ## 131350K .......... .......... .......... .......... .......... 97% 24.3M 0s
    ## 131400K .......... .......... .......... .......... .......... 97% 70.9M 0s
    ## 131450K .......... .......... .......... .......... .......... 97% 70.1M 0s
    ## 131500K .......... .......... .......... .......... .......... 97% 64.3M 0s
    ## 131550K .......... .......... .......... .......... .......... 97% 87.2M 0s
    ## 131600K .......... .......... .......... .......... .......... 97% 86.9M 0s
    ## 131650K .......... .......... .......... .......... .......... 97% 74.1M 0s
    ## 131700K .......... .......... .......... .......... .......... 97% 72.2M 0s
    ## 131750K .......... .......... .......... .......... .......... 97% 69.9M 0s
    ## 131800K .......... .......... .......... .......... .......... 97% 73.7M 0s
    ## 131850K .......... .......... .......... .......... .......... 97% 27.6M 0s
    ## 131900K .......... .......... .......... .......... .......... 97% 76.8M 0s
    ## 131950K .......... .......... .......... .......... .......... 97% 81.5M 0s
    ## 132000K .......... .......... .......... .......... .......... 98% 73.2M 0s
    ## 132050K .......... .......... .......... .......... .......... 98% 65.5M 0s
    ## 132100K .......... .......... .......... .......... .......... 98% 66.5M 0s
    ## 132150K .......... .......... .......... .......... .......... 98% 86.2M 0s
    ## 132200K .......... .......... .......... .......... .......... 98% 88.6M 0s
    ## 132250K .......... .......... .......... .......... .......... 98% 54.0M 0s
    ## 132300K .......... .......... .......... .......... .......... 98% 67.3M 0s
    ## 132350K .......... .......... .......... .......... .......... 98% 70.3M 0s
    ## 132400K .......... .......... .......... .......... .......... 98% 75.9M 0s
    ## 132450K .......... .......... .......... .......... .......... 98% 52.3M 0s
    ## 132500K .......... .......... .......... .......... .......... 98% 74.8M 0s
    ## 132550K .......... .......... .......... .......... .......... 98% 71.9M 0s
    ## 132600K .......... .......... .......... .......... .......... 98% 52.0M 0s
    ## 132650K .......... .......... .......... .......... .......... 98% 79.5M 0s
    ## 132700K .......... .......... .......... .......... .......... 98% 64.4M 0s
    ## 132750K .......... .......... .......... .......... .......... 98% 47.8M 0s
    ## 132800K .......... .......... .......... .......... .......... 98% 63.3M 0s
    ## 132850K .......... .......... .......... .......... .......... 98% 80.6M 0s
    ## 132900K .......... .......... .......... .......... .......... 98% 60.1M 0s
    ## 132950K .......... .......... .......... .......... .......... 98% 67.9M 0s
    ## 133000K .......... .......... .......... .......... .......... 98% 66.7M 0s
    ## 133050K .......... .......... .......... .......... .......... 98% 61.5M 0s
    ## 133100K .......... .......... .......... .......... .......... 98% 70.5M 0s
    ## 133150K .......... .......... .......... .......... .......... 98% 73.2M 0s
    ## 133200K .......... .......... .......... .......... .......... 98% 68.2M 0s
    ## 133250K .......... .......... .......... .......... .......... 98% 66.0M 0s
    ## 133300K .......... .......... .......... .......... .......... 98% 65.9M 0s
    ## 133350K .......... .......... .......... .......... .......... 99% 62.6M 0s
    ## 133400K .......... .......... .......... .......... .......... 99% 73.8M 0s
    ## 133450K .......... .......... .......... .......... .......... 99% 66.9M 0s
    ## 133500K .......... .......... .......... .......... .......... 99% 65.5M 0s
    ## 133550K .......... .......... .......... .......... .......... 99% 76.8M 0s
    ## 133600K .......... .......... .......... .......... .......... 99% 56.7M 0s
    ## 133650K .......... .......... .......... .......... .......... 99% 61.2M 0s
    ## 133700K .......... .......... .......... .......... .......... 99% 70.2M 0s
    ## 133750K .......... .......... .......... .......... .......... 99% 74.0M 0s
    ## 133800K .......... .......... .......... .......... .......... 99% 60.1M 0s
    ## 133850K .......... .......... .......... .......... .......... 99% 63.2M 0s
    ## 133900K .......... .......... .......... .......... .......... 99% 67.2M 0s
    ## 133950K .......... .......... .......... .......... .......... 99% 55.7M 0s
    ## 134000K .......... .......... .......... .......... .......... 99% 63.7M 0s
    ## 134050K .......... .......... .......... .......... .......... 99% 69.9M 0s
    ## 134100K .......... .......... .......... .......... .......... 99% 58.2M 0s
    ## 134150K .......... .......... .......... .......... .......... 99% 77.8M 0s
    ## 134200K .......... .......... .......... .......... .......... 99% 68.4M 0s
    ## 134250K .......... .......... .......... .......... .......... 99% 67.9M 0s
    ## 134300K .......... .......... .......... .......... .......... 99% 69.8M 0s
    ## 134350K .......... .......... .......... .......... .......... 99% 95.2M 0s
    ## 134400K .......... .......... .......... .......... .......... 99% 71.1M 0s
    ## 134450K .......... .......... .......... .......... .......... 99% 77.9M 0s
    ## 134500K .......... .......... .......... .......... .......... 99% 66.3M 0s
    ## 134550K .......... .......... .......... .......... .......... 99% 71.5M 0s
    ## 134600K .......... .......... .......... .......... .......... 99% 85.3M 0s
    ## 134650K .......... .......... .......... .......... .......... 99% 93.9M 0s
    ## 134700K .......... .......... .......... ..........           100% 97.5M=2.6s
    ## 
    ## 2020-12-04 14:26:22 (50.2 MB/s) - ‘silva_nr99_v138_train_set.fa.gz.3’ saved [137973851/137973851]

Téléchargement des données sur Silva

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/Dada2_Phyloseq_CC1/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

Assignation des séquences étudiées à la classification. Cette
assignation se fait par rapport à une base de séquence de référence
ayant une classification connue.

## Téléchargement des données Silva permmetant l’assignement d’espèces aux des séquences analysées

``` bash
wget https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
```

    ## --2020-12-04 14:29:01--  https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 81840166 (78M) [application/octet-stream]
    ## Saving to: ‘silva_species_assignment_v138.fa.gz.3’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 14.9M 5s
    ##     50K .......... .......... .......... .......... ..........  0% 13.5M 6s
    ##    100K .......... .......... .......... .......... ..........  0% 15.9M 5s
    ##    150K .......... .......... .......... .......... ..........  0% 13.5M 5s
    ##    200K .......... .......... .......... .......... ..........  0% 26.6M 5s
    ##    250K .......... .......... .......... .......... ..........  0% 12.7M 5s
    ##    300K .......... .......... .......... .......... ..........  0% 33.9M 5s
    ##    350K .......... .......... .......... .......... ..........  0% 21.0M 5s
    ##    400K .......... .......... .......... .......... ..........  0% 41.4M 4s
    ##    450K .......... .......... .......... .......... ..........  0% 19.3M 4s
    ##    500K .......... .......... .......... .......... ..........  0% 49.3M 4s
    ##    550K .......... .......... .......... .......... ..........  0% 17.4M 4s
    ##    600K .......... .......... .......... .......... ..........  0% 63.4M 4s
    ##    650K .......... .......... .......... .......... ..........  0% 16.1M 4s
    ##    700K .......... .......... .......... .......... ..........  0%  129M 4s
    ##    750K .......... .......... .......... .......... ..........  1% 12.1M 4s
    ##    800K .......... .......... .......... .......... ..........  1%  141M 4s
    ##    850K .......... .......... .......... .......... ..........  1% 10.4M 4s
    ##    900K .......... .......... .......... .......... ..........  1% 15.7M 4s
    ##    950K .......... .......... .......... .......... ..........  1% 83.5M 4s
    ##   1000K .......... .......... .......... .......... ..........  1% 93.7M 4s
    ##   1050K .......... .......... .......... .......... ..........  1% 15.7M 4s
    ##   1100K .......... .......... .......... .......... ..........  1% 13.4M 4s
    ##   1150K .......... .......... .......... .......... ..........  1% 67.4M 4s
    ##   1200K .......... .......... .......... .......... ..........  1%  161M 4s
    ##   1250K .......... .......... .......... .......... ..........  1% 13.6M 4s
    ##   1300K .......... .......... .......... .......... ..........  1% 17.7M 4s
    ##   1350K .......... .......... .......... .......... ..........  1% 36.3M 4s
    ##   1400K .......... .......... .......... .......... ..........  1% 12.4M 4s
    ##   1450K .......... .......... .......... .......... ..........  1% 95.0M 4s
    ##   1500K .......... .......... .......... .......... ..........  1%  133M 3s
    ##   1550K .......... .......... .......... .......... ..........  2% 16.5M 3s
    ##   1600K .......... .......... .......... .......... ..........  2% 13.9M 4s
    ##   1650K .......... .......... .......... .......... ..........  2% 86.6M 3s
    ##   1700K .......... .......... .......... .......... ..........  2%  151M 3s
    ##   1750K .......... .......... .......... .......... ..........  2% 18.1M 3s
    ##   1800K .......... .......... .......... .......... ..........  2% 52.5M 3s
    ##   1850K .......... .......... .......... .......... ..........  2% 24.8M 3s
    ##   1900K .......... .......... .......... .......... ..........  2% 40.5M 3s
    ##   1950K .......... .......... .......... .......... ..........  2% 29.6M 3s
    ##   2000K .......... .......... .......... .......... ..........  2% 21.7M 3s
    ##   2050K .......... .......... .......... .......... ..........  2% 90.0M 3s
    ##   2100K .......... .......... .......... .......... ..........  2% 13.7M 3s
    ##   2150K .......... .......... .......... .......... ..........  2%  176M 3s
    ##   2200K .......... .......... .......... .......... ..........  2% 98.4M 3s
    ##   2250K .......... .......... .......... .......... ..........  2% 14.4M 3s
    ##   2300K .......... .......... .......... .......... ..........  2% 68.9M 3s
    ##   2350K .......... .......... .......... .......... ..........  3%  127M 3s
    ##   2400K .......... .......... .......... .......... ..........  3% 17.7M 3s
    ##   2450K .......... .......... .......... .......... ..........  3% 42.8M 3s
    ##   2500K .......... .......... .......... .......... ..........  3% 24.6M 3s
    ##   2550K .......... .......... .......... .......... ..........  3% 35.4M 3s
    ##   2600K .......... .......... .......... .......... ..........  3%  135M 3s
    ##   2650K .......... .......... .......... .......... ..........  3% 26.6M 3s
    ##   2700K .......... .......... .......... .......... ..........  3% 24.4M 3s
    ##   2750K .......... .......... .......... .......... ..........  3% 24.4M 3s
    ##   2800K .......... .......... .......... .......... ..........  3%  173M 3s
    ##   2850K .......... .......... .......... .......... ..........  3% 30.7M 3s
    ##   2900K .......... .......... .......... .......... ..........  3% 29.0M 3s
    ##   2950K .......... .......... .......... .......... ..........  3% 24.6M 3s
    ##   3000K .......... .......... .......... .......... ..........  3% 33.3M 3s
    ##   3050K .......... .......... .......... .......... ..........  3% 36.0M 3s
    ##   3100K .......... .......... .......... .......... ..........  3%  133M 3s
    ##   3150K .......... .......... .......... .......... ..........  4% 22.3M 3s
    ##   3200K .......... .......... .......... .......... ..........  4% 36.0M 3s
    ##   3250K .......... .......... .......... .......... ..........  4% 21.3M 3s
    ##   3300K .......... .......... .......... .......... ..........  4%  166M 3s
    ##   3350K .......... .......... .......... .......... ..........  4% 43.8M 3s
    ##   3400K .......... .......... .......... .......... ..........  4% 22.0M 3s
    ##   3450K .......... .......... .......... .......... ..........  4% 46.6M 3s
    ##   3500K .......... .......... .......... .......... ..........  4% 27.2M 3s
    ##   3550K .......... .......... .......... .......... ..........  4% 85.5M 3s
    ##   3600K .......... .......... .......... .......... ..........  4% 32.7M 3s
    ##   3650K .......... .......... .......... .......... ..........  4% 34.4M 3s
    ##   3700K .......... .......... .......... .......... ..........  4% 30.5M 3s
    ##   3750K .......... .......... .......... .......... ..........  4% 29.2M 3s
    ##   3800K .......... .......... .......... .......... ..........  4%  162M 3s
    ##   3850K .......... .......... .......... .......... ..........  4% 31.4M 3s
    ##   3900K .......... .......... .......... .......... ..........  4% 29.9M 3s
    ##   3950K .......... .......... .......... .......... ..........  5% 38.1M 3s
    ##   4000K .......... .......... .......... .......... ..........  5% 91.0M 3s
    ##   4050K .......... .......... .......... .......... ..........  5% 28.2M 3s
    ##   4100K .......... .......... .......... .......... ..........  5% 25.8M 3s
    ##   4150K .......... .......... .......... .......... ..........  5% 40.5M 3s
    ##   4200K .......... .......... .......... .......... ..........  5% 24.3M 3s
    ##   4250K .......... .......... .......... .......... ..........  5%  140M 3s
    ##   4300K .......... .......... .......... .......... ..........  5% 32.1M 3s
    ##   4350K .......... .......... .......... .......... ..........  5% 22.6M 3s
    ##   4400K .......... .......... .......... .......... ..........  5% 52.0M 3s
    ##   4450K .......... .......... .......... .......... ..........  5%  153M 3s
    ##   4500K .......... .......... .......... .......... ..........  5% 23.7M 3s
    ##   4550K .......... .......... .......... .......... ..........  5% 78.8M 3s
    ##   4600K .......... .......... .......... .......... ..........  5% 56.7M 3s
    ##   4650K .......... .......... .......... .......... ..........  5% 21.8M 3s
    ##   4700K .......... .......... .......... .......... ..........  5% 74.1M 3s
    ##   4750K .......... .......... .......... .......... ..........  6%  136M 3s
    ##   4800K .......... .......... .......... .......... ..........  6% 61.6M 3s
    ##   4850K .......... .......... .......... .......... ..........  6% 29.6M 3s
    ##   4900K .......... .......... .......... .......... ..........  6% 28.2M 3s
    ##   4950K .......... .......... .......... .......... ..........  6%  165M 2s
    ##   5000K .......... .......... .......... .......... ..........  6%  189M 2s
    ##   5050K .......... .......... .......... .......... ..........  6% 34.2M 2s
    ##   5100K .......... .......... .......... .......... ..........  6% 30.0M 2s
    ##   5150K .......... .......... .......... .......... ..........  6% 55.6M 2s
    ##   5200K .......... .......... .......... .......... ..........  6%  117M 2s
    ##   5250K .......... .......... .......... .......... ..........  6% 57.4M 2s
    ##   5300K .......... .......... .......... .......... ..........  6% 34.8M 2s
    ##   5350K .......... .......... .......... .......... ..........  6% 29.4M 2s
    ##   5400K .......... .......... .......... .......... ..........  6% 27.8M 2s
    ##   5450K .......... .......... .......... .......... ..........  6%  123M 2s
    ##   5500K .......... .......... .......... .......... ..........  6%  140M 2s
    ##   5550K .......... .......... .......... .......... ..........  7% 62.8M 2s
    ##   5600K .......... .......... .......... .......... ..........  7% 21.7M 2s
    ##   5650K .......... .......... .......... .......... ..........  7% 71.7M 2s
    ##   5700K .......... .......... .......... .......... ..........  7% 58.9M 2s
    ##   5750K .......... .......... .......... .......... ..........  7% 79.1M 2s
    ##   5800K .......... .......... .......... .......... ..........  7% 52.1M 2s
    ##   5850K .......... .......... .......... .......... ..........  7% 26.8M 2s
    ##   5900K .......... .......... .......... .......... ..........  7%  135M 2s
    ##   5950K .......... .......... .......... .......... ..........  7% 31.9M 2s
    ##   6000K .......... .......... .......... .......... ..........  7% 52.8M 2s
    ##   6050K .......... .......... .......... .......... ..........  7% 45.9M 2s
    ##   6100K .......... .......... .......... .......... ..........  7%  121M 2s
    ##   6150K .......... .......... .......... .......... ..........  7% 34.2M 2s
    ##   6200K .......... .......... .......... .......... ..........  7% 45.1M 2s
    ##   6250K .......... .......... .......... .......... ..........  7% 51.9M 2s
    ##   6300K .......... .......... .......... .......... ..........  7% 34.8M 2s
    ##   6350K .......... .......... .......... .......... ..........  8%  155M 2s
    ##   6400K .......... .......... .......... .......... ..........  8% 56.6M 2s
    ##   6450K .......... .......... .......... .......... ..........  8% 47.5M 2s
    ##   6500K .......... .......... .......... .......... ..........  8% 33.7M 2s
    ##   6550K .......... .......... .......... .......... ..........  8% 62.4M 2s
    ##   6600K .......... .......... .......... .......... ..........  8% 72.3M 2s
    ##   6650K .......... .......... .......... .......... ..........  8% 55.1M 2s
    ##   6700K .......... .......... .......... .......... ..........  8% 33.9M 2s
    ##   6750K .......... .......... .......... .......... ..........  8% 42.8M 2s
    ##   6800K .......... .......... .......... .......... ..........  8% 64.1M 2s
    ##   6850K .......... .......... .......... .......... ..........  8% 91.5M 2s
    ##   6900K .......... .......... .......... .......... ..........  8% 36.9M 2s
    ##   6950K .......... .......... .......... .......... ..........  8% 38.6M 2s
    ##   7000K .......... .......... .......... .......... ..........  8% 50.9M 2s
    ##   7050K .......... .......... .......... .......... ..........  8% 43.2M 2s
    ##   7100K .......... .......... .......... .......... ..........  8%  151M 2s
    ##   7150K .......... .......... .......... .......... ..........  9% 39.1M 2s
    ##   7200K .......... .......... .......... .......... ..........  9% 41.0M 2s
    ##   7250K .......... .......... .......... .......... ..........  9% 30.2M 2s
    ##   7300K .......... .......... .......... .......... ..........  9%  161M 2s
    ##   7350K .......... .......... .......... .......... ..........  9% 66.9M 2s
    ##   7400K .......... .......... .......... .......... ..........  9% 50.4M 2s
    ##   7450K .......... .......... .......... .......... ..........  9% 31.2M 2s
    ##   7500K .......... .......... .......... .......... ..........  9% 45.4M 2s
    ##   7550K .......... .......... .......... .......... ..........  9%  149M 2s
    ##   7600K .......... .......... .......... .......... ..........  9% 5.44M 2s
    ##   7650K .......... .......... .......... .......... ..........  9% 49.4M 2s
    ##   7700K .......... .......... .......... .......... ..........  9% 76.1M 2s
    ##   7750K .......... .......... .......... .......... ..........  9% 51.6M 2s
    ##   7800K .......... .......... .......... .......... ..........  9% 59.6M 2s
    ##   7850K .......... .......... .......... .......... ..........  9% 29.3M 2s
    ##   7900K .......... .......... .......... .......... ..........  9% 84.7M 2s
    ##   7950K .......... .......... .......... .......... .......... 10% 44.2M 2s
    ##   8000K .......... .......... .......... .......... .......... 10%  136M 2s
    ##   8050K .......... .......... .......... .......... .......... 10% 6.04M 2s
    ##   8100K .......... .......... .......... .......... .......... 10% 62.7M 2s
    ##   8150K .......... .......... .......... .......... .......... 10% 98.7M 2s
    ##   8200K .......... .......... .......... .......... .......... 10%  168M 2s
    ##   8250K .......... .......... .......... .......... .......... 10% 27.7M 2s
    ##   8300K .......... .......... .......... .......... .......... 10% 28.5M 2s
    ##   8350K .......... .......... .......... .......... .......... 10%  139M 2s
    ##   8400K .......... .......... .......... .......... .......... 10% 43.0M 2s
    ##   8450K .......... .......... .......... .......... .......... 10% 31.6M 2s
    ##   8500K .......... .......... .......... .......... .......... 10%  134M 2s
    ##   8550K .......... .......... .......... .......... .......... 10% 57.4M 2s
    ##   8600K .......... .......... .......... .......... .......... 10% 59.5M 2s
    ##   8650K .......... .......... .......... .......... .......... 10% 32.1M 2s
    ##   8700K .......... .......... .......... .......... .......... 10%  117M 2s
    ##   8750K .......... .......... .......... .......... .......... 11% 56.4M 2s
    ##   8800K .......... .......... .......... .......... .......... 11% 86.9M 2s
    ##   8850K .......... .......... .......... .......... .......... 11% 26.7M 2s
    ##   8900K .......... .......... .......... .......... .......... 11% 56.7M 2s
    ##   8950K .......... .......... .......... .......... .......... 11%  131M 2s
    ##   9000K .......... .......... .......... .......... .......... 11% 69.1M 2s
    ##   9050K .......... .......... .......... .......... .......... 11% 28.3M 2s
    ##   9100K .......... .......... .......... .......... .......... 11% 47.0M 2s
    ##   9150K .......... .......... .......... .......... .......... 11% 83.9M 2s
    ##   9200K .......... .......... .......... .......... .......... 11%  123M 2s
    ##   9250K .......... .......... .......... .......... .......... 11% 32.2M 2s
    ##   9300K .......... .......... .......... .......... .......... 11% 22.5M 2s
    ##   9350K .......... .......... .......... .......... .......... 11%  131M 2s
    ##   9400K .......... .......... .......... .......... .......... 11%  132M 2s
    ##   9450K .......... .......... .......... .......... .......... 11%  154M 2s
    ##   9500K .......... .......... .......... .......... .......... 11% 16.4M 2s
    ##   9550K .......... .......... .......... .......... .......... 12%  186M 2s
    ##   9600K .......... .......... .......... .......... .......... 12%  137M 2s
    ##   9650K .......... .......... .......... .......... .......... 12%  201M 2s
    ##   9700K .......... .......... .......... .......... .......... 12% 17.4M 2s
    ##   9750K .......... .......... .......... .......... .......... 12%  161M 2s
    ##   9800K .......... .......... .......... .......... .......... 12%  207M 2s
    ##   9850K .......... .......... .......... .......... .......... 12% 93.3M 2s
    ##   9900K .......... .......... .......... .......... .......... 12% 9.82M 2s
    ##   9950K .......... .......... .......... .......... .......... 12%  164M 2s
    ##  10000K .......... .......... .......... .......... .......... 12%  140M 2s
    ##  10050K .......... .......... .......... .......... .......... 12% 18.9M 2s
    ##  10100K .......... .......... .......... .......... .......... 12%  139M 2s
    ##  10150K .......... .......... .......... .......... .......... 12%  140M 2s
    ##  10200K .......... .......... .......... .......... .......... 12% 89.4M 2s
    ##  10250K .......... .......... .......... .......... .......... 12% 20.1M 2s
    ##  10300K .......... .......... .......... .......... .......... 12% 7.53M 2s
    ##  10350K .......... .......... .......... .......... .......... 13%  136M 2s
    ##  10400K .......... .......... .......... .......... .......... 13% 51.2M 2s
    ##  10450K .......... .......... .......... .......... .......... 13%  127M 2s
    ##  10500K .......... .......... .......... .......... .......... 13% 30.2M 2s
    ##  10550K .......... .......... .......... .......... .......... 13% 22.0M 2s
    ##  10600K .......... .......... .......... .......... .......... 13%  221M 2s
    ##  10650K .......... .......... .......... .......... .......... 13%  191M 2s
    ##  10700K .......... .......... .......... .......... .......... 13% 72.7M 2s
    ##  10750K .......... .......... .......... .......... .......... 13% 17.1M 2s
    ##  10800K .......... .......... .......... .......... .......... 13% 62.2M 2s
    ##  10850K .......... .......... .......... .......... .......... 13%  169M 2s
    ##  10900K .......... .......... .......... .......... .......... 13% 6.27M 2s
    ##  10950K .......... .......... .......... .......... .......... 13% 66.1M 2s
    ##  11000K .......... .......... .......... .......... .......... 13% 77.8M 2s
    ##  11050K .......... .......... .......... .......... .......... 13%  106M 2s
    ##  11100K .......... .......... .......... .......... .......... 13% 27.1M 2s
    ##  11150K .......... .......... .......... .......... .......... 14% 88.2M 2s
    ##  11200K .......... .......... .......... .......... .......... 14% 30.0M 2s
    ##  11250K .......... .......... .......... .......... .......... 14% 36.9M 2s
    ##  11300K .......... .......... .......... .......... .......... 14%  172M 2s
    ##  11350K .......... .......... .......... .......... .......... 14%  140M 2s
    ##  11400K .......... .......... .......... .......... .......... 14% 33.0M 2s
    ##  11450K .......... .......... .......... .......... .......... 14% 28.4M 2s
    ##  11500K .......... .......... .......... .......... .......... 14%  105M 2s
    ##  11550K .......... .......... .......... .......... .......... 14%  137M 2s
    ##  11600K .......... .......... .......... .......... .......... 14% 28.5M 2s
    ##  11650K .......... .......... .......... .......... .......... 14% 41.1M 2s
    ##  11700K .......... .......... .......... .......... .......... 14% 76.2M 2s
    ##  11750K .......... .......... .......... .......... .......... 14% 39.2M 2s
    ##  11800K .......... .......... .......... .......... .......... 14% 71.2M 2s
    ##  11850K .......... .......... .......... .......... .......... 14% 33.7M 2s
    ##  11900K .......... .......... .......... .......... .......... 14% 74.6M 2s
    ##  11950K .......... .......... .......... .......... .......... 15% 30.9M 2s
    ##  12000K .......... .......... .......... .......... .......... 15%  182M 2s
    ##  12050K .......... .......... .......... .......... .......... 15% 22.8M 2s
    ##  12100K .......... .......... .......... .......... .......... 15%  168M 2s
    ##  12150K .......... .......... .......... .......... .......... 15% 61.9M 2s
    ##  12200K .......... .......... .......... .......... .......... 15% 24.8M 2s
    ##  12250K .......... .......... .......... .......... .......... 15% 42.5M 2s
    ##  12300K .......... .......... .......... .......... .......... 15%  151M 2s
    ##  12350K .......... .......... .......... .......... .......... 15% 69.0M 2s
    ##  12400K .......... .......... .......... .......... .......... 15% 51.4M 2s
    ##  12450K .......... .......... .......... .......... .......... 15% 29.8M 2s
    ##  12500K .......... .......... .......... .......... .......... 15%  125M 2s
    ##  12550K .......... .......... .......... .......... .......... 15% 58.1M 2s
    ##  12600K .......... .......... .......... .......... .......... 15% 61.7M 2s
    ##  12650K .......... .......... .......... .......... .......... 15% 42.2M 2s
    ##  12700K .......... .......... .......... .......... .......... 15% 45.8M 2s
    ##  12750K .......... .......... .......... .......... .......... 16% 92.5M 2s
    ##  12800K .......... .......... .......... .......... .......... 16% 39.6M 2s
    ##  12850K .......... .......... .......... .......... .......... 16%  102M 2s
    ##  12900K .......... .......... .......... .......... .......... 16% 24.9M 2s
    ##  12950K .......... .......... .......... .......... .......... 16%  173M 2s
    ##  13000K .......... .......... .......... .......... .......... 16% 23.9M 2s
    ##  13050K .......... .......... .......... .......... .......... 16%  162M 2s
    ##  13100K .......... .......... .......... .......... .......... 16% 51.1M 2s
    ##  13150K .......... .......... .......... .......... .......... 16% 16.3M 2s
    ##  13200K .......... .......... .......... .......... .......... 16%  140M 2s
    ##  13250K .......... .......... .......... .......... .......... 16%  189M 2s
    ##  13300K .......... .......... .......... .......... .......... 16% 48.8M 2s
    ##  13350K .......... .......... .......... .......... .......... 16% 24.4M 2s
    ##  13400K .......... .......... .......... .......... .......... 16%  138M 2s
    ##  13450K .......... .......... .......... .......... .......... 16%  171M 2s
    ##  13500K .......... .......... .......... .......... .......... 16% 50.5M 2s
    ##  13550K .......... .......... .......... .......... .......... 17% 28.2M 2s
    ##  13600K .......... .......... .......... .......... .......... 17% 54.3M 2s
    ##  13650K .......... .......... .......... .......... .......... 17% 66.0M 2s
    ##  13700K .......... .......... .......... .......... .......... 17% 88.0M 2s
    ##  13750K .......... .......... .......... .......... .......... 17% 25.6M 2s
    ##  13800K .......... .......... .......... .......... .......... 17% 60.6M 2s
    ##  13850K .......... .......... .......... .......... .......... 17% 43.3M 2s
    ##  13900K .......... .......... .......... .......... .......... 17%  109M 2s
    ##  13950K .......... .......... .......... .......... .......... 17% 22.0M 2s
    ##  14000K .......... .......... .......... .......... .......... 17%  163M 2s
    ##  14050K .......... .......... .......... .......... .......... 17% 45.4M 2s
    ##  14100K .......... .......... .......... .......... .......... 17% 23.8M 2s
    ##  14150K .......... .......... .......... .......... .......... 17% 58.7M 2s
    ##  14200K .......... .......... .......... .......... .......... 17% 58.3M 2s
    ##  14250K .......... .......... .......... .......... .......... 17% 63.5M 2s
    ##  14300K .......... .......... .......... .......... .......... 17% 36.7M 2s
    ##  14350K .......... .......... .......... .......... .......... 18%  173M 2s
    ##  14400K .......... .......... .......... .......... .......... 18% 38.4M 2s
    ##  14450K .......... .......... .......... .......... .......... 18% 55.3M 2s
    ##  14500K .......... .......... .......... .......... .......... 18% 39.9M 2s
    ##  14550K .......... .......... .......... .......... .......... 18% 11.9M 2s
    ##  14600K .......... .......... .......... .......... .......... 18%  122M 2s
    ##  14650K .......... .......... .......... .......... .......... 18% 43.1M 2s
    ##  14700K .......... .......... .......... .......... .......... 18%  125M 2s
    ##  14750K .......... .......... .......... .......... .......... 18% 25.2M 2s
    ##  14800K .......... .......... .......... .......... .......... 18% 24.7M 2s
    ##  14850K .......... .......... .......... .......... .......... 18%  140M 2s
    ##  14900K .......... .......... .......... .......... .......... 18%  122M 2s
    ##  14950K .......... .......... .......... .......... .......... 18% 6.76M 2s
    ##  15000K .......... .......... .......... .......... .......... 18% 62.4M 2s
    ##  15050K .......... .......... .......... .......... .......... 18% 69.6M 2s
    ##  15100K .......... .......... .......... .......... .......... 18%  158M 2s
    ##  15150K .......... .......... .......... .......... .......... 19% 33.5M 2s
    ##  15200K .......... .......... .......... .......... .......... 19% 28.9M 2s
    ##  15250K .......... .......... .......... .......... .......... 19%  155M 2s
    ##  15300K .......... .......... .......... .......... .......... 19%  143M 2s
    ##  15350K .......... .......... .......... .......... .......... 19%  120M 2s
    ##  15400K .......... .......... .......... .......... .......... 19% 18.0M 2s
    ##  15450K .......... .......... .......... .......... .......... 19% 45.2M 2s
    ##  15500K .......... .......... .......... .......... .......... 19%  132M 2s
    ##  15550K .......... .......... .......... .......... .......... 19%  129M 2s
    ##  15600K .......... .......... .......... .......... .......... 19% 28.1M 2s
    ##  15650K .......... .......... .......... .......... .......... 19% 28.2M 2s
    ##  15700K .......... .......... .......... .......... .......... 19%  132M 2s
    ##  15750K .......... .......... .......... .......... .......... 19% 74.1M 2s
    ##  15800K .......... .......... .......... .......... .......... 19%  132M 2s
    ##  15850K .......... .......... .......... .......... .......... 19% 20.1M 2s
    ##  15900K .......... .......... .......... .......... .......... 19% 85.3M 2s
    ##  15950K .......... .......... .......... .......... .......... 20% 57.6M 2s
    ##  16000K .......... .......... .......... .......... .......... 20% 26.0M 2s
    ##  16050K .......... .......... .......... .......... .......... 20%  135M 2s
    ##  16100K .......... .......... .......... .......... .......... 20% 35.6M 2s
    ##  16150K .......... .......... .......... .......... .......... 20%  110M 2s
    ##  16200K .......... .......... .......... .......... .......... 20% 29.7M 2s
    ##  16250K .......... .......... .......... .......... .......... 20% 10.4M 2s
    ##  16300K .......... .......... .......... .......... .......... 20%  138M 2s
    ##  16350K .......... .......... .......... .......... .......... 20%  107M 2s
    ##  16400K .......... .......... .......... .......... .......... 20%  122M 2s
    ##  16450K .......... .......... .......... .......... .......... 20% 22.3M 2s
    ##  16500K .......... .......... .......... .......... .......... 20% 85.0M 2s
    ##  16550K .......... .......... .......... .......... .......... 20% 61.6M 2s
    ##  16600K .......... .......... .......... .......... .......... 20% 69.8M 2s
    ##  16650K .......... .......... .......... .......... .......... 20% 34.8M 2s
    ##  16700K .......... .......... .......... .......... .......... 20%  161M 2s
    ##  16750K .......... .......... .......... .......... .......... 21% 40.4M 2s
    ##  16800K .......... .......... .......... .......... .......... 21% 15.0M 2s
    ##  16850K .......... .......... .......... .......... .......... 21% 49.3M 2s
    ##  16900K .......... .......... .......... .......... .......... 21% 65.1M 2s
    ##  16950K .......... .......... .......... .......... .......... 21% 25.7M 2s
    ##  17000K .......... .......... .......... .......... .......... 21%  165M 2s
    ##  17050K .......... .......... .......... .......... .......... 21% 38.9M 2s
    ##  17100K .......... .......... .......... .......... .......... 21%  133M 2s
    ##  17150K .......... .......... .......... .......... .......... 21% 46.7M 2s
    ##  17200K .......... .......... .......... .......... .......... 21% 5.89M 2s
    ##  17250K .......... .......... .......... .......... .......... 21%  114M 2s
    ##  17300K .......... .......... .......... .......... .......... 21%  156M 2s
    ##  17350K .......... .......... .......... .......... .......... 21%  163M 2s
    ##  17400K .......... .......... .......... .......... .......... 21%  138M 2s
    ##  17450K .......... .......... .......... .......... .......... 21% 20.6M 2s
    ##  17500K .......... .......... .......... .......... .......... 21% 49.0M 2s
    ##  17550K .......... .......... .......... .......... .......... 22%  111M 2s
    ##  17600K .......... .......... .......... .......... .......... 22%  103M 2s
    ##  17650K .......... .......... .......... .......... .......... 22%  150M 2s
    ##  17700K .......... .......... .......... .......... .......... 22% 35.7M 2s
    ##  17750K .......... .......... .......... .......... .......... 22% 42.2M 2s
    ##  17800K .......... .......... .......... .......... .......... 22% 53.1M 2s
    ##  17850K .......... .......... .......... .......... .......... 22% 41.0M 2s
    ##  17900K .......... .......... .......... .......... .......... 22%  179M 2s
    ##  17950K .......... .......... .......... .......... .......... 22% 37.9M 2s
    ##  18000K .......... .......... .......... .......... .......... 22% 31.8M 2s
    ##  18050K .......... .......... .......... .......... .......... 22% 63.4M 2s
    ##  18100K .......... .......... .......... .......... .......... 22% 49.1M 2s
    ##  18150K .......... .......... .......... .......... .......... 22%  154M 2s
    ##  18200K .......... .......... .......... .......... .......... 22% 29.4M 2s
    ##  18250K .......... .......... .......... .......... .......... 22% 83.7M 2s
    ##  18300K .......... .......... .......... .......... .......... 22% 32.7M 2s
    ##  18350K .......... .......... .......... .......... .......... 23% 44.5M 2s
    ##  18400K .......... .......... .......... .......... .......... 23% 78.4M 2s
    ##  18450K .......... .......... .......... .......... .......... 23% 54.0M 2s
    ##  18500K .......... .......... .......... .......... .......... 23% 8.24M 2s
    ##  18550K .......... .......... .......... .......... .......... 23% 56.2M 2s
    ##  18600K .......... .......... .......... .......... .......... 23% 91.1M 2s
    ##  18650K .......... .......... .......... .......... .......... 23%  161M 2s
    ##  18700K .......... .......... .......... .......... .......... 23% 25.1M 2s
    ##  18750K .......... .......... .......... .......... .......... 23% 50.7M 2s
    ##  18800K .......... .......... .......... .......... .......... 23% 71.8M 2s
    ##  18850K .......... .......... .......... .......... .......... 23%  145M 2s
    ##  18900K .......... .......... .......... .......... .......... 23% 31.0M 2s
    ##  18950K .......... .......... .......... .......... .......... 23% 33.2M 2s
    ##  19000K .......... .......... .......... .......... .......... 23% 99.2M 2s
    ##  19050K .......... .......... .......... .......... .......... 23% 26.2M 2s
    ##  19100K .......... .......... .......... .......... .......... 23%  135M 2s
    ##  19150K .......... .......... .......... .......... .......... 24% 36.4M 2s
    ##  19200K .......... .......... .......... .......... .......... 24%  134M 2s
    ##  19250K .......... .......... .......... .......... .......... 24% 25.5M 2s
    ##  19300K .......... .......... .......... .......... .......... 24%  203M 2s
    ##  19350K .......... .......... .......... .......... .......... 24% 32.3M 2s
    ##  19400K .......... .......... .......... .......... .......... 24%  109M 2s
    ##  19450K .......... .......... .......... .......... .......... 24% 21.8M 2s
    ##  19500K .......... .......... .......... .......... .......... 24% 88.6M 2s
    ##  19550K .......... .......... .......... .......... .......... 24%  107M 2s
    ##  19600K .......... .......... .......... .......... .......... 24% 75.5M 2s
    ##  19650K .......... .......... .......... .......... .......... 24% 38.4M 2s
    ##  19700K .......... .......... .......... .......... .......... 24% 29.1M 2s
    ##  19750K .......... .......... .......... .......... .......... 24% 61.0M 2s
    ##  19800K .......... .......... .......... .......... .......... 24% 15.6M 2s
    ##  19850K .......... .......... .......... .......... .......... 24%  132M 2s
    ##  19900K .......... .......... .......... .......... .......... 24% 38.5M 2s
    ##  19950K .......... .......... .......... .......... .......... 25%  159M 2s
    ##  20000K .......... .......... .......... .......... .......... 25%  185M 2s
    ##  20050K .......... .......... .......... .......... .......... 25% 30.6M 2s
    ##  20100K .......... .......... .......... .......... .......... 25% 27.2M 2s
    ##  20150K .......... .......... .......... .......... .......... 25%  159M 2s
    ##  20200K .......... .......... .......... .......... .......... 25% 59.9M 2s
    ##  20250K .......... .......... .......... .......... .......... 25% 20.2M 2s
    ##  20300K .......... .......... .......... .......... .......... 25% 98.2M 2s
    ##  20350K .......... .......... .......... .......... .......... 25%  140M 2s
    ##  20400K .......... .......... .......... .......... .......... 25% 57.6M 2s
    ##  20450K .......... .......... .......... .......... .......... 25% 25.1M 2s
    ##  20500K .......... .......... .......... .......... .......... 25%  109M 2s
    ##  20550K .......... .......... .......... .......... .......... 25% 86.3M 2s
    ##  20600K .......... .......... .......... .......... .......... 25% 58.0M 2s
    ##  20650K .......... .......... .......... .......... .......... 25% 33.7M 2s
    ##  20700K .......... .......... .......... .......... .......... 25% 21.5M 2s
    ##  20750K .......... .......... .......... .......... .......... 26%  139M 2s
    ##  20800K .......... .......... .......... .......... .......... 26%  194M 2s
    ##  20850K .......... .......... .......... .......... .......... 26% 53.6M 2s
    ##  20900K .......... .......... .......... .......... .......... 26% 10.3M 2s
    ##  20950K .......... .......... .......... .......... .......... 26%  118M 2s
    ##  21000K .......... .......... .......... .......... .......... 26% 59.8M 2s
    ##  21050K .......... .......... .......... .......... .......... 26%  118M 2s
    ##  21100K .......... .......... .......... .......... .......... 26% 27.4M 2s
    ##  21150K .......... .......... .......... .......... .......... 26% 47.9M 2s
    ##  21200K .......... .......... .......... .......... .......... 26%  135M 2s
    ##  21250K .......... .......... .......... .......... .......... 26% 60.3M 2s
    ##  21300K .......... .......... .......... .......... .......... 26% 33.6M 1s
    ##  21350K .......... .......... .......... .......... .......... 26% 24.6M 2s
    ##  21400K .......... .......... .......... .......... .......... 26%  165M 1s
    ##  21450K .......... .......... .......... .......... .......... 26%  149M 1s
    ##  21500K .......... .......... .......... .......... .......... 26%  156M 1s
    ##  21550K .......... .......... .......... .......... .......... 27% 6.83M 2s
    ##  21600K .......... .......... .......... .......... .......... 27% 59.8M 2s
    ##  21650K .......... .......... .......... .......... .......... 27% 77.1M 1s
    ##  21700K .......... .......... .......... .......... .......... 27%  136M 1s
    ##  21750K .......... .......... .......... .......... .......... 27% 29.9M 1s
    ##  21800K .......... .......... .......... .......... .......... 27% 22.8M 1s
    ##  21850K .......... .......... .......... .......... .......... 27% 86.4M 1s
    ##  21900K .......... .......... .......... .......... .......... 27% 32.0M 1s
    ##  21950K .......... .......... .......... .......... .......... 27%  160M 1s
    ##  22000K .......... .......... .......... .......... .......... 27% 79.0M 1s
    ##  22050K .......... .......... .......... .......... .......... 27% 28.6M 1s
    ##  22100K .......... .......... .......... .......... .......... 27% 30.9M 1s
    ##  22150K .......... .......... .......... .......... .......... 27% 7.93M 1s
    ##  22200K .......... .......... .......... .......... .......... 27%  161M 1s
    ##  22250K .......... .......... .......... .......... .......... 27%  161M 1s
    ##  22300K .......... .......... .......... .......... .......... 27%  192M 1s
    ##  22350K .......... .......... .......... .......... .......... 28%  164M 1s
    ##  22400K .......... .......... .......... .......... .......... 28% 19.6M 1s
    ##  22450K .......... .......... .......... .......... .......... 28% 55.2M 1s
    ##  22500K .......... .......... .......... .......... .......... 28% 77.8M 1s
    ##  22550K .......... .......... .......... .......... .......... 28% 85.3M 1s
    ##  22600K .......... .......... .......... .......... .......... 28%  125M 1s
    ##  22650K .......... .......... .......... .......... .......... 28%  158M 1s
    ##  22700K .......... .......... .......... .......... .......... 28% 10.7M 1s
    ##  22750K .......... .......... .......... .......... .......... 28%  112M 1s
    ##  22800K .......... .......... .......... .......... .......... 28%  127M 1s
    ##  22850K .......... .......... .......... .......... .......... 28%  160M 1s
    ##  22900K .......... .......... .......... .......... .......... 28%  190M 1s
    ##  22950K .......... .......... .......... .......... .......... 28% 14.7M 1s
    ##  23000K .......... .......... .......... .......... .......... 28% 15.2M 1s
    ##  23050K .......... .......... .......... .......... .......... 28%  152M 1s
    ##  23100K .......... .......... .......... .......... .......... 28% 96.2M 1s
    ##  23150K .......... .......... .......... .......... .......... 29%  127M 1s
    ##  23200K .......... .......... .......... .......... .......... 29% 34.9M 1s
    ##  23250K .......... .......... .......... .......... .......... 29% 51.7M 1s
    ##  23300K .......... .......... .......... .......... .......... 29% 44.3M 1s
    ##  23350K .......... .......... .......... .......... .......... 29% 79.0M 1s
    ##  23400K .......... .......... .......... .......... .......... 29%  118M 1s
    ##  23450K .......... .......... .......... .......... .......... 29% 48.8M 1s
    ##  23500K .......... .......... .......... .......... .......... 29% 56.7M 1s
    ##  23550K .......... .......... .......... .......... .......... 29%  130M 1s
    ##  23600K .......... .......... .......... .......... .......... 29% 39.1M 1s
    ##  23650K .......... .......... .......... .......... .......... 29% 71.4M 1s
    ##  23700K .......... .......... .......... .......... .......... 29% 47.2M 1s
    ##  23750K .......... .......... .......... .......... .......... 29% 60.8M 1s
    ##  23800K .......... .......... .......... .......... .......... 29% 48.1M 1s
    ##  23850K .......... .......... .......... .......... .......... 29% 93.5M 1s
    ##  23900K .......... .......... .......... .......... .......... 29% 7.29M 1s
    ##  23950K .......... .......... .......... .......... .......... 30% 43.7M 1s
    ##  24000K .......... .......... .......... .......... .......... 30%  131M 1s
    ##  24050K .......... .......... .......... .......... .......... 30% 54.9M 1s
    ##  24100K .......... .......... .......... .......... .......... 30%  152M 1s
    ##  24150K .......... .......... .......... .......... .......... 30% 41.2M 1s
    ##  24200K .......... .......... .......... .......... .......... 30%  103M 1s
    ##  24250K .......... .......... .......... .......... .......... 30% 34.4M 1s
    ##  24300K .......... .......... .......... .......... .......... 30% 65.2M 1s
    ##  24350K .......... .......... .......... .......... .......... 30%  128M 1s
    ##  24400K .......... .......... .......... .......... .......... 30% 60.6M 1s
    ##  24450K .......... .......... .......... .......... .......... 30%  109M 1s
    ##  24500K .......... .......... .......... .......... .......... 30% 14.2M 1s
    ##  24550K .......... .......... .......... .......... .......... 30% 33.0M 1s
    ##  24600K .......... .......... .......... .......... .......... 30%  120M 1s
    ##  24650K .......... .......... .......... .......... .......... 30%  163M 1s
    ##  24700K .......... .......... .......... .......... .......... 30% 86.8M 1s
    ##  24750K .......... .......... .......... .......... .......... 31% 66.3M 1s
    ##  24800K .......... .......... .......... .......... .......... 31% 29.0M 1s
    ##  24850K .......... .......... .......... .......... .......... 31% 5.00M 1s
    ##  24900K .......... .......... .......... .......... .......... 31% 57.1M 1s
    ##  24950K .......... .......... .......... .......... .......... 31% 49.8M 1s
    ##  25000K .......... .......... .......... .......... .......... 31%  140M 1s
    ##  25050K .......... .......... .......... .......... .......... 31%  158M 1s
    ##  25100K .......... .......... .......... .......... .......... 31% 48.4M 1s
    ##  25150K .......... .......... .......... .......... .......... 31% 24.5M 1s
    ##  25200K .......... .......... .......... .......... .......... 31% 48.1M 1s
    ##  25250K .......... .......... .......... .......... .......... 31% 66.8M 1s
    ##  25300K .......... .......... .......... .......... .......... 31%  138M 1s
    ##  25350K .......... .......... .......... .......... .......... 31% 95.7M 1s
    ##  25400K .......... .......... .......... .......... .......... 31%  191M 1s
    ##  25450K .......... .......... .......... .......... .......... 31% 14.1M 1s
    ##  25500K .......... .......... .......... .......... .......... 31% 9.89M 1s
    ##  25550K .......... .......... .......... .......... .......... 32% 59.0M 1s
    ##  25600K .......... .......... .......... .......... .......... 32%  137M 1s
    ##  25650K .......... .......... .......... .......... .......... 32%  120M 1s
    ##  25700K .......... .......... .......... .......... .......... 32%  143M 1s
    ##  25750K .......... .......... .......... .......... .......... 32% 41.4M 1s
    ##  25800K .......... .......... .......... .......... .......... 32% 17.7M 1s
    ##  25850K .......... .......... .......... .......... .......... 32%  150M 1s
    ##  25900K .......... .......... .......... .......... .......... 32%  188M 1s
    ##  25950K .......... .......... .......... .......... .......... 32%  121M 1s
    ##  26000K .......... .......... .......... .......... .......... 32% 75.6M 1s
    ##  26050K .......... .......... .......... .......... .......... 32% 15.4M 1s
    ##  26100K .......... .......... .......... .......... .......... 32%  126M 1s
    ##  26150K .......... .......... .......... .......... .......... 32% 34.9M 1s
    ##  26200K .......... .......... .......... .......... .......... 32%  142M 1s
    ##  26250K .......... .......... .......... .......... .......... 32%  107M 1s
    ##  26300K .......... .......... .......... .......... .......... 32% 45.4M 1s
    ##  26350K .......... .......... .......... .......... .......... 33% 11.6M 1s
    ##  26400K .......... .......... .......... .......... .......... 33%  175M 1s
    ##  26450K .......... .......... .......... .......... .......... 33% 43.8M 1s
    ##  26500K .......... .......... .......... .......... .......... 33% 81.5M 1s
    ##  26550K .......... .......... .......... .......... .......... 33%  171M 1s
    ##  26600K .......... .......... .......... .......... .......... 33%  196M 1s
    ##  26650K .......... .......... .......... .......... .......... 33% 37.7M 1s
    ##  26700K .......... .......... .......... .......... .......... 33% 45.0M 1s
    ##  26750K .......... .......... .......... .......... .......... 33% 24.2M 1s
    ##  26800K .......... .......... .......... .......... .......... 33%  111M 1s
    ##  26850K .......... .......... .......... .......... .......... 33%  163M 1s
    ##  26900K .......... .......... .......... .......... .......... 33%  114M 1s
    ##  26950K .......... .......... .......... .......... .......... 33% 85.9M 1s
    ##  27000K .......... .......... .......... .......... .......... 33% 21.8M 1s
    ##  27050K .......... .......... .......... .......... .......... 33% 73.2M 1s
    ##  27100K .......... .......... .......... .......... .......... 33% 61.5M 1s
    ##  27150K .......... .......... .......... .......... .......... 34% 84.8M 1s
    ##  27200K .......... .......... .......... .......... .......... 34%  131M 1s
    ##  27250K .......... .......... .......... .......... .......... 34% 42.5M 1s
    ##  27300K .......... .......... .......... .......... .......... 34% 32.2M 1s
    ##  27350K .......... .......... .......... .......... .......... 34%  130M 1s
    ##  27400K .......... .......... .......... .......... .......... 34% 45.8M 1s
    ##  27450K .......... .......... .......... .......... .......... 34%  148M 1s
    ##  27500K .......... .......... .......... .......... .......... 34% 33.0M 1s
    ##  27550K .......... .......... .......... .......... .......... 34%  193M 1s
    ##  27600K .......... .......... .......... .......... .......... 34%  117M 1s
    ##  27650K .......... .......... .......... .......... .......... 34% 25.9M 1s
    ##  27700K .......... .......... .......... .......... .......... 34% 68.7M 1s
    ##  27750K .......... .......... .......... .......... .......... 34% 50.5M 1s
    ##  27800K .......... .......... .......... .......... .......... 34%  134M 1s
    ##  27850K .......... .......... .......... .......... .......... 34% 83.3M 1s
    ##  27900K .......... .......... .......... .......... .......... 34% 60.8M 1s
    ##  27950K .......... .......... .......... .......... .......... 35% 42.7M 1s
    ##  28000K .......... .......... .......... .......... .......... 35% 64.0M 1s
    ##  28050K .......... .......... .......... .......... .......... 35% 63.3M 1s
    ##  28100K .......... .......... .......... .......... .......... 35%  138M 1s
    ##  28150K .......... .......... .......... .......... .......... 35% 51.1M 1s
    ##  28200K .......... .......... .......... .......... .......... 35% 33.5M 1s
    ##  28250K .......... .......... .......... .......... .......... 35% 59.0M 1s
    ##  28300K .......... .......... .......... .......... .......... 35%  154M 1s
    ##  28350K .......... .......... .......... .......... .......... 35% 69.6M 1s
    ##  28400K .......... .......... .......... .......... .......... 35% 84.5M 1s
    ##  28450K .......... .......... .......... .......... .......... 35% 17.4M 1s
    ##  28500K .......... .......... .......... .......... .......... 35% 35.0M 1s
    ##  28550K .......... .......... .......... .......... .......... 35% 39.7M 1s
    ##  28600K .......... .......... .......... .......... .......... 35%  150M 1s
    ##  28650K .......... .......... .......... .......... .......... 35%  133M 1s
    ##  28700K .......... .......... .......... .......... .......... 35% 82.0M 1s
    ##  28750K .......... .......... .......... .......... .......... 36% 38.4M 1s
    ##  28800K .......... .......... .......... .......... .......... 36% 34.8M 1s
    ##  28850K .......... .......... .......... .......... .......... 36% 44.0M 1s
    ##  28900K .......... .......... .......... .......... .......... 36% 48.9M 1s
    ##  28950K .......... .......... .......... .......... .......... 36%  149M 1s
    ##  29000K .......... .......... .......... .......... .......... 36%  199M 1s
    ##  29050K .......... .......... .......... .......... .......... 36% 41.1M 1s
    ##  29100K .......... .......... .......... .......... .......... 36% 17.7M 1s
    ##  29150K .......... .......... .......... .......... .......... 36%  108M 1s
    ##  29200K .......... .......... .......... .......... .......... 36%  178M 1s
    ##  29250K .......... .......... .......... .......... .......... 36%  174M 1s
    ##  29300K .......... .......... .......... .......... .......... 36% 20.7M 1s
    ##  29350K .......... .......... .......... .......... .......... 36% 14.5M 1s
    ##  29400K .......... .......... .......... .......... .......... 36%  155M 1s
    ##  29450K .......... .......... .......... .......... .......... 36%  151M 1s
    ##  29500K .......... .......... .......... .......... .......... 36%  190M 1s
    ##  29550K .......... .......... .......... .......... .......... 37%  175M 1s
    ##  29600K .......... .......... .......... .......... .......... 37% 15.1M 1s
    ##  29650K .......... .......... .......... .......... .......... 37%  147M 1s
    ##  29700K .......... .......... .......... .......... .......... 37% 45.0M 1s
    ##  29750K .......... .......... .......... .......... .......... 37%  134M 1s
    ##  29800K .......... .......... .......... .......... .......... 37%  162M 1s
    ##  29850K .......... .......... .......... .......... .......... 37% 37.3M 1s
    ##  29900K .......... .......... .......... .......... .......... 37% 41.6M 1s
    ##  29950K .......... .......... .......... .......... .......... 37% 79.5M 1s
    ##  30000K .......... .......... .......... .......... .......... 37% 57.8M 1s
    ##  30050K .......... .......... .......... .......... .......... 37%  132M 1s
    ##  30100K .......... .......... .......... .......... .......... 37% 63.7M 1s
    ##  30150K .......... .......... .......... .......... .......... 37% 36.7M 1s
    ##  30200K .......... .......... .......... .......... .......... 37%  161M 1s
    ##  30250K .......... .......... .......... .......... .......... 37% 16.1M 1s
    ##  30300K .......... .......... .......... .......... .......... 37% 61.8M 1s
    ##  30350K .......... .......... .......... .......... .......... 38%  152M 1s
    ##  30400K .......... .......... .......... .......... .......... 38%  176M 1s
    ##  30450K .......... .......... .......... .......... .......... 38%  173M 1s
    ##  30500K .......... .......... .......... .......... .......... 38% 15.6M 1s
    ##  30550K .......... .......... .......... .......... .......... 38% 42.2M 1s
    ##  30600K .......... .......... .......... .......... .......... 38% 63.7M 1s
    ##  30650K .......... .......... .......... .......... .......... 38%  130M 1s
    ##  30700K .......... .......... .......... .......... .......... 38%  133M 1s
    ##  30750K .......... .......... .......... .......... .......... 38% 79.3M 1s
    ##  30800K .......... .......... .......... .......... .......... 38% 31.7M 1s
    ##  30850K .......... .......... .......... .......... .......... 38%  181M 1s
    ##  30900K .......... .......... .......... .......... .......... 38% 54.0M 1s
    ##  30950K .......... .......... .......... .......... .......... 38% 17.9M 1s
    ##  31000K .......... .......... .......... .......... .......... 38%  147M 1s
    ##  31050K .......... .......... .......... .......... .......... 38% 73.9M 1s
    ##  31100K .......... .......... .......... .......... .......... 38%  143M 1s
    ##  31150K .......... .......... .......... .......... .......... 39%  172M 1s
    ##  31200K .......... .......... .......... .......... .......... 39% 25.7M 1s
    ##  31250K .......... .......... .......... .......... .......... 39% 38.4M 1s
    ##  31300K .......... .......... .......... .......... .......... 39% 55.7M 1s
    ##  31350K .......... .......... .......... .......... .......... 39% 83.7M 1s
    ##  31400K .......... .......... .......... .......... .......... 39%  127M 1s
    ##  31450K .......... .......... .......... .......... .......... 39% 10.0M 1s
    ##  31500K .......... .......... .......... .......... .......... 39% 31.0M 1s
    ##  31550K .......... .......... .......... .......... .......... 39% 60.8M 1s
    ##  31600K .......... .......... .......... .......... .......... 39% 93.6M 1s
    ##  31650K .......... .......... .......... .......... .......... 39%  106M 1s
    ##  31700K .......... .......... .......... .......... .......... 39%  161M 1s
    ##  31750K .......... .......... .......... .......... .......... 39% 6.19M 1s
    ##  31800K .......... .......... .......... .......... .......... 39%  117M 1s
    ##  31850K .......... .......... .......... .......... .......... 39% 48.4M 1s
    ##  31900K .......... .......... .......... .......... .......... 39%  125M 1s
    ##  31950K .......... .......... .......... .......... .......... 40%  140M 1s
    ##  32000K .......... .......... .......... .......... .......... 40% 47.4M 1s
    ##  32050K .......... .......... .......... .......... .......... 40% 76.9M 1s
    ##  32100K .......... .......... .......... .......... .......... 40% 13.1M 1s
    ##  32150K .......... .......... .......... .......... .......... 40%  152M 1s
    ##  32200K .......... .......... .......... .......... .......... 40%  170M 1s
    ##  32250K .......... .......... .......... .......... .......... 40%  162M 1s
    ##  32300K .......... .......... .......... .......... .......... 40%  212M 1s
    ##  32350K .......... .......... .......... .......... .......... 40% 18.6M 1s
    ##  32400K .......... .......... .......... .......... .......... 40%  110M 1s
    ##  32450K .......... .......... .......... .......... .......... 40% 48.1M 1s
    ##  32500K .......... .......... .......... .......... .......... 40%  129M 1s
    ##  32550K .......... .......... .......... .......... .......... 40%  160M 1s
    ##  32600K .......... .......... .......... .......... .......... 40% 27.8M 1s
    ##  32650K .......... .......... .......... .......... .......... 40%  135M 1s
    ##  32700K .......... .......... .......... .......... .......... 40% 17.6M 1s
    ##  32750K .......... .......... .......... .......... .......... 41%  131M 1s
    ##  32800K .......... .......... .......... .......... .......... 41%  174M 1s
    ##  32850K .......... .......... .......... .......... .......... 41%  151M 1s
    ##  32900K .......... .......... .......... .......... .......... 41%  166M 1s
    ##  32950K .......... .......... .......... .......... .......... 41% 54.0M 1s
    ##  33000K .......... .......... .......... .......... .......... 41% 27.1M 1s
    ##  33050K .......... .......... .......... .......... .......... 41% 54.6M 1s
    ##  33100K .......... .......... .......... .......... .......... 41% 82.1M 1s
    ##  33150K .......... .......... .......... .......... .......... 41% 81.1M 1s
    ##  33200K .......... .......... .......... .......... .......... 41% 26.0M 1s
    ##  33250K .......... .......... .......... .......... .......... 41%  174M 1s
    ##  33300K .......... .......... .......... .......... .......... 41% 79.0M 1s
    ##  33350K .......... .......... .......... .......... .......... 41% 29.6M 1s
    ##  33400K .......... .......... .......... .......... .......... 41%  164M 1s
    ##  33450K .......... .......... .......... .......... .......... 41%  173M 1s
    ##  33500K .......... .......... .......... .......... .......... 41% 29.0M 1s
    ##  33550K .......... .......... .......... .......... .......... 42% 55.2M 1s
    ##  33600K .......... .......... .......... .......... .......... 42%  144M 1s
    ##  33650K .......... .......... .......... .......... .......... 42% 46.8M 1s
    ##  33700K .......... .......... .......... .......... .......... 42%  103M 1s
    ##  33750K .......... .......... .......... .......... .......... 42%  135M 1s
    ##  33800K .......... .......... .......... .......... .......... 42% 22.0M 1s
    ##  33850K .......... .......... .......... .......... .......... 42% 38.9M 1s
    ##  33900K .......... .......... .......... .......... .......... 42% 86.4M 1s
    ##  33950K .......... .......... .......... .......... .......... 42% 54.3M 1s
    ##  34000K .......... .......... .......... .......... .......... 42%  117M 1s
    ##  34050K .......... .......... .......... .......... .......... 42% 66.7M 1s
    ##  34100K .......... .......... .......... .......... .......... 42% 55.4M 1s
    ##  34150K .......... .......... .......... .......... .......... 42% 24.0M 1s
    ##  34200K .......... .......... .......... .......... .......... 42%  174M 1s
    ##  34250K .......... .......... .......... .......... .......... 42%  111M 1s
    ##  34300K .......... .......... .......... .......... .......... 42%  158M 1s
    ##  34350K .......... .......... .......... .......... .......... 43% 49.8M 1s
    ##  34400K .......... .......... .......... .......... .......... 43% 15.4M 1s
    ##  34450K .......... .......... .......... .......... .......... 43%  135M 1s
    ##  34500K .......... .......... .......... .......... .......... 43%  163M 1s
    ##  34550K .......... .......... .......... .......... .......... 43% 55.1M 1s
    ##  34600K .......... .......... .......... .......... .......... 43%  179M 1s
    ##  34650K .......... .......... .......... .......... .......... 43% 25.5M 1s
    ##  34700K .......... .......... .......... .......... .......... 43%  121M 1s
    ##  34750K .......... .......... .......... .......... .......... 43% 20.5M 1s
    ##  34800K .......... .......... .......... .......... .......... 43%  138M 1s
    ##  34850K .......... .......... .......... .......... .......... 43%  137M 1s
    ##  34900K .......... .......... .......... .......... .......... 43%  151M 1s
    ##  34950K .......... .......... .......... .......... .......... 43% 60.0M 1s
    ##  35000K .......... .......... .......... .......... .......... 43% 28.8M 1s
    ##  35050K .......... .......... .......... .......... .......... 43%  143M 1s
    ##  35100K .......... .......... .......... .......... .......... 43% 67.4M 1s
    ##  35150K .......... .......... .......... .......... .......... 44%  131M 1s
    ##  35200K .......... .......... .......... .......... .......... 44% 55.0M 1s
    ##  35250K .......... .......... .......... .......... .......... 44% 13.5M 1s
    ##  35300K .......... .......... .......... .......... .......... 44%  180M 1s
    ##  35350K .......... .......... .......... .......... .......... 44%  174M 1s
    ##  35400K .......... .......... .......... .......... .......... 44%  162M 1s
    ##  35450K .......... .......... .......... .......... .......... 44%  186M 1s
    ##  35500K .......... .......... .......... .......... .......... 44% 13.9M 1s
    ##  35550K .......... .......... .......... .......... .......... 44% 65.5M 1s
    ##  35600K .......... .......... .......... .......... .......... 44% 63.3M 1s
    ##  35650K .......... .......... .......... .......... .......... 44%  139M 1s
    ##  35700K .......... .......... .......... .......... .......... 44%  130M 1s
    ##  35750K .......... .......... .......... .......... .......... 44%  175M 1s
    ##  35800K .......... .......... .......... .......... .......... 44% 28.6M 1s
    ##  35850K .......... .......... .......... .......... .......... 44% 28.8M 1s
    ##  35900K .......... .......... .......... .......... .......... 44%  115M 1s
    ##  35950K .......... .......... .......... .......... .......... 45%  170M 1s
    ##  36000K .......... .......... .......... .......... .......... 45% 25.1M 1s
    ##  36050K .......... .......... .......... .......... .......... 45%  194M 1s
    ##  36100K .......... .......... .......... .......... .......... 45%  195M 1s
    ##  36150K .......... .......... .......... .......... .......... 45% 44.2M 1s
    ##  36200K .......... .......... .......... .......... .......... 45% 70.7M 1s
    ##  36250K .......... .......... .......... .......... .......... 45% 55.9M 1s
    ##  36300K .......... .......... .......... .......... .......... 45% 4.89M 1s
    ##  36350K .......... .......... .......... .......... .......... 45%  121M 1s
    ##  36400K .......... .......... .......... .......... .......... 45% 46.7M 1s
    ##  36450K .......... .......... .......... .......... .......... 45% 86.3M 1s
    ##  36500K .......... .......... .......... .......... .......... 45%  137M 1s
    ##  36550K .......... .......... .......... .......... .......... 45%  160M 1s
    ##  36600K .......... .......... .......... .......... .......... 45% 53.4M 1s
    ##  36650K .......... .......... .......... .......... .......... 45% 40.4M 1s
    ##  36700K .......... .......... .......... .......... .......... 45% 54.6M 1s
    ##  36750K .......... .......... .......... .......... .......... 46% 81.9M 1s
    ##  36800K .......... .......... .......... .......... .......... 46%  139M 1s
    ##  36850K .......... .......... .......... .......... .......... 46% 62.2M 1s
    ##  36900K .......... .......... .......... .......... .......... 46% 37.8M 1s
    ##  36950K .......... .......... .......... .......... .......... 46% 17.1M 1s
    ##  37000K .......... .......... .......... .......... .......... 46%  187M 1s
    ##  37050K .......... .......... .......... .......... .......... 46% 53.9M 1s
    ##  37100K .......... .......... .......... .......... .......... 46% 62.1M 1s
    ##  37150K .......... .......... .......... .......... .......... 46%  122M 1s
    ##  37200K .......... .......... .......... .......... .......... 46%  194M 1s
    ##  37250K .......... .......... .......... .......... .......... 46% 19.4M 1s
    ##  37300K .......... .......... .......... .......... .......... 46%  198M 1s
    ##  37350K .......... .......... .......... .......... .......... 46% 41.4M 1s
    ##  37400K .......... .......... .......... .......... .......... 46% 64.3M 1s
    ##  37450K .......... .......... .......... .......... .......... 46%  106M 1s
    ##  37500K .......... .......... .......... .......... .......... 46%  149M 1s
    ##  37550K .......... .......... .......... .......... .......... 47%  168M 1s
    ##  37600K .......... .......... .......... .......... .......... 47% 32.1M 1s
    ##  37650K .......... .......... .......... .......... .......... 47% 47.6M 1s
    ##  37700K .......... .......... .......... .......... .......... 47% 47.1M 1s
    ##  37750K .......... .......... .......... .......... .......... 47%  109M 1s
    ##  37800K .......... .......... .......... .......... .......... 47% 81.8M 1s
    ##  37850K .......... .......... .......... .......... .......... 47% 38.3M 1s
    ##  37900K .......... .......... .......... .......... .......... 47% 72.0M 1s
    ##  37950K .......... .......... .......... .......... .......... 47% 55.5M 1s
    ##  38000K .......... .......... .......... .......... .......... 47%  133M 1s
    ##  38050K .......... .......... .......... .......... .......... 47%  104M 1s
    ##  38100K .......... .......... .......... .......... .......... 47% 42.8M 1s
    ##  38150K .......... .......... .......... .......... .......... 47% 15.4M 1s
    ##  38200K .......... .......... .......... .......... .......... 47%  126M 1s
    ##  38250K .......... .......... .......... .......... .......... 47%  149M 1s
    ##  38300K .......... .......... .......... .......... .......... 47%  176M 1s
    ##  38350K .......... .......... .......... .......... .......... 48%  128M 1s
    ##  38400K .......... .......... .......... .......... .......... 48%  154M 1s
    ##  38450K .......... .......... .......... .......... .......... 48% 24.3M 1s
    ##  38500K .......... .......... .......... .......... .......... 48% 17.2M 1s
    ##  38550K .......... .......... .......... .......... .......... 48%  171M 1s
    ##  38600K .......... .......... .......... .......... .......... 48%  205M 1s
    ##  38650K .......... .......... .......... .......... .......... 48%  194M 1s
    ##  38700K .......... .......... .......... .......... .......... 48%  178M 1s
    ##  38750K .......... .......... .......... .......... .......... 48% 13.7M 1s
    ##  38800K .......... .......... .......... .......... .......... 48% 79.7M 1s
    ##  38850K .......... .......... .......... .......... .......... 48% 44.5M 1s
    ##  38900K .......... .......... .......... .......... .......... 48%  110M 1s
    ##  38950K .......... .......... .......... .......... .......... 48%  143M 1s
    ##  39000K .......... .......... .......... .......... .......... 48% 71.8M 1s
    ##  39050K .......... .......... .......... .......... .......... 48% 30.1M 1s
    ##  39100K .......... .......... .......... .......... .......... 48% 74.1M 1s
    ##  39150K .......... .......... .......... .......... .......... 49% 82.1M 1s
    ##  39200K .......... .......... .......... .......... .......... 49% 67.8M 1s
    ##  39250K .......... .......... .......... .......... .......... 49% 40.1M 1s
    ##  39300K .......... .......... .......... .......... .......... 49%  162M 1s
    ##  39350K .......... .......... .......... .......... .......... 49% 22.7M 1s
    ##  39400K .......... .......... .......... .......... .......... 49%  163M 1s
    ##  39450K .......... .......... .......... .......... .......... 49% 29.8M 1s
    ##  39500K .......... .......... .......... .......... .......... 49%  210M 1s
    ##  39550K .......... .......... .......... .......... .......... 49%  163M 1s
    ##  39600K .......... .......... .......... .......... .......... 49%  165M 1s
    ##  39650K .......... .......... .......... .......... .......... 49% 28.8M 1s
    ##  39700K .......... .......... .......... .......... .......... 49% 32.1M 1s
    ##  39750K .......... .......... .......... .......... .......... 49%  155M 1s
    ##  39800K .......... .......... .......... .......... .......... 49% 42.9M 1s
    ##  39850K .......... .......... .......... .......... .......... 49% 89.4M 1s
    ##  39900K .......... .......... .......... .......... .......... 49%  136M 1s
    ##  39950K .......... .......... .......... .......... .......... 50% 47.9M 1s
    ##  40000K .......... .......... .......... .......... .......... 50% 84.8M 1s
    ##  40050K .......... .......... .......... .......... .......... 50% 78.2M 1s
    ##  40100K .......... .......... .......... .......... .......... 50% 83.2M 1s
    ##  40150K .......... .......... .......... .......... .......... 50% 45.1M 1s
    ##  40200K .......... .......... .......... .......... .......... 50% 71.1M 1s
    ##  40250K .......... .......... .......... .......... .......... 50% 54.7M 1s
    ##  40300K .......... .......... .......... .......... .......... 50% 64.4M 1s
    ##  40350K .......... .......... .......... .......... .......... 50%  112M 1s
    ##  40400K .......... .......... .......... .......... .......... 50% 65.8M 1s
    ##  40450K .......... .......... .......... .......... .......... 50% 61.1M 1s
    ##  40500K .......... .......... .......... .......... .......... 50%  137M 1s
    ##  40550K .......... .......... .......... .......... .......... 50% 17.0M 1s
    ##  40600K .......... .......... .......... .......... .......... 50%  140M 1s
    ##  40650K .......... .......... .......... .......... .......... 50% 60.9M 1s
    ##  40700K .......... .......... .......... .......... .......... 50%  128M 1s
    ##  40750K .......... .......... .......... .......... .......... 51% 24.3M 1s
    ##  40800K .......... .......... .......... .......... .......... 51% 41.4M 1s
    ##  40850K .......... .......... .......... .......... .......... 51% 69.2M 1s
    ##  40900K .......... .......... .......... .......... .......... 51% 38.9M 1s
    ##  40950K .......... .......... .......... .......... .......... 51%  180M 1s
    ##  41000K .......... .......... .......... .......... .......... 51%  209M 1s
    ##  41050K .......... .......... .......... .......... .......... 51% 78.0M 1s
    ##  41100K .......... .......... .......... .......... .......... 51% 52.2M 1s
    ##  41150K .......... .......... .......... .......... .......... 51% 35.6M 1s
    ##  41200K .......... .......... .......... .......... .......... 51% 93.1M 1s
    ##  41250K .......... .......... .......... .......... .......... 51% 39.8M 1s
    ##  41300K .......... .......... .......... .......... .......... 51%  162M 1s
    ##  41350K .......... .......... .......... .......... .......... 51% 71.1M 1s
    ##  41400K .......... .......... .......... .......... .......... 51% 40.2M 1s
    ##  41450K .......... .......... .......... .......... .......... 51%  131M 1s
    ##  41500K .......... .......... .......... .......... .......... 51%  116M 1s
    ##  41550K .......... .......... .......... .......... .......... 52% 25.6M 1s
    ##  41600K .......... .......... .......... .......... .......... 52%  216M 1s
    ##  41650K .......... .......... .......... .......... .......... 52% 39.4M 1s
    ##  41700K .......... .......... .......... .......... .......... 52%  154M 1s
    ##  41750K .......... .......... .......... .......... .......... 52%  188M 1s
    ##  41800K .......... .......... .......... .......... .......... 52% 95.3M 1s
    ##  41850K .......... .......... .......... .......... .......... 52% 43.4M 1s
    ##  41900K .......... .......... .......... .......... .......... 52% 46.8M 1s
    ##  41950K .......... .......... .......... .......... .......... 52% 44.5M 1s
    ##  42000K .......... .......... .......... .......... .......... 52%  133M 1s
    ##  42050K .......... .......... .......... .......... .......... 52% 75.8M 1s
    ##  42100K .......... .......... .......... .......... .......... 52% 67.1M 1s
    ##  42150K .......... .......... .......... .......... .......... 52% 62.4M 1s
    ##  42200K .......... .......... .......... .......... .......... 52%  166M 1s
    ##  42250K .......... .......... .......... .......... .......... 52% 33.4M 1s
    ##  42300K .......... .......... .......... .......... .......... 52%  144M 1s
    ##  42350K .......... .......... .......... .......... .......... 53% 64.6M 1s
    ##  42400K .......... .......... .......... .......... .......... 53% 55.7M 1s
    ##  42450K .......... .......... .......... .......... .......... 53% 34.2M 1s
    ##  42500K .......... .......... .......... .......... .......... 53% 68.7M 1s
    ##  42550K .......... .......... .......... .......... .......... 53% 92.9M 1s
    ##  42600K .......... .......... .......... .......... .......... 53% 64.1M 1s
    ##  42650K .......... .......... .......... .......... .......... 53% 72.2M 1s
    ##  42700K .......... .......... .......... .......... .......... 53% 74.0M 1s
    ##  42750K .......... .......... .......... .......... .......... 53% 98.0M 1s
    ##  42800K .......... .......... .......... .......... .......... 53% 99.4M 1s
    ##  42850K .......... .......... .......... .......... .......... 53% 6.45M 1s
    ##  42900K .......... .......... .......... .......... .......... 53% 27.5M 1s
    ##  42950K .......... .......... .......... .......... .......... 53%  170M 1s
    ##  43000K .......... .......... .......... .......... .......... 53% 62.9M 1s
    ##  43050K .......... .......... .......... .......... .......... 53%  116M 1s
    ##  43100K .......... .......... .......... .......... .......... 53% 92.9M 1s
    ##  43150K .......... .......... .......... .......... .......... 54%  138M 1s
    ##  43200K .......... .......... .......... .......... .......... 54% 27.1M 1s
    ##  43250K .......... .......... .......... .......... .......... 54% 90.9M 1s
    ##  43300K .......... .......... .......... .......... .......... 54% 40.3M 1s
    ##  43350K .......... .......... .......... .......... .......... 54% 96.5M 1s
    ##  43400K .......... .......... .......... .......... .......... 54%  153M 1s
    ##  43450K .......... .......... .......... .......... .......... 54% 35.3M 1s
    ##  43500K .......... .......... .......... .......... .......... 54% 53.3M 1s
    ##  43550K .......... .......... .......... .......... .......... 54% 50.1M 1s
    ##  43600K .......... .......... .......... .......... .......... 54%  163M 1s
    ##  43650K .......... .......... .......... .......... .......... 54% 74.6M 1s
    ##  43700K .......... .......... .......... .......... .......... 54% 79.0M 1s
    ##  43750K .......... .......... .......... .......... .......... 54% 40.9M 1s
    ##  43800K .......... .......... .......... .......... .......... 54% 17.9M 1s
    ##  43850K .......... .......... .......... .......... .......... 54%  160M 1s
    ##  43900K .......... .......... .......... .......... .......... 54%  191M 1s
    ##  43950K .......... .......... .......... .......... .......... 55%  183M 1s
    ##  44000K .......... .......... .......... .......... .......... 55% 7.81M 1s
    ##  44050K .......... .......... .......... .......... .......... 55%  148M 1s
    ##  44100K .......... .......... .......... .......... .......... 55%  191M 1s
    ##  44150K .......... .......... .......... .......... .......... 55%  171M 1s
    ##  44200K .......... .......... .......... .......... .......... 55%  197M 1s
    ##  44250K .......... .......... .......... .......... .......... 55%  185M 1s
    ##  44300K .......... .......... .......... .......... .......... 55% 18.6M 1s
    ##  44350K .......... .......... .......... .......... .......... 55% 74.9M 1s
    ##  44400K .......... .......... .......... .......... .......... 55% 66.1M 1s
    ##  44450K .......... .......... .......... .......... .......... 55% 74.3M 1s
    ##  44500K .......... .......... .......... .......... .......... 55%  138M 1s
    ##  44550K .......... .......... .......... .......... .......... 55% 13.9M 1s
    ##  44600K .......... .......... .......... .......... .......... 55%  144M 1s
    ##  44650K .......... .......... .......... .......... .......... 55%  164M 1s
    ##  44700K .......... .......... .......... .......... .......... 55%  148M 1s
    ##  44750K .......... .......... .......... .......... .......... 56%  179M 1s
    ##  44800K .......... .......... .......... .......... .......... 56%  216M 1s
    ##  44850K .......... .......... .......... .......... .......... 56% 14.7M 1s
    ##  44900K .......... .......... .......... .......... .......... 56%  155M 1s
    ##  44950K .......... .......... .......... .......... .......... 56% 43.1M 1s
    ##  45000K .......... .......... .......... .......... .......... 56%  163M 1s
    ##  45050K .......... .......... .......... .......... .......... 56%  136M 1s
    ##  45100K .......... .......... .......... .......... .......... 56%  204M 1s
    ##  45150K .......... .......... .......... .......... .......... 56% 17.9M 1s
    ##  45200K .......... .......... .......... .......... .......... 56% 37.6M 1s
    ##  45250K .......... .......... .......... .......... .......... 56%  133M 1s
    ##  45300K .......... .......... .......... .......... .......... 56% 44.1M 1s
    ##  45350K .......... .......... .......... .......... .......... 56%  140M 1s
    ##  45400K .......... .......... .......... .......... .......... 56%  153M 1s
    ##  45450K .......... .......... .......... .......... .......... 56% 34.2M 1s
    ##  45500K .......... .......... .......... .......... .......... 56%  131M 1s
    ##  45550K .......... .......... .......... .......... .......... 57%  110M 1s
    ##  45600K .......... .......... .......... .......... .......... 57% 31.4M 1s
    ##  45650K .......... .......... .......... .......... .......... 57% 73.9M 1s
    ##  45700K .......... .......... .......... .......... .......... 57% 35.2M 1s
    ##  45750K .......... .......... .......... .......... .......... 57%  148M 1s
    ##  45800K .......... .......... .......... .......... .......... 57% 23.4M 1s
    ##  45850K .......... .......... .......... .......... .......... 57%  188M 1s
    ##  45900K .......... .......... .......... .......... .......... 57%  183M 1s
    ##  45950K .......... .......... .......... .......... .......... 57% 93.6M 1s
    ##  46000K .......... .......... .......... .......... .......... 57%  159M 1s
    ##  46050K .......... .......... .......... .......... .......... 57%  188M 1s
    ##  46100K .......... .......... .......... .......... .......... 57% 11.4M 1s
    ##  46150K .......... .......... .......... .......... .......... 57%  188M 1s
    ##  46200K .......... .......... .......... .......... .......... 57%  178M 1s
    ##  46250K .......... .......... .......... .......... .......... 57%  114M 1s
    ##  46300K .......... .......... .......... .......... .......... 57%  134M 1s
    ##  46350K .......... .......... .......... .......... .......... 58%  170M 1s
    ##  46400K .......... .......... .......... .......... .......... 58% 31.3M 1s
    ##  46450K .......... .......... .......... .......... .......... 58% 48.8M 1s
    ##  46500K .......... .......... .......... .......... .......... 58% 75.8M 1s
    ##  46550K .......... .......... .......... .......... .......... 58% 47.5M 1s
    ##  46600K .......... .......... .......... .......... .......... 58% 52.2M 1s
    ##  46650K .......... .......... .......... .......... .......... 58%  146M 1s
    ##  46700K .......... .......... .......... .......... .......... 58% 72.2M 1s
    ##  46750K .......... .......... .......... .......... .......... 58%  197M 1s
    ##  46800K .......... .......... .......... .......... .......... 58% 32.2M 1s
    ##  46850K .......... .......... .......... .......... .......... 58%  104M 1s
    ##  46900K .......... .......... .......... .......... .......... 58% 31.1M 1s
    ##  46950K .......... .......... .......... .......... .......... 58%  166M 1s
    ##  47000K .......... .......... .......... .......... .......... 58%  156M 1s
    ##  47050K .......... .......... .......... .......... .......... 58%  154M 1s
    ##  47100K .......... .......... .......... .......... .......... 58%  113M 1s
    ##  47150K .......... .......... .......... .......... .......... 59% 88.8M 1s
    ##  47200K .......... .......... .......... .......... .......... 59% 20.3M 1s
    ##  47250K .......... .......... .......... .......... .......... 59% 29.1M 1s
    ##  47300K .......... .......... .......... .......... .......... 59% 78.4M 1s
    ##  47350K .......... .......... .......... .......... .......... 59%  132M 1s
    ##  47400K .......... .......... .......... .......... .......... 59%  197M 1s
    ##  47450K .......... .......... .......... .......... .......... 59% 42.7M 1s
    ##  47500K .......... .......... .......... .......... .......... 59%  202M 1s
    ##  47550K .......... .......... .......... .......... .......... 59% 32.1M 1s
    ##  47600K .......... .......... .......... .......... .......... 59%  105M 1s
    ##  47650K .......... .......... .......... .......... .......... 59%  195M 1s
    ##  47700K .......... .......... .......... .......... .......... 59% 32.9M 1s
    ##  47750K .......... .......... .......... .......... .......... 59% 49.1M 1s
    ##  47800K .......... .......... .......... .......... .......... 59% 74.9M 1s
    ##  47850K .......... .......... .......... .......... .......... 59%  146M 1s
    ##  47900K .......... .......... .......... .......... .......... 59% 30.0M 1s
    ##  47950K .......... .......... .......... .......... .......... 60% 40.8M 1s
    ##  48000K .......... .......... .......... .......... .......... 60%  148M 1s
    ##  48050K .......... .......... .......... .......... .......... 60% 87.3M 1s
    ##  48100K .......... .......... .......... .......... .......... 60%  141M 1s
    ##  48150K .......... .......... .......... .......... .......... 60%  177M 1s
    ##  48200K .......... .......... .......... .......... .......... 60% 22.0M 1s
    ##  48250K .......... .......... .......... .......... .......... 60%  173M 1s
    ##  48300K .......... .......... .......... .......... .......... 60%  162M 1s
    ##  48350K .......... .......... .......... .......... .......... 60%  169M 1s
    ##  48400K .......... .......... .......... .......... .......... 60%  123M 1s
    ##  48450K .......... .......... .......... .......... .......... 60% 28.8M 1s
    ##  48500K .......... .......... .......... .......... .......... 60% 32.1M 1s
    ##  48550K .......... .......... .......... .......... .......... 60%  143M 1s
    ##  48600K .......... .......... .......... .......... .......... 60%  167M 1s
    ##  48650K .......... .......... .......... .......... .......... 60%  157M 1s
    ##  48700K .......... .......... .......... .......... .......... 60%  196M 1s
    ##  48750K .......... .......... .......... .......... .......... 61% 19.3M 1s
    ##  48800K .......... .......... .......... .......... .......... 61%  167M 1s
    ##  48850K .......... .......... .......... .......... .......... 61%  111M 1s
    ##  48900K .......... .......... .......... .......... .......... 61% 34.5M 1s
    ##  48950K .......... .......... .......... .......... .......... 61%  139M 1s
    ##  49000K .......... .......... .......... .......... .......... 61%  154M 1s
    ##  49050K .......... .......... .......... .......... .......... 61%  176M 1s
    ##  49100K .......... .......... .......... .......... .......... 61% 61.4M 1s
    ##  49150K .......... .......... .......... .......... .......... 61% 21.0M 1s
    ##  49200K .......... .......... .......... .......... .......... 61% 86.6M 1s
    ##  49250K .......... .......... .......... .......... .......... 61% 54.1M 1s
    ##  49300K .......... .......... .......... .......... .......... 61%  108M 1s
    ##  49350K .......... .......... .......... .......... .......... 61%  140M 1s
    ##  49400K .......... .......... .......... .......... .......... 61% 72.6M 1s
    ##  49450K .......... .......... .......... .......... .......... 61% 40.5M 1s
    ##  49500K .......... .......... .......... .......... .......... 61%  173M 1s
    ##  49550K .......... .......... .......... .......... .......... 62% 10.1M 1s
    ##  49600K .......... .......... .......... .......... .......... 62% 62.2M 1s
    ##  49650K .......... .......... .......... .......... .......... 62% 76.7M 1s
    ##  49700K .......... .......... .......... .......... .......... 62%  102M 1s
    ##  49750K .......... .......... .......... .......... .......... 62%  138M 1s
    ##  49800K .......... .......... .......... .......... .......... 62%  157M 1s
    ##  49850K .......... .......... .......... .......... .......... 62% 40.6M 1s
    ##  49900K .......... .......... .......... .......... .......... 62% 28.5M 1s
    ##  49950K .......... .......... .......... .......... .......... 62%  123M 1s
    ##  50000K .......... .......... .......... .......... .......... 62% 49.8M 1s
    ##  50050K .......... .......... .......... .......... .......... 62%  122M 1s
    ##  50100K .......... .......... .......... .......... .......... 62%  124M 1s
    ##  50150K .......... .......... .......... .......... .......... 62% 31.2M 1s
    ##  50200K .......... .......... .......... .......... .......... 62%  197M 1s
    ##  50250K .......... .......... .......... .......... .......... 62%  199M 1s
    ##  50300K .......... .......... .......... .......... .......... 62% 49.5M 1s
    ##  50350K .......... .......... .......... .......... .......... 63% 37.9M 1s
    ##  50400K .......... .......... .......... .......... .......... 63%  195M 1s
    ##  50450K .......... .......... .......... .......... .......... 63% 45.3M 1s
    ##  50500K .......... .......... .......... .......... .......... 63% 30.2M 1s
    ##  50550K .......... .......... .......... .......... .......... 63%  137M 1s
    ##  50600K .......... .......... .......... .......... .......... 63%  180M 1s
    ##  50650K .......... .......... .......... .......... .......... 63% 22.8M 1s
    ##  50700K .......... .......... .......... .......... .......... 63%  183M 1s
    ##  50750K .......... .......... .......... .......... .......... 63%  160M 1s
    ##  50800K .......... .......... .......... .......... .......... 63%  217M 1s
    ##  50850K .......... .......... .......... .......... .......... 63%  198M 1s
    ##  50900K .......... .......... .......... .......... .......... 63% 31.3M 1s
    ##  50950K .......... .......... .......... .......... .......... 63% 82.5M 1s
    ##  51000K .......... .......... .......... .......... .......... 63% 29.3M 1s
    ##  51050K .......... .......... .......... .......... .......... 63% 66.4M 1s
    ##  51100K .......... .......... .......... .......... .......... 63%  150M 1s
    ##  51150K .......... .......... .......... .......... .......... 64%  148M 1s
    ##  51200K .......... .......... .......... .......... .......... 64% 28.1M 1s
    ##  51250K .......... .......... .......... .......... .......... 64%  205M 1s
    ##  51300K .......... .......... .......... .......... .......... 64% 18.7M 1s
    ##  51350K .......... .......... .......... .......... .......... 64% 57.3M 1s
    ##  51400K .......... .......... .......... .......... .......... 64%  172M 1s
    ##  51450K .......... .......... .......... .......... .......... 64%  169M 1s
    ##  51500K .......... .......... .......... .......... .......... 64%  209M 1s
    ##  51550K .......... .......... .......... .......... .......... 64% 29.9M 1s
    ##  51600K .......... .......... .......... .......... .......... 64% 27.6M 1s
    ##  51650K .......... .......... .......... .......... .......... 64%  156M 1s
    ##  51700K .......... .......... .......... .......... .......... 64% 29.1M 1s
    ##  51750K .......... .......... .......... .......... .......... 64% 43.5M 1s
    ##  51800K .......... .......... .......... .......... .......... 64%  163M 1s
    ##  51850K .......... .......... .......... .......... .......... 64%  174M 1s
    ##  51900K .......... .......... .......... .......... .......... 65%  245M 1s
    ##  51950K .......... .......... .......... .......... .......... 65%  172M 1s
    ##  52000K .......... .......... .......... .......... .......... 65% 50.7M 1s
    ##  52050K .......... .......... .......... .......... .......... 65%  101M 1s
    ##  52100K .......... .......... .......... .......... .......... 65% 30.4M 1s
    ##  52150K .......... .......... .......... .......... .......... 65%  115M 1s
    ##  52200K .......... .......... .......... .......... .......... 65%  231M 1s
    ##  52250K .......... .......... .......... .......... .......... 65% 33.9M 1s
    ##  52300K .......... .......... .......... .......... .......... 65%  184M 1s
    ##  52350K .......... .......... .......... .......... .......... 65% 33.9M 1s
    ##  52400K .......... .......... .......... .......... .......... 65%  179M 1s
    ##  52450K .......... .......... .......... .......... .......... 65%  179M 1s
    ##  52500K .......... .......... .......... .......... .......... 65%  187M 1s
    ##  52550K .......... .......... .......... .......... .......... 65% 34.0M 1s
    ##  52600K .......... .......... .......... .......... .......... 65%  142M 1s
    ##  52650K .......... .......... .......... .......... .......... 65%  181M 1s
    ##  52700K .......... .......... .......... .......... .......... 66% 57.7M 1s
    ##  52750K .......... .......... .......... .......... .......... 66% 94.8M 1s
    ##  52800K .......... .......... .......... .......... .......... 66%  205M 1s
    ##  52850K .......... .......... .......... .......... .......... 66%  105M 1s
    ##  52900K .......... .......... .......... .......... .......... 66% 34.4M 1s
    ##  52950K .......... .......... .......... .......... .......... 66%  104M 1s
    ##  53000K .......... .......... .......... .......... .......... 66% 63.1M 1s
    ##  53050K .......... .......... .......... .......... .......... 66% 24.8M 1s
    ##  53100K .......... .......... .......... .......... .......... 66%  161M 1s
    ##  53150K .......... .......... .......... .......... .......... 66%  180M 1s
    ##  53200K .......... .......... .......... .......... .......... 66% 58.5M 1s
    ##  53250K .......... .......... .......... .......... .......... 66%  129M 1s
    ##  53300K .......... .......... .......... .......... .......... 66%  244M 1s
    ##  53350K .......... .......... .......... .......... .......... 66%  181M 1s
    ##  53400K .......... .......... .......... .......... .......... 66% 39.0M 1s
    ##  53450K .......... .......... .......... .......... .......... 66% 27.7M 1s
    ##  53500K .......... .......... .......... .......... .......... 67%  133M 1s
    ##  53550K .......... .......... .......... .......... .......... 67%  188M 1s
    ##  53600K .......... .......... .......... .......... .......... 67%  191M 1s
    ##  53650K .......... .......... .......... .......... .......... 67%  148M 1s
    ##  53700K .......... .......... .......... .......... .......... 67% 60.8M 1s
    ##  53750K .......... .......... .......... .......... .......... 67% 25.5M 1s
    ##  53800K .......... .......... .......... .......... .......... 67% 39.1M 1s
    ##  53850K .......... .......... .......... .......... .......... 67%  159M 1s
    ##  53900K .......... .......... .......... .......... .......... 67% 32.8M 1s
    ##  53950K .......... .......... .......... .......... .......... 67% 96.3M 1s
    ##  54000K .......... .......... .......... .......... .......... 67%  124M 1s
    ##  54050K .......... .......... .......... .......... .......... 67%  136M 1s
    ##  54100K .......... .......... .......... .......... .......... 67% 79.4M 1s
    ##  54150K .......... .......... .......... .......... .......... 67% 31.2M 1s
    ##  54200K .......... .......... .......... .......... .......... 67% 86.9M 1s
    ##  54250K .......... .......... .......... .......... .......... 67% 82.5M 1s
    ##  54300K .......... .......... .......... .......... .......... 68% 50.1M 1s
    ##  54350K .......... .......... .......... .......... .......... 68% 65.0M 1s
    ##  54400K .......... .......... .......... .......... .......... 68%  163M 1s
    ##  54450K .......... .......... .......... .......... .......... 68%  130M 1s
    ##  54500K .......... .......... .......... .......... .......... 68%  112M 1s
    ##  54550K .......... .......... .......... .......... .......... 68% 13.4M 1s
    ##  54600K .......... .......... .......... .......... .......... 68%  135M 1s
    ##  54650K .......... .......... .......... .......... .......... 68%  124M 1s
    ##  54700K .......... .......... .......... .......... .......... 68%  181M 1s
    ##  54750K .......... .......... .......... .......... .......... 68%  176M 1s
    ##  54800K .......... .......... .......... .......... .......... 68%  218M 1s
    ##  54850K .......... .......... .......... .......... .......... 68% 19.6M 1s
    ##  54900K .......... .......... .......... .......... .......... 68% 81.1M 1s
    ##  54950K .......... .......... .......... .......... .......... 68% 32.6M 1s
    ##  55000K .......... .......... .......... .......... .......... 68%  174M 1s
    ##  55050K .......... .......... .......... .......... .......... 68%  173M 1s
    ##  55100K .......... .......... .......... .......... .......... 69%  160M 1s
    ##  55150K .......... .......... .......... .......... .......... 69%  143M 1s
    ##  55200K .......... .......... .......... .......... .......... 69%  173M 1s
    ##  55250K .......... .......... .......... .......... .......... 69% 26.3M 1s
    ##  55300K .......... .......... .......... .......... .......... 69%  132M 1s
    ##  55350K .......... .......... .......... .......... .......... 69%  175M 1s
    ##  55400K .......... .......... .......... .......... .......... 69%  111M 1s
    ##  55450K .......... .......... .......... .......... .......... 69% 35.1M 1s
    ##  55500K .......... .......... .......... .......... .......... 69%  112M 1s
    ##  55550K .......... .......... .......... .......... .......... 69% 16.1M 1s
    ##  55600K .......... .......... .......... .......... .......... 69%  204M 1s
    ##  55650K .......... .......... .......... .......... .......... 69%  164M 1s
    ##  55700K .......... .......... .......... .......... .......... 69%  214M 1s
    ##  55750K .......... .......... .......... .......... .......... 69%  176M 1s
    ##  55800K .......... .......... .......... .......... .......... 69% 20.4M 1s
    ##  55850K .......... .......... .......... .......... .......... 69% 37.1M 1s
    ##  55900K .......... .......... .......... .......... .......... 70% 56.2M 1s
    ##  55950K .......... .......... .......... .......... .......... 70% 61.4M 1s
    ##  56000K .......... .......... .......... .......... .......... 70%  129M 1s
    ##  56050K .......... .......... .......... .......... .......... 70%  130M 1s
    ##  56100K .......... .......... .......... .......... .......... 70%  167M 1s
    ##  56150K .......... .......... .......... .......... .......... 70%  119M 1s
    ##  56200K .......... .......... .......... .......... .......... 70% 29.4M 1s
    ##  56250K .......... .......... .......... .......... .......... 70% 36.0M 1s
    ##  56300K .......... .......... .......... .......... .......... 70% 66.1M 1s
    ##  56350K .......... .......... .......... .......... .......... 70% 42.4M 1s
    ##  56400K .......... .......... .......... .......... .......... 70%  139M 1s
    ##  56450K .......... .......... .......... .......... .......... 70%  108M 1s
    ##  56500K .......... .......... .......... .......... .......... 70% 82.6M 1s
    ##  56550K .......... .......... .......... .......... .......... 70%  107M 1s
    ##  56600K .......... .......... .......... .......... .......... 70%  106M 0s
    ##  56650K .......... .......... .......... .......... .......... 70% 16.7M 0s
    ##  56700K .......... .......... .......... .......... .......... 71%  144M 0s
    ##  56750K .......... .......... .......... .......... .......... 71% 19.4M 0s
    ##  56800K .......... .......... .......... .......... .......... 71%  207M 0s
    ##  56850K .......... .......... .......... .......... .......... 71%  171M 0s
    ##  56900K .......... .......... .......... .......... .......... 71%  209M 0s
    ##  56950K .......... .......... .......... .......... .......... 71% 21.9M 0s
    ##  57000K .......... .......... .......... .......... .......... 71% 30.7M 0s
    ##  57050K .......... .......... .......... .......... .......... 71% 57.6M 0s
    ##  57100K .......... .......... .......... .......... .......... 71% 87.6M 0s
    ##  57150K .......... .......... .......... .......... .......... 71%  129M 0s
    ##  57200K .......... .......... .......... .......... .......... 71%  142M 0s
    ##  57250K .......... .......... .......... .......... .......... 71%  110M 0s
    ##  57300K .......... .......... .......... .......... .......... 71%  187M 0s
    ##  57350K .......... .......... .......... .......... .......... 71% 36.9M 0s
    ##  57400K .......... .......... .......... .......... .......... 71% 31.0M 0s
    ##  57450K .......... .......... .......... .......... .......... 71% 90.3M 0s
    ##  57500K .......... .......... .......... .......... .......... 72%  104M 0s
    ##  57550K .......... .......... .......... .......... .......... 72% 27.2M 0s
    ##  57600K .......... .......... .......... .......... .......... 72%  174M 0s
    ##  57650K .......... .......... .......... .......... .......... 72%  187M 0s
    ##  57700K .......... .......... .......... .......... .......... 72% 31.1M 0s
    ##  57750K .......... .......... .......... .......... .......... 72%  153M 0s
    ##  57800K .......... .......... .......... .......... .......... 72%  206M 0s
    ##  57850K .......... .......... .......... .......... .......... 72% 53.5M 0s
    ##  57900K .......... .......... .......... .......... .......... 72% 82.9M 0s
    ##  57950K .......... .......... .......... .......... .......... 72% 82.5M 0s
    ##  58000K .......... .......... .......... .......... .......... 72% 68.8M 0s
    ##  58050K .......... .......... .......... .......... .......... 72% 50.3M 0s
    ##  58100K .......... .......... .......... .......... .......... 72% 96.9M 0s
    ##  58150K .......... .......... .......... .......... .......... 72% 54.9M 0s
    ##  58200K .......... .......... .......... .......... .......... 72% 91.7M 0s
    ##  58250K .......... .......... .......... .......... .......... 72% 85.6M 0s
    ##  58300K .......... .......... .......... .......... .......... 73% 88.7M 0s
    ##  58350K .......... .......... .......... .......... .......... 73% 31.2M 0s
    ##  58400K .......... .......... .......... .......... .......... 73% 84.6M 0s
    ##  58450K .......... .......... .......... .......... .......... 73% 86.9M 0s
    ##  58500K .......... .......... .......... .......... .......... 73% 43.5M 0s
    ##  58550K .......... .......... .......... .......... .......... 73% 49.7M 0s
    ##  58600K .......... .......... .......... .......... .......... 73% 96.5M 0s
    ##  58650K .......... .......... .......... .......... .......... 73% 87.7M 0s
    ##  58700K .......... .......... .......... .......... .......... 73% 28.3M 0s
    ##  58750K .......... .......... .......... .......... .......... 73% 85.2M 0s
    ##  58800K .......... .......... .......... .......... .......... 73% 81.4M 0s
    ##  58850K .......... .......... .......... .......... .......... 73% 86.4M 0s
    ##  58900K .......... .......... .......... .......... .......... 73% 87.0M 0s
    ##  58950K .......... .......... .......... .......... .......... 73%  114M 0s
    ##  59000K .......... .......... .......... .......... .......... 73% 97.2M 0s
    ##  59050K .......... .......... .......... .......... .......... 73% 71.2M 0s
    ##  59100K .......... .......... .......... .......... .......... 74% 91.3M 0s
    ##  59150K .......... .......... .......... .......... .......... 74% 82.9M 0s
    ##  59200K .......... .......... .......... .......... .......... 74% 80.9M 0s
    ##  59250K .......... .......... .......... .......... .......... 74%  101M 0s
    ##  59300K .......... .......... .......... .......... .......... 74% 99.2M 0s
    ##  59350K .......... .......... .......... .......... .......... 74% 54.3M 0s
    ##  59400K .......... .......... .......... .......... .......... 74% 97.5M 0s
    ##  59450K .......... .......... .......... .......... .......... 74% 51.4M 0s
    ##  59500K .......... .......... .......... .......... .......... 74% 42.8M 0s
    ##  59550K .......... .......... .......... .......... .......... 74% 56.4M 0s
    ##  59600K .......... .......... .......... .......... .......... 74% 40.6M 0s
    ##  59650K .......... .......... .......... .......... .......... 74% 74.2M 0s
    ##  59700K .......... .......... .......... .......... .......... 74% 30.5M 0s
    ##  59750K .......... .......... .......... .......... .......... 74% 97.9M 0s
    ##  59800K .......... .......... .......... .......... .......... 74% 24.5M 0s
    ##  59850K .......... .......... .......... .......... .......... 74% 79.9M 0s
    ##  59900K .......... .......... .......... .......... .......... 75%  105M 0s
    ##  59950K .......... .......... .......... .......... .......... 75%  104M 0s
    ##  60000K .......... .......... .......... .......... .......... 75%  124M 0s
    ##  60050K .......... .......... .......... .......... .......... 75% 19.6M 0s
    ##  60100K .......... .......... .......... .......... .......... 75% 44.0M 0s
    ##  60150K .......... .......... .......... .......... .......... 75% 42.2M 0s
    ##  60200K .......... .......... .......... .......... .......... 75% 50.4M 0s
    ##  60250K .......... .......... .......... .......... .......... 75% 76.6M 0s
    ##  60300K .......... .......... .......... .......... .......... 75%  101M 0s
    ##  60350K .......... .......... .......... .......... .......... 75% 87.9M 0s
    ##  60400K .......... .......... .......... .......... .......... 75%  117M 0s
    ##  60450K .......... .......... .......... .......... .......... 75% 5.83M 0s
    ##  60500K .......... .......... .......... .......... .......... 75%  118M 0s
    ##  60550K .......... .......... .......... .......... .......... 75%  109M 0s
    ##  60600K .......... .......... .......... .......... .......... 75% 98.3M 0s
    ##  60650K .......... .......... .......... .......... .......... 75%  117M 0s
    ##  60700K .......... .......... .......... .......... .......... 76%  113M 0s
    ##  60750K .......... .......... .......... .......... .......... 76%  103M 0s
    ##  60800K .......... .......... .......... .......... .......... 76%  112M 0s
    ##  60850K .......... .......... .......... .......... .......... 76% 62.2M 0s
    ##  60900K .......... .......... .......... .......... .......... 76% 81.1M 0s
    ##  60950K .......... .......... .......... .......... .......... 76% 76.0M 0s
    ##  61000K .......... .......... .......... .......... .......... 76% 88.3M 0s
    ##  61050K .......... .......... .......... .......... .......... 76% 93.9M 0s
    ##  61100K .......... .......... .......... .......... .......... 76%  110M 0s
    ##  61150K .......... .......... .......... .......... .......... 76%  102M 0s
    ##  61200K .......... .......... .......... .......... .......... 76% 26.6M 0s
    ##  61250K .......... .......... .......... .......... .......... 76% 84.1M 0s
    ##  61300K .......... .......... .......... .......... .......... 76%  119M 0s
    ##  61350K .......... .......... .......... .......... .......... 76%  102M 0s
    ##  61400K .......... .......... .......... .......... .......... 76%  115M 0s
    ##  61450K .......... .......... .......... .......... .......... 76%  121M 0s
    ##  61500K .......... .......... .......... .......... .......... 77% 98.4M 0s
    ##  61550K .......... .......... .......... .......... .......... 77% 29.4M 0s
    ##  61600K .......... .......... .......... .......... .......... 77%  117M 0s
    ##  61650K .......... .......... .......... .......... .......... 77% 56.4M 0s
    ##  61700K .......... .......... .......... .......... .......... 77% 99.1M 0s
    ##  61750K .......... .......... .......... .......... .......... 77%  143M 0s
    ##  61800K .......... .......... .......... .......... .......... 77% 27.1M 0s
    ##  61850K .......... .......... .......... .......... .......... 77% 93.5M 0s
    ##  61900K .......... .......... .......... .......... .......... 77% 94.9M 0s
    ##  61950K .......... .......... .......... .......... .......... 77% 82.8M 0s
    ##  62000K .......... .......... .......... .......... .......... 77%  113M 0s
    ##  62050K .......... .......... .......... .......... .......... 77%  139M 0s
    ##  62100K .......... .......... .......... .......... .......... 77% 35.1M 0s
    ##  62150K .......... .......... .......... .......... .......... 77%  121M 0s
    ##  62200K .......... .......... .......... .......... .......... 77% 42.3M 0s
    ##  62250K .......... .......... .......... .......... .......... 77% 80.8M 0s
    ##  62300K .......... .......... .......... .......... .......... 78% 90.6M 0s
    ##  62350K .......... .......... .......... .......... .......... 78% 9.82M 0s
    ##  62400K .......... .......... .......... .......... .......... 78% 86.9M 0s
    ##  62450K .......... .......... .......... .......... .......... 78%  131M 0s
    ##  62500K .......... .......... .......... .......... .......... 78%  102M 0s
    ##  62550K .......... .......... .......... .......... .......... 78%  121M 0s
    ##  62600K .......... .......... .......... .......... .......... 78%  130M 0s
    ##  62650K .......... .......... .......... .......... .......... 78%  138M 0s
    ##  62700K .......... .......... .......... .......... .......... 78%  117M 0s
    ##  62750K .......... .......... .......... .......... .......... 78%  137M 0s
    ##  62800K .......... .......... .......... .......... .......... 78%  141M 0s
    ##  62850K .......... .......... .......... .......... .......... 78% 54.6M 0s
    ##  62900K .......... .......... .......... .......... .......... 78% 30.1M 0s
    ##  62950K .......... .......... .......... .......... .......... 78%  161M 0s
    ##  63000K .......... .......... .......... .......... .......... 78% 92.5M 0s
    ##  63050K .......... .......... .......... .......... .......... 78%  116M 0s
    ##  63100K .......... .......... .......... .......... .......... 79%  153M 0s
    ##  63150K .......... .......... .......... .......... .......... 79% 57.5M 0s
    ##  63200K .......... .......... .......... .......... .......... 79% 53.9M 0s
    ##  63250K .......... .......... .......... .......... .......... 79%  155M 0s
    ##  63300K .......... .......... .......... .......... .......... 79% 12.3M 0s
    ##  63350K .......... .......... .......... .......... .......... 79%  139M 0s
    ##  63400K .......... .......... .......... .......... .......... 79%  124M 0s
    ##  63450K .......... .......... .......... .......... .......... 79%  136M 0s
    ##  63500K .......... .......... .......... .......... .......... 79%  131M 0s
    ##  63550K .......... .......... .......... .......... .......... 79%  140M 0s
    ##  63600K .......... .......... .......... .......... .......... 79% 37.4M 0s
    ##  63650K .......... .......... .......... .......... .......... 79% 45.5M 0s
    ##  63700K .......... .......... .......... .......... .......... 79%  125M 0s
    ##  63750K .......... .......... .......... .......... .......... 79% 38.2M 0s
    ##  63800K .......... .......... .......... .......... .......... 79%  109M 0s
    ##  63850K .......... .......... .......... .......... .......... 79%  127M 0s
    ##  63900K .......... .......... .......... .......... .......... 80% 76.3M 0s
    ##  63950K .......... .......... .......... .......... .......... 80%  124M 0s
    ##  64000K .......... .......... .......... .......... .......... 80%  102M 0s
    ##  64050K .......... .......... .......... .......... .......... 80% 16.9M 0s
    ##  64100K .......... .......... .......... .......... .......... 80%  103M 0s
    ##  64150K .......... .......... .......... .......... .......... 80% 29.4M 0s
    ##  64200K .......... .......... .......... .......... .......... 80%  107M 0s
    ##  64250K .......... .......... .......... .......... .......... 80%  135M 0s
    ##  64300K .......... .......... .......... .......... .......... 80%  160M 0s
    ##  64350K .......... .......... .......... .......... .......... 80% 19.3M 0s
    ##  64400K .......... .......... .......... .......... .......... 80%  117M 0s
    ##  64450K .......... .......... .......... .......... .......... 80% 47.4M 0s
    ##  64500K .......... .......... .......... .......... .......... 80% 67.8M 0s
    ##  64550K .......... .......... .......... .......... .......... 80% 94.0M 0s
    ##  64600K .......... .......... .......... .......... .......... 80%  127M 0s
    ##  64650K .......... .......... .......... .......... .......... 80%  151M 0s
    ##  64700K .......... .......... .......... .......... .......... 81% 30.2M 0s
    ##  64750K .......... .......... .......... .......... .......... 81%  116M 0s
    ##  64800K .......... .......... .......... .......... .......... 81% 77.5M 0s
    ##  64850K .......... .......... .......... .......... .......... 81%  124M 0s
    ##  64900K .......... .......... .......... .......... .......... 81% 53.4M 0s
    ##  64950K .......... .......... .......... .......... .......... 81%  134M 0s
    ##  65000K .......... .......... .......... .......... .......... 81% 79.2M 0s
    ##  65050K .......... .......... .......... .......... .......... 81%  117M 0s
    ##  65100K .......... .......... .......... .......... .......... 81% 28.2M 0s
    ##  65150K .......... .......... .......... .......... .......... 81% 81.1M 0s
    ##  65200K .......... .......... .......... .......... .......... 81%  104M 0s
    ##  65250K .......... .......... .......... .......... .......... 81% 52.2M 0s
    ##  65300K .......... .......... .......... .......... .......... 81% 95.8M 0s
    ##  65350K .......... .......... .......... .......... .......... 81%  146M 0s
    ##  65400K .......... .......... .......... .......... .......... 81% 34.9M 0s
    ##  65450K .......... .......... .......... .......... .......... 81%  117M 0s
    ##  65500K .......... .......... .......... .......... .......... 82% 66.3M 0s
    ##  65550K .......... .......... .......... .......... .......... 82% 39.2M 0s
    ##  65600K .......... .......... .......... .......... .......... 82% 87.7M 0s
    ##  65650K .......... .......... .......... .......... .......... 82%  112M 0s
    ##  65700K .......... .......... .......... .......... .......... 82% 38.8M 0s
    ##  65750K .......... .......... .......... .......... .......... 82% 85.7M 0s
    ##  65800K .......... .......... .......... .......... .......... 82%  133M 0s
    ##  65850K .......... .......... .......... .......... .......... 82% 81.2M 0s
    ##  65900K .......... .......... .......... .......... .......... 82% 69.8M 0s
    ##  65950K .......... .......... .......... .......... .......... 82% 54.2M 0s
    ##  66000K .......... .......... .......... .......... .......... 82%  123M 0s
    ##  66050K .......... .......... .......... .......... .......... 82% 40.1M 0s
    ##  66100K .......... .......... .......... .......... .......... 82% 85.8M 0s
    ##  66150K .......... .......... .......... .......... .......... 82% 43.9M 0s
    ##  66200K .......... .......... .......... .......... .......... 82% 92.1M 0s
    ##  66250K .......... .......... .......... .......... .......... 82%  142M 0s
    ##  66300K .......... .......... .......... .......... .......... 83% 77.7M 0s
    ##  66350K .......... .......... .......... .......... .......... 83% 62.8M 0s
    ##  66400K .......... .......... .......... .......... .......... 83% 79.8M 0s
    ##  66450K .......... .......... .......... .......... .......... 83% 64.2M 0s
    ##  66500K .......... .......... .......... .......... .......... 83% 96.5M 0s
    ##  66550K .......... .......... .......... .......... .......... 83% 52.3M 0s
    ##  66600K .......... .......... .......... .......... .......... 83% 66.2M 0s
    ##  66650K .......... .......... .......... .......... .......... 83% 41.1M 0s
    ##  66700K .......... .......... .......... .......... .......... 83% 69.6M 0s
    ##  66750K .......... .......... .......... .......... .......... 83% 88.3M 0s
    ##  66800K .......... .......... .......... .......... .......... 83%  124M 0s
    ##  66850K .......... .......... .......... .......... .......... 83% 69.3M 0s
    ##  66900K .......... .......... .......... .......... .......... 83% 63.4M 0s
    ##  66950K .......... .......... .......... .......... .......... 83% 76.2M 0s
    ##  67000K .......... .......... .......... .......... .......... 83%  105M 0s
    ##  67050K .......... .......... .......... .......... .......... 83% 23.6M 0s
    ##  67100K .......... .......... .......... .......... .......... 84% 62.5M 0s
    ##  67150K .......... .......... .......... .......... .......... 84% 52.8M 0s
    ##  67200K .......... .......... .......... .......... .......... 84% 92.9M 0s
    ##  67250K .......... .......... .......... .......... .......... 84% 13.3M 0s
    ##  67300K .......... .......... .......... .......... .......... 84% 59.4M 0s
    ##  67350K .......... .......... .......... .......... .......... 84% 49.5M 0s
    ##  67400K .......... .......... .......... .......... .......... 84% 73.4M 0s
    ##  67450K .......... .......... .......... .......... .......... 84% 82.8M 0s
    ##  67500K .......... .......... .......... .......... .......... 84% 82.7M 0s
    ##  67550K .......... .......... .......... .......... .......... 84% 86.0M 0s
    ##  67600K .......... .......... .......... .......... .......... 84% 74.3M 0s
    ##  67650K .......... .......... .......... .......... .......... 84% 86.0M 0s
    ##  67700K .......... .......... .......... .......... .......... 84% 66.6M 0s
    ##  67750K .......... .......... .......... .......... .......... 84% 97.0M 0s
    ##  67800K .......... .......... .......... .......... .......... 84% 70.8M 0s
    ##  67850K .......... .......... .......... .......... .......... 84% 85.8M 0s
    ##  67900K .......... .......... .......... .......... .......... 85% 74.3M 0s
    ##  67950K .......... .......... .......... .......... .......... 85% 68.7M 0s
    ##  68000K .......... .......... .......... .......... .......... 85% 91.2M 0s
    ##  68050K .......... .......... .......... .......... .......... 85% 49.3M 0s
    ##  68100K .......... .......... .......... .......... .......... 85% 73.5M 0s
    ##  68150K .......... .......... .......... .......... .......... 85% 75.3M 0s
    ##  68200K .......... .......... .......... .......... .......... 85% 89.9M 0s
    ##  68250K .......... .......... .......... .......... .......... 85% 86.9M 0s
    ##  68300K .......... .......... .......... .......... .......... 85% 79.3M 0s
    ##  68350K .......... .......... .......... .......... .......... 85% 71.2M 0s
    ##  68400K .......... .......... .......... .......... .......... 85% 94.9M 0s
    ##  68450K .......... .......... .......... .......... .......... 85% 28.8M 0s
    ##  68500K .......... .......... .......... .......... .......... 85% 69.8M 0s
    ##  68550K .......... .......... .......... .......... .......... 85% 99.9M 0s
    ##  68600K .......... .......... .......... .......... .......... 85% 83.6M 0s
    ##  68650K .......... .......... .......... .......... .......... 85%  112M 0s
    ##  68700K .......... .......... .......... .......... .......... 86% 83.8M 0s
    ##  68750K .......... .......... .......... .......... .......... 86%  107M 0s
    ##  68800K .......... .......... .......... .......... .......... 86%  100M 0s
    ##  68850K .......... .......... .......... .......... .......... 86% 15.4M 0s
    ##  68900K .......... .......... .......... .......... .......... 86% 71.7M 0s
    ##  68950K .......... .......... .......... .......... .......... 86% 67.3M 0s
    ##  69000K .......... .......... .......... .......... .......... 86% 94.9M 0s
    ##  69050K .......... .......... .......... .......... .......... 86% 94.4M 0s
    ##  69100K .......... .......... .......... .......... .......... 86% 75.9M 0s
    ##  69150K .......... .......... .......... .......... .......... 86%  105M 0s
    ##  69200K .......... .......... .......... .......... .......... 86% 57.7M 0s
    ##  69250K .......... .......... .......... .......... .......... 86% 59.3M 0s
    ##  69300K .......... .......... .......... .......... .......... 86% 92.7M 0s
    ##  69350K .......... .......... .......... .......... .......... 86% 63.9M 0s
    ##  69400K .......... .......... .......... .......... .......... 86% 68.6M 0s
    ##  69450K .......... .......... .......... .......... .......... 86%  113M 0s
    ##  69500K .......... .......... .......... .......... .......... 87% 41.8M 0s
    ##  69550K .......... .......... .......... .......... .......... 87% 84.0M 0s
    ##  69600K .......... .......... .......... .......... .......... 87%  110M 0s
    ##  69650K .......... .......... .......... .......... .......... 87% 99.5M 0s
    ##  69700K .......... .......... .......... .......... .......... 87% 67.6M 0s
    ##  69750K .......... .......... .......... .......... .......... 87%  114M 0s
    ##  69800K .......... .......... .......... .......... .......... 87% 98.1M 0s
    ##  69850K .......... .......... .......... .......... .......... 87% 26.4M 0s
    ##  69900K .......... .......... .......... .......... .......... 87%  114M 0s
    ##  69950K .......... .......... .......... .......... .......... 87% 21.7M 0s
    ##  70000K .......... .......... .......... .......... .......... 87% 85.3M 0s
    ##  70050K .......... .......... .......... .......... .......... 87% 92.3M 0s
    ##  70100K .......... .......... .......... .......... .......... 87%  100M 0s
    ##  70150K .......... .......... .......... .......... .......... 87% 97.6M 0s
    ##  70200K .......... .......... .......... .......... .......... 87% 22.8M 0s
    ##  70250K .......... .......... .......... .......... .......... 87% 85.7M 0s
    ##  70300K .......... .......... .......... .......... .......... 88% 71.4M 0s
    ##  70350K .......... .......... .......... .......... .......... 88% 57.2M 0s
    ##  70400K .......... .......... .......... .......... .......... 88% 68.6M 0s
    ##  70450K .......... .......... .......... .......... .......... 88% 79.7M 0s
    ##  70500K .......... .......... .......... .......... .......... 88%  101M 0s
    ##  70550K .......... .......... .......... .......... .......... 88% 15.4M 0s
    ##  70600K .......... .......... .......... .......... .......... 88% 72.9M 0s
    ##  70650K .......... .......... .......... .......... .......... 88%  139M 0s
    ##  70700K .......... .......... .......... .......... .......... 88%  105M 0s
    ##  70750K .......... .......... .......... .......... .......... 88%  143M 0s
    ##  70800K .......... .......... .......... .......... .......... 88%  119M 0s
    ##  70850K .......... .......... .......... .......... .......... 88%  146M 0s
    ##  70900K .......... .......... .......... .......... .......... 88% 65.2M 0s
    ##  70950K .......... .......... .......... .......... .......... 88% 27.7M 0s
    ##  71000K .......... .......... .......... .......... .......... 88% 65.7M 0s
    ##  71050K .......... .......... .......... .......... .......... 88% 75.6M 0s
    ##  71100K .......... .......... .......... .......... .......... 89% 79.5M 0s
    ##  71150K .......... .......... .......... .......... .......... 89%  109M 0s
    ##  71200K .......... .......... .......... .......... .......... 89% 40.7M 0s
    ##  71250K .......... .......... .......... .......... .......... 89% 89.8M 0s
    ##  71300K .......... .......... .......... .......... .......... 89%  117M 0s
    ##  71350K .......... .......... .......... .......... .......... 89% 12.7M 0s
    ##  71400K .......... .......... .......... .......... .......... 89% 98.9M 0s
    ##  71450K .......... .......... .......... .......... .......... 89%  137M 0s
    ##  71500K .......... .......... .......... .......... .......... 89%  119M 0s
    ##  71550K .......... .......... .......... .......... .......... 89%  118M 0s
    ##  71600K .......... .......... .......... .......... .......... 89%  105M 0s
    ##  71650K .......... .......... .......... .......... .......... 89%  130M 0s
    ##  71700K .......... .......... .......... .......... .......... 89% 82.3M 0s
    ##  71750K .......... .......... .......... .......... .......... 89% 31.8M 0s
    ##  71800K .......... .......... .......... .......... .......... 89% 34.7M 0s
    ##  71850K .......... .......... .......... .......... .......... 89% 97.5M 0s
    ##  71900K .......... .......... .......... .......... .......... 90%  110M 0s
    ##  71950K .......... .......... .......... .......... .......... 90%  134M 0s
    ##  72000K .......... .......... .......... .......... .......... 90%  116M 0s
    ##  72050K .......... .......... .......... .......... .......... 90% 33.3M 0s
    ##  72100K .......... .......... .......... .......... .......... 90% 94.7M 0s
    ##  72150K .......... .......... .......... .......... .......... 90% 89.9M 0s
    ##  72200K .......... .......... .......... .......... .......... 90% 91.6M 0s
    ##  72250K .......... .......... .......... .......... .......... 90% 60.0M 0s
    ##  72300K .......... .......... .......... .......... .......... 90% 42.5M 0s
    ##  72350K .......... .......... .......... .......... .......... 90%  118M 0s
    ##  72400K .......... .......... .......... .......... .......... 90% 91.6M 0s
    ##  72450K .......... .......... .......... .......... .......... 90%  127M 0s
    ##  72500K .......... .......... .......... .......... .......... 90% 98.6M 0s
    ##  72550K .......... .......... .......... .......... .......... 90%  112M 0s
    ##  72600K .......... .......... .......... .......... .......... 90%  115M 0s
    ##  72650K .......... .......... .......... .......... .......... 90%  128M 0s
    ##  72700K .......... .......... .......... .......... .......... 91% 32.1M 0s
    ##  72750K .......... .......... .......... .......... .......... 91% 48.7M 0s
    ##  72800K .......... .......... .......... .......... .......... 91% 88.6M 0s
    ##  72850K .......... .......... .......... .......... .......... 91% 86.2M 0s
    ##  72900K .......... .......... .......... .......... .......... 91% 89.1M 0s
    ##  72950K .......... .......... .......... .......... .......... 91% 62.2M 0s
    ##  73000K .......... .......... .......... .......... .......... 91% 24.9M 0s
    ##  73050K .......... .......... .......... .......... .......... 91%  127M 0s
    ##  73100K .......... .......... .......... .......... .......... 91%  103M 0s
    ##  73150K .......... .......... .......... .......... .......... 91%  119M 0s
    ##  73200K .......... .......... .......... .......... .......... 91% 29.0M 0s
    ##  73250K .......... .......... .......... .......... .......... 91%  127M 0s
    ##  73300K .......... .......... .......... .......... .......... 91% 88.3M 0s
    ##  73350K .......... .......... .......... .......... .......... 91%  133M 0s
    ##  73400K .......... .......... .......... .......... .......... 91%  107M 0s
    ##  73450K .......... .......... .......... .......... .......... 91% 80.5M 0s
    ##  73500K .......... .......... .......... .......... .......... 92% 97.9M 0s
    ##  73550K .......... .......... .......... .......... .......... 92%  132M 0s
    ##  73600K .......... .......... .......... .......... .......... 92% 42.5M 0s
    ##  73650K .......... .......... .......... .......... .......... 92% 74.7M 0s
    ##  73700K .......... .......... .......... .......... .......... 92% 89.8M 0s
    ##  73750K .......... .......... .......... .......... .......... 92% 59.6M 0s
    ##  73800K .......... .......... .......... .......... .......... 92%  103M 0s
    ##  73850K .......... .......... .......... .......... .......... 92%  116M 0s
    ##  73900K .......... .......... .......... .......... .......... 92% 71.3M 0s
    ##  73950K .......... .......... .......... .......... .......... 92%  101M 0s
    ##  74000K .......... .......... .......... .......... .......... 92% 93.7M 0s
    ##  74050K .......... .......... .......... .......... .......... 92% 64.5M 0s
    ##  74100K .......... .......... .......... .......... .......... 92% 54.4M 0s
    ##  74150K .......... .......... .......... .......... .......... 92% 44.6M 0s
    ##  74200K .......... .......... .......... .......... .......... 92% 33.4M 0s
    ##  74250K .......... .......... .......... .......... .......... 92%  135M 0s
    ##  74300K .......... .......... .......... .......... .......... 93%  108M 0s
    ##  74350K .......... .......... .......... .......... .......... 93%  121M 0s
    ##  74400K .......... .......... .......... .......... .......... 93% 67.6M 0s
    ##  74450K .......... .......... .......... .......... .......... 93%  121M 0s
    ##  74500K .......... .......... .......... .......... .......... 93% 49.6M 0s
    ##  74550K .......... .......... .......... .......... .......... 93% 58.9M 0s
    ##  74600K .......... .......... .......... .......... .......... 93% 31.4M 0s
    ##  74650K .......... .......... .......... .......... .......... 93%  102M 0s
    ##  74700K .......... .......... .......... .......... .......... 93% 24.8M 0s
    ##  74750K .......... .......... .......... .......... .......... 93%  141M 0s
    ##  74800K .......... .......... .......... .......... .......... 93% 36.2M 0s
    ##  74850K .......... .......... .......... .......... .......... 93%  123M 0s
    ##  74900K .......... .......... .......... .......... .......... 93%  118M 0s
    ##  74950K .......... .......... .......... .......... .......... 93%  114M 0s
    ##  75000K .......... .......... .......... .......... .......... 93%  122M 0s
    ##  75050K .......... .......... .......... .......... .......... 93%  137M 0s
    ##  75100K .......... .......... .......... .......... .......... 94% 83.2M 0s
    ##  75150K .......... .......... .......... .......... .......... 94%  109M 0s
    ##  75200K .......... .......... .......... .......... .......... 94%  109M 0s
    ##  75250K .......... .......... .......... .......... .......... 94% 29.9M 0s
    ##  75300K .......... .......... .......... .......... .......... 94% 70.5M 0s
    ##  75350K .......... .......... .......... .......... .......... 94% 22.7M 0s
    ##  75400K .......... .......... .......... .......... .......... 94%  107M 0s
    ##  75450K .......... .......... .......... .......... .......... 94%  132M 0s
    ##  75500K .......... .......... .......... .......... .......... 94% 99.3M 0s
    ##  75550K .......... .......... .......... .......... .......... 94%  100M 0s
    ##  75600K .......... .......... .......... .......... .......... 94%  137M 0s
    ##  75650K .......... .......... .......... .......... .......... 94%  121M 0s
    ##  75700K .......... .......... .......... .......... .......... 94% 79.8M 0s
    ##  75750K .......... .......... .......... .......... .......... 94% 41.3M 0s
    ##  75800K .......... .......... .......... .......... .......... 94% 61.5M 0s
    ##  75850K .......... .......... .......... .......... .......... 94% 57.2M 0s
    ##  75900K .......... .......... .......... .......... .......... 95% 87.2M 0s
    ##  75950K .......... .......... .......... .......... .......... 95% 25.8M 0s
    ##  76000K .......... .......... .......... .......... .......... 95%  104M 0s
    ##  76050K .......... .......... .......... .......... .......... 95% 41.0M 0s
    ##  76100K .......... .......... .......... .......... .......... 95% 80.2M 0s
    ##  76150K .......... .......... .......... .......... .......... 95%  110M 0s
    ##  76200K .......... .......... .......... .......... .......... 95%  144M 0s
    ##  76250K .......... .......... .......... .......... .......... 95% 51.3M 0s
    ##  76300K .......... .......... .......... .......... .......... 95% 73.3M 0s
    ##  76350K .......... .......... .......... .......... .......... 95%  124M 0s
    ##  76400K .......... .......... .......... .......... .......... 95%  114M 0s
    ##  76450K .......... .......... .......... .......... .......... 95% 46.9M 0s
    ##  76500K .......... .......... .......... .......... .......... 95%  115M 0s
    ##  76550K .......... .......... .......... .......... .......... 95% 82.9M 0s
    ##  76600K .......... .......... .......... .......... .......... 95% 40.2M 0s
    ##  76650K .......... .......... .......... .......... .......... 95% 54.2M 0s
    ##  76700K .......... .......... .......... .......... .......... 96%  110M 0s
    ##  76750K .......... .......... .......... .......... .......... 96% 60.3M 0s
    ##  76800K .......... .......... .......... .......... .......... 96%  107M 0s
    ##  76850K .......... .......... .......... .......... .......... 96%  128M 0s
    ##  76900K .......... .......... .......... .......... .......... 96% 98.7M 0s
    ##  76950K .......... .......... .......... .......... .......... 96%  130M 0s
    ##  77000K .......... .......... .......... .......... .......... 96% 91.3M 0s
    ##  77050K .......... .......... .......... .......... .......... 96% 42.5M 0s
    ##  77100K .......... .......... .......... .......... .......... 96% 84.5M 0s
    ##  77150K .......... .......... .......... .......... .......... 96% 71.6M 0s
    ##  77200K .......... .......... .......... .......... .......... 96% 38.0M 0s
    ##  77250K .......... .......... .......... .......... .......... 96%  116M 0s
    ##  77300K .......... .......... .......... .......... .......... 96% 78.1M 0s
    ##  77350K .......... .......... .......... .......... .......... 96%  112M 0s
    ##  77400K .......... .......... .......... .......... .......... 96% 89.9M 0s
    ##  77450K .......... .......... .......... .......... .......... 96% 57.3M 0s
    ##  77500K .......... .......... .......... .......... .......... 97% 81.1M 0s
    ##  77550K .......... .......... .......... .......... .......... 97%  115M 0s
    ##  77600K .......... .......... .......... .......... .......... 97% 45.8M 0s
    ##  77650K .......... .......... .......... .......... .......... 97% 29.5M 0s
    ##  77700K .......... .......... .......... .......... .......... 97% 83.0M 0s
    ##  77750K .......... .......... .......... .......... .......... 97%  128M 0s
    ##  77800K .......... .......... .......... .......... .......... 97% 95.8M 0s
    ##  77850K .......... .......... .......... .......... .......... 97% 64.6M 0s
    ##  77900K .......... .......... .......... .......... .......... 97%  109M 0s
    ##  77950K .......... .......... .......... .......... .......... 97% 87.8M 0s
    ##  78000K .......... .......... .......... .......... .......... 97% 87.3M 0s
    ##  78050K .......... .......... .......... .......... .......... 97%  105M 0s
    ##  78100K .......... .......... .......... .......... .......... 97% 42.4M 0s
    ##  78150K .......... .......... .......... .......... .......... 97%  103M 0s
    ##  78200K .......... .......... .......... .......... .......... 97% 63.0M 0s
    ##  78250K .......... .......... .......... .......... .......... 97% 56.0M 0s
    ##  78300K .......... .......... .......... .......... .......... 98% 53.1M 0s
    ##  78350K .......... .......... .......... .......... .......... 98%  110M 0s
    ##  78400K .......... .......... .......... .......... .......... 98%  130M 0s
    ##  78450K .......... .......... .......... .......... .......... 98% 64.4M 0s
    ##  78500K .......... .......... .......... .......... .......... 98% 98.4M 0s
    ##  78550K .......... .......... .......... .......... .......... 98% 95.1M 0s
    ##  78600K .......... .......... .......... .......... .......... 98% 85.5M 0s
    ##  78650K .......... .......... .......... .......... .......... 98% 57.8M 0s
    ##  78700K .......... .......... .......... .......... .......... 98% 27.3M 0s
    ##  78750K .......... .......... .......... .......... .......... 98% 71.5M 0s
    ##  78800K .......... .......... .......... .......... .......... 98% 39.1M 0s
    ##  78850K .......... .......... .......... .......... .......... 98% 16.5M 0s
    ##  78900K .......... .......... .......... .......... .......... 98%  123M 0s
    ##  78950K .......... .......... .......... .......... .......... 98%  153M 0s
    ##  79000K .......... .......... .......... .......... .......... 98%  135M 0s
    ##  79050K .......... .......... .......... .......... .......... 98%  109M 0s
    ##  79100K .......... .......... .......... .......... .......... 99%  103M 0s
    ##  79150K .......... .......... .......... .......... .......... 99%  126M 0s
    ##  79200K .......... .......... .......... .......... .......... 99%  104M 0s
    ##  79250K .......... .......... .......... .......... .......... 99% 47.4M 0s
    ##  79300K .......... .......... .......... .......... .......... 99% 99.4M 0s
    ##  79350K .......... .......... .......... .......... .......... 99%  116M 0s
    ##  79400K .......... .......... .......... .......... .......... 99% 37.1M 0s
    ##  79450K .......... .......... .......... .......... .......... 99% 93.7M 0s
    ##  79500K .......... .......... .......... .......... .......... 99% 99.4M 0s
    ##  79550K .......... .......... .......... .......... .......... 99%  117M 0s
    ##  79600K .......... .......... .......... .......... .......... 99%  116M 0s
    ##  79650K .......... .......... .......... .......... .......... 99% 94.0M 0s
    ##  79700K .......... .......... .......... .......... .......... 99%  143M 0s
    ##  79750K .......... .......... .......... .......... .......... 99% 97.3M 0s
    ##  79800K .......... .......... .......... .......... .......... 99% 36.4M 0s
    ##  79850K .......... .......... .......... .......... .......... 99%  134M 0s
    ##  79900K .......... .......... ..                              100%  136M=1.6s
    ## 
    ## 2020-12-04 14:29:04 (49.5 MB/s) - ‘silva_species_assignment_v138.fa.gz.3’ saved [81840166/81840166]

Importation des données d’espèces pour la taxonomie.

``` r
taxa <- addSpecies(taxa, "~/Dada2_Phyloseq_CC1/silva_species_assignment_v138.fa.gz")
```

Assignation des espèces à la taxonomie précédement créée.

``` r
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum         Class         Order           Family          
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      Genus         Species
    ## [1,] NA            NA     
    ## [2,] NA            NA     
    ## [3,] NA            NA     
    ## [4,] NA            NA     
    ## [5,] "Bacteroides" NA     
    ## [6,] NA            NA

L’assignation taxonomique donne le tableau ci-dessus. Elle donne
jusqu’au genre ce qu’elle a pu identifier des reads à analyser, en les
comparant avec les reads contenus dans la banque de Silva. Aucune espèce
n’a pu être identifiée à partir des reads.

## Evaluation de la précisison

``` r
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

    ## DADA2 inferred 20 sample sequences present in the Mock community.

“Mock community” est un ensemble contenant différentes souches connues
qui ont été séquencées. Cette communauté contient 20 souches
différentes. On utilise Dada2 sur cette communauté pour vérifier sa
précision. Ici, dada2 a déduit 20 échantillons avec les séquences, ce
qui est bien.

``` r
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

    ## Of those, 20 were exact matches to the expected reference sequences.

Parmi ces 20 échantillons, 20 correspondaient aux séquences de
références. Dada2 fonctionne bien.

# Premier pas sur Phyloseq (fin tuto Dada2)

## Importation des librairies et préparation des données

``` r
library(phyloseq)
library(Biostrings)
```

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
library(ggplot2)
```

``` r
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
```

On assigne les différentes données dans des fichiers.

``` r
library("phyloseq")
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps)
```

Construction de l’objet PhyloSeq : on associe la table des OTU, les
données des échantillons et de la table taxonomique à ps. La dernière
ligne permet d’enlever la communité artificielle (“mock”)

``` r
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 232 taxa and 19 samples ]
    ## sample_data() Sample Data:       [ 19 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 232 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 232 reference sequences ]

Visualisation de l’objet PhyloSeq

## Visualisation de l’alpha diversité

``` r
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
```

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->
Pas de grandes différences entre les échantillons contenus dans
prélévements “early” (entre eux) et les prélévemnts “late” (entre aux
aussi). Au fil des jours, la diversité n’a pas variée énormement.

## Ordination

``` r
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
```

    ## Run 0 stress 0.08574537 
    ## Run 1 stress 0.08574538 
    ## ... Procrustes: rmse 2.858155e-05  max resid 9.538668e-05 
    ## ... Similar to previous best
    ## Run 2 stress 0.08002299 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04283637  max resid 0.1433496 
    ## Run 3 stress 0.08574537 
    ## Run 4 stress 0.08002299 
    ## ... New best solution
    ## ... Procrustes: rmse 1.103104e-06  max resid 2.285813e-06 
    ## ... Similar to previous best
    ## Run 5 stress 0.08002299 
    ## ... Procrustes: rmse 1.137148e-06  max resid 2.334528e-06 
    ## ... Similar to previous best
    ## Run 6 stress 0.08942878 
    ## Run 7 stress 0.09421601 
    ## Run 8 stress 0.1216672 
    ## Run 9 stress 0.08942874 
    ## Run 10 stress 0.08002299 
    ## ... Procrustes: rmse 2.732367e-06  max resid 8.078333e-06 
    ## ... Similar to previous best
    ## Run 11 stress 0.08002299 
    ## ... Procrustes: rmse 1.50819e-06  max resid 2.757312e-06 
    ## ... Similar to previous best
    ## Run 12 stress 0.09421618 
    ## Run 13 stress 0.08002299 
    ## ... New best solution
    ## ... Procrustes: rmse 8.105378e-07  max resid 1.702962e-06 
    ## ... Similar to previous best
    ## Run 14 stress 0.08942865 
    ## Run 15 stress 0.08574537 
    ## Run 16 stress 0.08002299 
    ## ... Procrustes: rmse 1.070942e-05  max resid 3.737986e-05 
    ## ... Similar to previous best
    ## Run 17 stress 0.08942909 
    ## Run 18 stress 0.08942891 
    ## Run 19 stress 0.08942868 
    ## Run 20 stress 0.08002299 
    ## ... Procrustes: rmse 5.230355e-06  max resid 1.291396e-05 
    ## ... Similar to previous best
    ## *** Solution reached

On transforme les données pour qu’ils soient adaptés à l’utilisation de
Bray-Curtis, un indice de dissimilarité.

``` r
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")
```

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->
A partir de ce graphe, on peut une grande différence entre les
échantillons “early” et “late” (dans leur communauté microbienne)

## Représention des échantillons en “histogramme” (bar plot)

``` r
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
```

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->
On observe qu’il n’y a pas de grande différence entre les échantillons
“early” et “late” au niveau des familles de micro-organismes. Il n’y
pas de grandes différences sur la diversité ni sur les familles
bactériennes. Mais on observe une nette séparation entre les deux
échantillons. Cette différence est peut-être dû aux espèces présentes
dans les familles.
