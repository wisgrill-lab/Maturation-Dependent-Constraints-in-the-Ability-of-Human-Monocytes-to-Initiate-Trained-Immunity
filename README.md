

Code repository for the publication "Maturation-Dependent Constraints in the Ability of Human Monocytes to Initiate Trained Immunity"
================
<p><b>Analysis: Michael Eigenschink MD</b></p>
2024-08-17
<p><i>Medical University of Vienna</i></p>

<p>
  Michael Eigenschink<sup>1</sup>, Lukas Wisgrill<sup>1</sup><br></p>

  <p>
  <sup>1</sup>Division of Neonatology, Pediatric Intensive Care & Neuropaediatrics, Department of Pediatrics and Adolescent Medicine, Comprehensive Center for Pediatrics, Medical University of Vienna, Austria </p>




# Analysis of LC-MS data

## Import of data and initial wrangling

``` r
#Metabolomics (HILIC)
metabolome_HILIC <- rio::import(here::here("data","Results_Hilic_Manuel.csv"))

metabolome_HILIC <- metabolome_HILIC %>% mutate(Group = case_when(grepl("A", V1) ~ "Adult", grepl("P", V1) ~"Preterm", grepl("T", V1)~"Term"))
metabolome_HILIC <- metabolome_HILIC %>% mutate(Treatment = case_when(grepl("M1", V1) ~ "control", grepl("M2", V1) ~"Glucan"))


metabolome_HILIC <- metabolome_HILIC %>% mutate_at(c(2:49), as.numeric)
metabolome_HILIC <- metabolome_HILIC %>% mutate_at(c(50:51), as.factor)
metabolome_HILIC <- metabolome_HILIC %>% mutate_at(c(1), as.character)


#Metabolomics (Reverse Phase)
metabolome_RP <- rio::import(here::here("data","Results_RP_Manuel.csv"))

metabolome_RP <- metabolome_RP %>% mutate(Group = case_when(grepl("A", V1) ~ "Adult", grepl("P", V1) ~"Preterm", grepl("T", V1)~"Term"))
metabolome_RP <- metabolome_RP %>% mutate(Treatment = case_when(grepl("M1", V1) ~ "control", grepl("M2", V1) ~"Glucan"))

metabolome_RP <- metabolome_RP %>% mutate_at(c(2:20), as.numeric)
metabolome_RP <- metabolome_RP %>% mutate_at(c(21:22), as.factor)
metabolome_RP <- metabolome_RP %>% mutate_at(c(1), as.character)

#Which measurements are unique to reverse phase? (Mass-Spec specialists Univ. of Vienna: "HILIC gave more reliable results for compounds found in both columns")
metabolome_cols_RP <- metabolome_RP %>% dplyr::select(-V1, -Group, -Treatment) %>%  colnames()
metabolome_cols_HILIC <- metabolome_HILIC %>% dplyr::select(-V1, -Group, -Treatment) %>%  colnames()
difference <- setdiff(metabolome_cols_RP, metabolome_cols_HILIC)

metabolome_RP <- metabolome_RP %>% pivot_longer(cols = 2:20, names_to = "metabolites", values_to = "AUC")
metabolome_RP <- subset(metabolome_RP, (metabolome_RP$metabolites %in% difference))
metabolome_RP <- metabolome_RP %>%  pivot_wider(names_from = "metabolites", values_from = "AUC")
metabolome_RP <- metabolome_RP %>% dplyr::select(-Group, -Treatment)

#how many unique metabolites were identified?
vec1 <- colnames(metabolome_HILIC)
vec2 <- colnames(metabolome_RP)

vecvec <- c(vec1, vec2)
length(unique(vecvec))
```

    ## [1] 56

``` r
#merge data
metabolome <- metabolome_HILIC %>% left_join(metabolome_RP, by = "V1")
metabolome <- metabolome %>% dplyr::select(1:49, 52:56, 50, 51)

#NA_distances -> drop spermidine
#metabolome <- metabolome_HILIC %>% dplyr::select(-Spermidine)
```

### Figure 3A: Multilevel PCA

``` r
metabolome$ID <- c("1", "1", "2", "2", "3", "3", "4","4", "5", "5", "6", "6", "7", "7", "8", "8","9", "9", "10", "10", "11", "11", "12", "12")
rownames(metabolome) <- metabolome$V1 
metabolome$ID <- as.factor(metabolome$ID)
metabolome <- metabolome %>% mutate(index = paste0(Group, Treatment, SEP = ""))
metabolome_2 <- metabolome %>% dplyr::select(-V1)

metabolome$index <- as.factor(metabolome$index)
metabolome_tr <- metabolome_2 [1:53]

pca.result <- pca(metabolome_tr, multilevel = metabolome_2$ID, scale = TRUE)
```

    ## Splitting the variation for 1 level factor.

``` r
bipot_multilevel_PCA <- biplot(pca.result, cutoff =0.85, pch = metabolome$Treatment, ind.names = TRUE, 
       group = metabolome$Group, ellipse = TRUE, col.per.group = c("#b2c5b3", "#CC6677", "#6699CC"), ind.names.size = 4, pch.size   
= 2.75, var.arrow.size  = 0.5, var.arrow.length = 0.3, vline=TRUE, hline    = TRUE)


graph2svg(bipot_multilevel_PCA, file = here::here("plots","bipot_multilevel_PCA"),  width = 4.5, height = 2.45)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/bipot_multilevel_PCA.svg

### Figure 3B: Boxplots

``` r
metabolome_2$indiv <- c("1", "1", "2","2", "3","3", "4","4", "1","1", "2", "2", "3", "3", "1", "1", "2", "2", "3", "3", "4", "4", "5", "5")
metabolome_2$indiv <- as.factor(metabolome_2$indiv)
 
metabolome_2 <- metabolome_2 %>% dplyr::select(-ID, -index)
metabolome_3 <- metabolome_2 %>%  pivot_longer(col = 1:53, names_to = "metabolites", values_to = "AUC")
 
metabolome_3 <- metabolome_3 %>% mutate(new_col = paste0(Group, Treatment, SEP = ""))
metabolome_3$new_col <- as.factor(metabolome_3$new_col)
 
metabolites_Figure <- ggplot(metabolome_3, aes(x = new_col, y = AUC, fill = Group)) + scale_x_discrete(limits = c("Adultcontrol", "AdultGlucan", "Termcontrol", "TermGlucan", "Pretermcontrol", "PretermGlucan")) +
   geom_boxplot()  +
   geom_line(aes(group = interaction(Group, indiv), col = Group),
             alpha = 0.4) + geom_point(aes(group = interaction(Group, indiv), col = Group), alpha = 0.6) + 
   scale_fill_manual(values=c("#97C591", "#ce8793", "#649CD1")) +
   scale_color_manual(values=c("#97C591", "#ce8793", "#649CD1")) +
   facet_wrap(~metabolites, scales ="free") + ylab("AUC") + xlab("Condition") + theme_bw() +
  theme(axis.text=element_text(size=8), 
        axis.title = element_text(size = 10,  face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 10))

#all names of metabolites from figure
vector_metabo_Figure <- c("L-Ornithine", "Arginine", "Proline", 
     "Tryptophan", "Dihydroxyacetonephosphate", 
     "Serine", "Glycine", "Methionine", "Pyruvate",
     "Glutamate", "Glutamine","alpha-Ketoglutarate",
     "Succinate", "Fumarate", "Malate", "Aspartate", "Lactate")

metabolome_supplement <- subset(metabolome_3,!(metabolome_3$metabolites %in% vector_metabo_Figure))

metabolites_Figure_supplement <- ggplot(metabolome_supplement, aes(x = new_col, y = AUC, fill = Group)) + scale_x_discrete(limits = c("Adultcontrol", "AdultGlucan", "Termcontrol", "TermGlucan", "Pretermcontrol", "PretermGlucan")) +
   geom_boxplot()  +
   geom_line(aes(group = interaction(Group, indiv), col = Group),
             alpha = 0.4) + geom_point(aes(group = interaction(Group, indiv), col = Group), alpha = 0.6) + 
   scale_fill_manual(values=c("#97C591", "#ce8793", "#649CD1")) +
   scale_color_manual(values=c("#97C591", "#ce8793", "#649CD1")) +
   facet_wrap(~metabolites, scales ="free") + ylab("AUC") + xlab("Condition") + theme_bw() +
  theme(axis.text=element_text(size=8), 
        axis.title = element_text(size = 10,  face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 10))

graph2svg(metabolites_Figure_supplement, file = here::here("plots","metabolites_Figure_supplement"),  width = 10, height = 9)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/metabolites_Figure_supplement.svg

## Statistical analysis: Figure 2

``` r
metabolome_2 <- metabolome_2 %>%  pivot_longer(col = 1:53, names_to = "metabolites", values_to = "AUC")
metabo_analyst <- subset(metabolome_2, metabolites %in% vector_metabo_Figure)

x <- metabo_analyst %>% 
  group_by(Group, Treatment, metabolites) %>%
  identify_outliers(AUC)

x <- x %>% filter(is.extreme == TRUE)
exclude <- x$metabolites %>% unique()

metabo_analyst <- metabo_analyst %>% mutate(exclude = case_when(metabolites == "Methionine" & Group == "Term" & indiv == 4 ~ "exclude",
                                                                        metabolites == "Fumarate" & Group == "Term" & indiv == 4 ~ "exclude",
                                                                        metabolites == "Glutamine" & Group == "Term" & indiv == 4 ~ "exclude",
                                                                        metabolites == "Glycine" & Group == "Term" & indiv == 4 ~ "exclude",
                                                                        metabolites == "L-Ornithine" & Group == "Term" & indiv == 4 ~ "exclude",
                                                                        metabolites == "Proline" & Group == "Term" & indiv == 4 ~ "exclude",
                                                                        metabolites == "Pyruvate" & Group == "Term" & indiv == 4 ~ "exclude",
                                                                        TRUE ~ "keep"))

metabolome_filtered <- metabo_analyst %>% filter(exclude == "keep")
metabo_analyst <- subset(metabolome_filtered, metabolites %in% vector_metabo_Figure)

#sanity check
x <- metabo_analyst %>% 
  group_by(Group, Treatment, metabolites) %>%
  identify_outliers(AUC) %>% ungroup()

#check normality
metabo_analyst <- metabo_analyst %>% mutate(log_AUC = log10(AUC))
norm.assum <- metabo_analyst %>%
  group_by(Group, Treatment, metabolites) %>%
  shapiro_test(log_AUC)

#Ornithine
anova <- metabo_analyst %>% filter(metabolites == "L-Ornithine") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)
a <- get_anova_table(res.aov)


#nothing

#Arginine
anova <- metabo_analyst %>% filter(metabolites == "Arginine") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)
a <- get_anova_table(res.aov)

# Group effect, p = 0.031
# Group differences across treatment
test2 <- anova %>%
  group_by(Treatment) %>%
  pairwise_t_test(
    AUC ~ Group, paired = FALSE,
    p.adjust.method = "BH"
    )
print(test2)
a <- test2
#no significant comparisons

#Proline
anova <- metabo_analyst %>% filter(metabolites == "Proline") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)
a <- get_anova_table(res.aov)

#nothing significant

#Tryptophan
anova <- metabo_analyst %>% filter(metabolites == "Tryptophan") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)
a <- get_anova_table(res.aov)
#nothing significant

#Serine
anova <- metabo_analyst %>% filter(metabolites == "Serine") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)
a <- get_anova_table(res.aov)

#nothing significant

#Glycine
anova <- metabo_analyst %>% filter(metabolites == "Glycine") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)
a <- get_anova_table(res.aov)

#nothing significant

#Methionine
anova <- metabo_analyst %>% filter(metabolites == "Methionine") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)
a <- get_anova_table(res.aov)


# Treatment effect, p = 0.010
# Effect of treatment across groups
test <- anova %>%
  group_by(Group) %>%
  pairwise_t_test(
    AUC ~ Treatment, paired = TRUE,
    p.adjust.method = "BH")
print(test)
#no significant comparisons
a <- test

#Pyruvate
anova <- metabo_analyst %>% filter(metabolites == "Pyruvate") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)
a <- get_anova_table(res.aov)


#nothing significant

#Glutamate
anova <- metabo_analyst %>% filter(metabolites == "Glutamate") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)
a <- get_anova_table(res.aov)


# Treatment effect, p = 0.032
# Effect of treatment across groups
test <- anova %>%
  group_by(Group) %>%
  pairwise_t_test(
    AUC ~ Treatment, paired = TRUE,
    p.adjust.method = "BH")
print(test)
a <- test
# Significant treatment effect in preterm infants, p = 0.042 after BH

#Glutamine
anova <- metabo_analyst %>% filter(metabolites == "Glutamine") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)
a <- get_anova_table(res.aov)

# Treatment effect, p = 0.036
# Effect of treatment across groups
test <- anova %>%
  group_by(Group) %>%
  pairwise_t_test(
    AUC ~ Treatment, paired = TRUE,
    p.adjust.method = "BH")
print(test)

a <- test

# Significant treatment effect in adults, p = 0.000257 after BH

#alpha-Ketoglutarate
anova <- metabo_analyst %>% filter(metabolites == "alpha-Ketoglutarate") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)

a <- get_anova_table(res.aov)
#nothing significant

#Succinate
anova <- metabo_analyst %>% filter(metabolites == "Succinate") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)

a <- get_anova_table(res.aov)

#nothing significant

#Fumarate
anova <- metabo_analyst %>% filter(metabolites == "Fumarate") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)

a <- get_anova_table(res.aov)
#nothing significant

#Malate
anova <- metabo_analyst %>% filter(metabolites == "Malate") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)

a <- get_anova_table(res.aov)
#Interaction Treatment:Group significant, p = 0.038

# Effect of treatment across groups
test <- anova %>%
  group_by(Group) %>%
  pairwise_t_test(
    AUC ~ Treatment, paired = TRUE,
    p.adjust.method = "BH")
print(test)
# nothing significant
a <- test
# Group differences across treatment
test2 <- anova %>%
  group_by(Treatment) %>%
  pairwise_t_test(
    AUC ~ Group, paired = FALSE,
    p.adjust.method = "BH"
    )
print(test2)
a <- test2
# No significant differences remain after BH-correction

#Aspartate
anova <- metabo_analyst %>% filter(metabolites == "Aspartate") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)

a <- get_anova_table(res.aov)
#Treatment effect, p = 0.039

# Effect of treatment across groups
test <- anova %>%
  group_by(Group) %>%
  pairwise_t_test(
    AUC ~ Treatment, paired = TRUE,
    p.adjust.method = "BH")
print(test)

a <- test
# Significant treatment effect in adults, p = 0.018 after BH
# Significant treatment effect in preterm infants, p = 0.005 after BH

#Lactate
anova <- metabo_analyst %>% filter(metabolites == "Lactate") %>% as.data.frame() 
anova$indiv <- as.factor(anova$indiv)
res.aov <- anova_test(data = anova, dv = log_AUC, wid = indiv, within = c(Treatment, Group))
get_anova_table(res.aov)
a <- get_anova_table(res.aov)

#Treatment effect, p = 0.011

# Effect of treatment across groups
test <- anova %>%
  group_by(Group) %>%
  pairwise_t_test(
    AUC ~ Treatment, paired = TRUE,
    p.adjust.method = "BH")
print(test)

#Significant treatment effect in adults, p = 0.026 after BH
```

# Paired-End RNAseq analysis

## Import of count matrix & metadata + DGEList

``` r
metadata <- readRDS(here::here("data","metadata.rds"))
se <- readRDS(here::here("data","se.rds"))

#add metadata
colData(se) <- DataFrame(metadata)
colData(se)

#assign rownames explicitely
rownames(colData(se)) <- colData(se)$SampleName

#convert data to DGEList
x <- SE2DGEList(se)
```

## Check library size

``` r
barplot(x$samples$lib.size, names=colnames(x), las=2)
```

![](Trained_immunity_GIT_files/figure-gfm/setup22-1.png)<!-- -->

``` r
#different library sizes across samples
```

## Exclude sample that has been identified as an outlier

``` r
x <- x[,!x$samples$Index == "11"]
```

### Wrangling, filtering & normalization

``` r
#wrangling - recode as factor
x[["samples"]][["SampleName"]] <- x[["samples"]][["SampleName"]] %>% as.factor()
x[["samples"]][["Group"]] <- x[["samples"]][["Group"]] %>% as.factor()
x[["samples"]][["Treatment"]] <- x[["samples"]][["Treatment"]] %>% as.factor()

#wrangling - add new grouping variable to the dataset
x[["samples"]] <- x[["samples"]] %>% mutate(new_col = paste0(Group, Treatment, SEP = ""))
x[["samples"]][["new_col"]] <- as.factor(x[["samples"]][["new_col"]])

#use filterByExpr to filter all genes <10 in at least n = smalles groupsize 
keep.exprs <- edgeR::filterByExpr(x, group=x[["samples"]][["new_col"]])
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
```

    ## [1] 16881    26

``` r
#adjust for different library sizes across samples using TMM
x <- calcNormFactors(x, method = "TMM")
```

## Limma analysis for nested multilevel before-after design

### Creation of model matrix & dupcor + voom iteration

``` r
#in line with recommendations in the "Limma" workflow and recommendations by the limma authors use vomm + dupcor and iterate twice
#define factors for model and blocking
Treat <- factor(x[["samples"]][["new_col"]])
Index <- factor(x[["samples"]][["Index"]])

#design matrix
design = model.matrix( ~ 0 + Treat)
colnames(design) <- levels(Treat)

#voom first round
vobj_tmp = voom(x, design, plot=TRUE)
```

![](Trained_immunity_GIT_files/figure-gfm/setup25-1.png)<!-- -->

``` r
#dupcor first round
dupcor <- duplicateCorrelation(vobj_tmp,design,block=Index)

#voom second round
vobj = voom(x, design, plot=TRUE, block=Index, correlation=dupcor$consensus)
```

![](Trained_immunity_GIT_files/figure-gfm/setup25-2.png)<!-- -->

``` r
#dupcor second round
dupcor <- duplicateCorrelation(vobj, design, block=Index)

#lmFit
fitDupCor <- lmFit(vobj, design, block=Index, correlation=dupcor$consensus)
```

### Contrasts

``` r
cm <- makeContrasts(
  Adult_Effect = adultbetaglucan-adultmock,
  Preterm_Effect = pretermbetaglucan-pretermmock,
  Term_Effect = termbetaglucan-termmock,
  MOCK_AdTe = termmock-adultmock,
  MOCK_AdPre = pretermmock-adultmock,
  MOCK_PreTe = pretermmock-termmock,
  GLUC_AdTe = termbetaglucan-adultbetaglucan,
  GLUC_AdPre = pretermbetaglucan-adultbetaglucan,
  GLUC_PreTe = pretermbetaglucan-termbetaglucan,
  levels=design)

fit <- contrasts.fit(fitDupCor, cm)
```

### Moderated t-statistics

``` r
#eBayes
fitDupCor <- eBayes( fit )

#summary
summary(decideTests(fitDupCor, p.value = 0.05))
```

    ##        Adult_Effect Preterm_Effect Term_Effect MOCK_AdTe MOCK_AdPre MOCK_PreTe
    ## Down            168              0           0        33        118         13
    ## NotSig        16484          16878       16881     16817      16622      16860
    ## Up              229              3           0        31        141          8
    ##        GLUC_AdTe GLUC_AdPre GLUC_PreTe
    ## Down          34        313          0
    ## NotSig     16795      16185      16881
    ## Up            52        383          0

``` r
#extract tables
Adult_Effect = topTable(fitDupCor, coef="Adult_Effect", number = Inf)
Preterm_Effect = topTable(fitDupCor, coef="Preterm_Effect", number = Inf)
Term_Effect = topTable(fitDupCor, coef="Term_Effect", number = Inf)
MOCK_AdTe = topTable(fitDupCor, coef="MOCK_AdTe", number = Inf)
MOCK_AdPre = topTable(fitDupCor, coef="MOCK_AdPre", number = Inf)
MOCK_PreTe = topTable(fitDupCor, coef="MOCK_PreTe", number = Inf)
GLUC_AdTe = topTable(fitDupCor, coef="GLUC_AdTe", number = Inf)
GLUC_AdPre = topTable(fitDupCor, coef="GLUC_AdPre", number = Inf)
GLUC_PreTe = topTable(fitDupCor, coef="GLUC_PreTe", number = Inf)

#use cutoffs
Adult_Effect = subset(Adult_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
Preterm_Effect = subset(Preterm_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
Term_Effect = subset(Term_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
MOCK_AdTe = subset(MOCK_AdTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
MOCK_AdPre = subset(MOCK_AdPre, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
MOCK_PreTe = subset(MOCK_PreTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
GLUC_AdTe = subset(GLUC_AdTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
GLUC_AdPre= subset(GLUC_AdPre, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
GLUC_PreTe = subset(GLUC_PreTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
```

# Analysis of differential gene expression

## Figure 3A: Barchart of differentially expressed genes (up vs. downregulated)

``` r
#could all be looped!

#set back comparisons to baseline (in order to avoid wrong p-value filtering/logFC cutoffs!)
Adult_Effect = topTable(fitDupCor, coef="Adult_Effect", number = Inf)
Preterm_Effect = topTable(fitDupCor, coef="Preterm_Effect", number = Inf)
Term_Effect = topTable(fitDupCor, coef="Term_Effect", number = Inf)
MOCK_AdTe = topTable(fitDupCor, coef="MOCK_AdTe", number = Inf)
MOCK_AdPre = topTable(fitDupCor, coef="MOCK_AdPre", number = Inf)
MOCK_PreTe = topTable(fitDupCor, coef="MOCK_PreTe", number = Inf)
GLUC_AdTe = topTable(fitDupCor, coef="GLUC_AdTe", number = Inf)
GLUC_AdPre = topTable(fitDupCor, coef="GLUC_AdPre", number = Inf)
GLUC_PreTe = topTable(fitDupCor, coef="GLUC_PreTe", number = Inf)

Adult_Effect = subset(Adult_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #37
Preterm_Effect = subset(Preterm_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #3
Term_Effect = subset(Term_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #1
MOCK_AdTe = subset(MOCK_AdTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #59
MOCK_AdPre = subset(MOCK_AdPre, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #86
MOCK_PreTe = subset(MOCK_PreTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #0
GLUC_AdTe = subset(GLUC_AdTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #83
GLUC_AdPre= subset(GLUC_AdPre, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #354
GLUC_PreTe = subset(GLUC_PreTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #0

Adult_Effect_venn = Adult_Effect
Adult_Effect_venn$genes <- rownames(Adult_Effect_venn)
Adult_Effect_venn <- c(Adult_Effect_venn$genes)

Preterm_Effect_venn = Preterm_Effect
Preterm_Effect_venn$genes <- rownames(Preterm_Effect_venn)
Preterm_Effect_venn <- c(Preterm_Effect_venn$genes)

Term_Effect_venn = Term_Effect 
Term_Effect_venn$genes <- rownames(Term_Effect_venn)
Term_Effect_venn <- c(Term_Effect_venn$genes)

MOCK_AdTe_venn = MOCK_AdTe
MOCK_AdTe_venn$genes <- rownames(MOCK_AdTe_venn)
MOCK_AdTe_venn <- c(MOCK_AdTe_venn$genes)

MOCK_AdPre_venn = MOCK_AdPre
MOCK_AdPre_venn$genes <- rownames(MOCK_AdPre_venn)
MOCK_AdPre_venn <- c(MOCK_AdPre_venn$genes)

MOCK_PreTe_venn = MOCK_PreTe
MOCK_PreTe_venn$genes <- rownames(MOCK_PreTe_venn)
MOCK_PreTe_venn <- c(MOCK_PreTe_venn$genes)

GLUC_AdTe_venn = GLUC_AdTe
GLUC_AdTe_venn$genes <- rownames(GLUC_AdTe_venn)
GLUC_AdTe_venn <- c(GLUC_AdTe_venn$genes)

GLUC_AdPre_venn = GLUC_AdPre
GLUC_AdPre_venn$genes <- rownames(GLUC_AdPre_venn)
GLUC_AdPre_venn <- c(GLUC_AdPre_venn$genes)

GLUC_PreTe_venn = GLUC_PreTe
GLUC_PreTe_venn$genes <- rownames(GLUC_PreTe_venn)
GLUC_PreTe_venn <- c(GLUC_PreTe_venn$genes)

#generate a dataframe for the barchart
a <- data.frame(matrix(nrow = 1, ncol = 8)) 

#get count of all DEG included in respective  comparisons
a$Adult_Effect <- as.numeric(length(Adult_Effect_venn))
a$Preterm_Effect <- as.numeric(length(Preterm_Effect_venn))
a$Term_Effect <- as.numeric(length(Term_Effect_venn))
a$GlucAdPre <- as.numeric(length(GLUC_AdPre_venn))
a$GLUC_AdTe <- as.numeric(length(GLUC_AdTe_venn))
a$GLUC_PreTe <- as.numeric(length(GLUC_PreTe_venn))
a$MOCK_AdPre <- as.numeric(length(MOCK_AdPre_venn))
a$MOCK_AdTe <- as.numeric(length(MOCK_AdTe_venn))
a$MOCK_PreTe <- as.numeric(length(MOCK_PreTe_venn))
a <- a[, -c(1:8)]

a <- pivot_longer(a, names_to = "Group", values_to = "genes", cols = 1:9)

#filter for upregulated genes
Adult_Effect = topTable(fitDupCor, coef="Adult_Effect", number = Inf)
Preterm_Effect = topTable(fitDupCor, coef="Preterm_Effect", number = Inf)
Term_Effect = topTable(fitDupCor, coef="Term_Effect", number = Inf)
MOCK_AdTe = topTable(fitDupCor, coef="MOCK_AdTe", number = Inf)
MOCK_AdPre = topTable(fitDupCor, coef="MOCK_AdPre", number = Inf)
MOCK_PreTe = topTable(fitDupCor, coef="MOCK_PreTe", number = Inf)
GLUC_AdTe = topTable(fitDupCor, coef="GLUC_AdTe", number = Inf)
GLUC_AdPre = topTable(fitDupCor, coef="GLUC_AdPre", number = Inf)
GLUC_PreTe = topTable(fitDupCor, coef="GLUC_PreTe", number = Inf)

Adult_Effect = subset(Adult_Effect, adj.P.Val < 0.05)# & (logFC > 0.58)) #37
Preterm_Effect = subset(Preterm_Effect, adj.P.Val < 0.05)# & (logFC > 0.58)) #3
Term_Effect = subset(Term_Effect, adj.P.Val < 0.05)# & (logFC > 0.58)) #1
MOCK_AdTe = subset(MOCK_AdTe, adj.P.Val < 0.05)# & (logFC > 0.58)) #59
MOCK_AdPre = subset(MOCK_AdPre, adj.P.Val < 0.05)# & (logFC > 0.58)) #86
MOCK_PreTe = subset(MOCK_PreTe, adj.P.Val < 0.05)# & (logFC > 0.58)) #0
GLUC_AdTe = subset(GLUC_AdTe, adj.P.Val < 0.05)# & (logFC > 0.58)) #83
GLUC_AdPre= subset(GLUC_AdPre, adj.P.Val < 0.05)# & (logFC > 0.58)) #354
GLUC_PreTe = subset(GLUC_PreTe, adj.P.Val < 0.05)# & (logFC > 0.58)) #0

Adult_Effect_venn = Adult_Effect
Adult_Effect_venn$genes <- rownames(Adult_Effect_venn)
Adult_Effect_venn <- c(Adult_Effect_venn$genes)

Preterm_Effect_venn = Preterm_Effect
Preterm_Effect_venn$genes <- rownames(Preterm_Effect_venn)
Preterm_Effect_venn <- c(Preterm_Effect_venn$genes)

Term_Effect_venn = Term_Effect 
Term_Effect_venn$genes <- rownames(Term_Effect_venn)
Term_Effect_venn <- c(Term_Effect_venn$genes)

MOCK_AdTe_venn = MOCK_AdTe
MOCK_AdTe_venn$genes <- rownames(MOCK_AdTe_venn)
MOCK_AdTe_venn <- c(MOCK_AdTe_venn$genes)

MOCK_AdPre_venn = MOCK_AdPre
MOCK_AdPre_venn$genes <- rownames(MOCK_AdPre_venn)
MOCK_AdPre_venn <- c(MOCK_AdPre_venn$genes)

MOCK_PreTe_venn = MOCK_PreTe
MOCK_PreTe_venn$genes <- rownames(MOCK_PreTe_venn)
MOCK_PreTe_venn <- c(MOCK_PreTe_venn$genes)

GLUC_AdTe_venn = GLUC_AdTe
GLUC_AdTe_venn$genes <- rownames(GLUC_AdTe_venn)
GLUC_AdTe_venn <- c(GLUC_AdTe_venn$genes)

GLUC_AdPre_venn = GLUC_AdPre
GLUC_AdPre_venn$genes <- rownames(GLUC_AdPre_venn)
GLUC_AdPre_venn <- c(GLUC_AdPre_venn$genes)

GLUC_PreTe_venn = GLUC_PreTe
GLUC_PreTe_venn$genes <- rownames(GLUC_PreTe_venn)
GLUC_PreTe_venn <- c(GLUC_PreTe_venn$genes)

#generate a dataframe for the barchart
b <- data.frame(matrix(nrow = 1, ncol = 8)) 

#get count of all upregulated DEG included in respective comparisons
b$Adult_Effect <- as.numeric(length(Adult_Effect_venn))
b$Preterm_Effect <- as.numeric(length(Preterm_Effect_venn))
b$Term_Effect <- as.numeric(length(Term_Effect_venn))
b$GlucAdPre <- as.numeric(length(GLUC_AdPre_venn))
b$GLUC_AdTe <- as.numeric(length(GLUC_AdTe_venn))
b$GLUC_PreTe <- as.numeric(length(GLUC_PreTe_venn))
b$MOCK_AdPre <- as.numeric(length(MOCK_AdPre_venn))
b$MOCK_AdTe <- as.numeric(length(MOCK_AdTe_venn))
b$MOCK_PreTe <- as.numeric(length(MOCK_PreTe_venn))
b <- b[, -c(1:8)]

b <- pivot_longer(b, names_to = "Group", values_to = "genes_positive", cols = 1:9)


#filter for downregulated genes
Adult_Effect = topTable(fitDupCor, coef="Adult_Effect", number = Inf)
Preterm_Effect = topTable(fitDupCor, coef="Preterm_Effect", number = Inf)
Term_Effect = topTable(fitDupCor, coef="Term_Effect", number = Inf)
MOCK_AdTe = topTable(fitDupCor, coef="MOCK_AdTe", number = Inf)
MOCK_AdPre = topTable(fitDupCor, coef="MOCK_AdPre", number = Inf)
MOCK_PreTe = topTable(fitDupCor, coef="MOCK_PreTe", number = Inf)
GLUC_AdTe = topTable(fitDupCor, coef="GLUC_AdTe", number = Inf)
GLUC_AdPre = topTable(fitDupCor, coef="GLUC_AdPre", number = Inf)
GLUC_PreTe = topTable(fitDupCor, coef="GLUC_PreTe", number = Inf)

Adult_Effect = subset(Adult_Effect, adj.P.Val < 0.05 & (logFC < -0.58)) #37
Preterm_Effect = subset(Preterm_Effect, adj.P.Val < 0.05 & (logFC < -0.58)) #3
Term_Effect = subset(Term_Effect, adj.P.Val < 0.05 & (logFC < -0.58)) #1
MOCK_AdTe = subset(MOCK_AdTe, adj.P.Val < 0.05 & (logFC < -0.58)) #59
MOCK_AdPre = subset(MOCK_AdPre, adj.P.Val < 0.05 & (logFC < -0.58)) #86
MOCK_PreTe = subset(MOCK_PreTe, adj.P.Val < 0.05 & (logFC < -0.58)) #0
GLUC_AdTe = subset(GLUC_AdTe, adj.P.Val < 0.05 & (logFC < -0.58)) #83
GLUC_AdPre= subset(GLUC_AdPre, adj.P.Val < 0.05 & (logFC < -0.58)) #354
GLUC_PreTe = subset(GLUC_PreTe, adj.P.Val < 0.05 & (logFC < -0.58)) #0


Adult_Effect_venn = Adult_Effect
Adult_Effect_venn$genes <- rownames(Adult_Effect_venn)
Adult_Effect_venn <- c(Adult_Effect_venn$genes)

Preterm_Effect_venn = Preterm_Effect
Preterm_Effect_venn$genes <- rownames(Preterm_Effect_venn)
Preterm_Effect_venn <- c(Preterm_Effect_venn$genes)

Term_Effect_venn = Term_Effect 
Term_Effect_venn$genes <- rownames(Term_Effect_venn)
Term_Effect_venn <- c(Term_Effect_venn$genes)

MOCK_AdTe_venn = MOCK_AdTe
MOCK_AdTe_venn$genes <- rownames(MOCK_AdTe_venn)
MOCK_AdTe_venn <- c(MOCK_AdTe_venn$genes)

MOCK_AdPre_venn = MOCK_AdPre
MOCK_AdPre_venn$genes <- rownames(MOCK_AdPre_venn)
MOCK_AdPre_venn <- c(MOCK_AdPre_venn$genes)

MOCK_PreTe_venn = MOCK_PreTe
MOCK_PreTe_venn$genes <- rownames(MOCK_PreTe_venn)
MOCK_PreTe_venn <- c(MOCK_PreTe_venn$genes)

GLUC_AdTe_venn = GLUC_AdTe
GLUC_AdTe_venn$genes <- rownames(GLUC_AdTe_venn)
GLUC_AdTe_venn <- c(GLUC_AdTe_venn$genes)

GLUC_AdPre_venn = GLUC_AdPre
GLUC_AdPre_venn$genes <- rownames(GLUC_AdPre_venn)
GLUC_AdPre_venn <- c(GLUC_AdPre_venn$genes)

GLUC_PreTe_venn = GLUC_PreTe
GLUC_PreTe_venn$genes <- rownames(GLUC_PreTe_venn)
GLUC_PreTe_venn <- c(GLUC_PreTe_venn$genes)

#generate a dataframe for the barchart
c <- data.frame(matrix(nrow = 1, ncol = 8)) 

#get count of all downregulated DEG included in respective comparisons
c$Adult_Effect <- as.numeric(length(Adult_Effect_venn))
c$Preterm_Effect <- as.numeric(length(Preterm_Effect_venn))
c$Term_Effect <- as.numeric(length(Term_Effect_venn))
c$GlucAdPre <- as.numeric(length(GLUC_AdPre_venn))
c$GLUC_AdTe <- as.numeric(length(GLUC_AdTe_venn))
c$GLUC_PreTe <- as.numeric(length(GLUC_PreTe_venn))
c$MOCK_AdPre <- as.numeric(length(MOCK_AdPre_venn))
c$MOCK_AdTe <- as.numeric(length(MOCK_AdTe_venn))
c$MOCK_PreTe <- as.numeric(length(MOCK_PreTe_venn))
c <- c[, -c(1:8)]

c <- pivot_longer(c, names_to = "Group", values_to = "genes_negative", cols = 1:9)


#merge the dfs for the barchart
a <- a %>% left_join(b, by = "Group")
a <- a %>% left_join(c, by = "Group")
a <- a %>% mutate(check = a$genes_positive+a$genes_negative)
a <- a %>% mutate(pos = a$genes_positive+a$genes_negative)

a <- a %>% dplyr::select(-check, -genes_positive)

b <- a %>% pivot_longer(names_to = "condition", values_to = "value2", cols=3:4)
b <- b %>% dplyr::select(Group, value2)
a <- a %>% mutate(genes_positive = genes-genes_negative)
a <- a %>% dplyr::select(Group, genes_negative, genes_positive)
a <- a %>% pivot_longer(names_to = "condition", values_to = "Count", cols=2:3)
a$value2 <- b$value2

a$condition <- as.factor(a$condition)
a$Group <- as.factor(a$Group)

positions <- c("Adult_Effect", "Term_Effect", "Preterm_Effect", "MOCK_AdTe", "MOCK_AdPre", "MOCK_PreTe", "GLUC_AdTe","GlucAdPre","GLUC_PreTe")

#create Figure 3A
barchart_diff_exp_genes <- ggplot(a, aes(x = Group, y = Count, fill = condition)) + 
  scale_fill_manual(values = c("#7eafd1", "#de8793"), labels = c("negative", "positive")) + 
  theme_bw()  +  
  geom_col(colour = "black", width = 0.7, position = position_stack(reverse = TRUE)) + 
  geom_text(aes(y = (value2), label = Count), colour = "black", size = 5, vjust = 1.5) + scale_x_discrete(limits = positions) +
  theme(axis.text=element_text(size=8), 
        axis.title = element_text(size = 10,  face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 10))

#print Figure 3A
print(barchart_diff_exp_genes)
```

![](Trained_immunity_GIT_files/figure-gfm/setup28-1.png)<!-- -->

``` r
#export Figure 3A as svg file
graph2svg(barchart_diff_exp_genes, file = here::here("plots","barchart_diff_exp_genes"),  width = 4, height = 2.5)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/barchart_diff_exp_genes.svg

### Set back to baseline + export for supplementary information

``` r
#export for supplementary information
#create map between ENSEMBLE and ENTREZID
map <- AnnotationDbi::select(org.Hs.eg.db,
                             columns = c("ENTREZID",
                                         "ENSEMBL"),
                             keys = keys(org.Hs.eg.db, keytype = "ENTREZID")) %>%
  drop_na
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
#create annotation for features
fdata <- AnnotationDbi::select(org.Hs.eg.db,
                               columns = c("ENTREZID",
                                           "ENSEMBL",
                                           "SYMBOL",
                                           "GENENAME"),
                               keys = keys(org.Hs.eg.db, keytype = "ENTREZID")) %>%
  group_by(ENTREZID) %>%
  summarize(ENSEMBL = paste(unique(ENSEMBL), collapse = ", "),
            SYMBOL = paste(unique(SYMBOL), collapse = ", "),
            GENENAME = paste(unique(GENENAME), collapse = ", ")) %>%
  mutate(rowname = ENTREZID) %>%
  column_to_rownames
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
#eBayes
fitDupCor <- eBayes( fit )

#summary
summary(decideTests(fitDupCor, p.value = 0.05))
```

    ##        Adult_Effect Preterm_Effect Term_Effect MOCK_AdTe MOCK_AdPre MOCK_PreTe
    ## Down            168              0           0        33        118         13
    ## NotSig        16484          16878       16881     16817      16622      16860
    ## Up              229              3           0        31        141          8
    ##        GLUC_AdTe GLUC_AdPre GLUC_PreTe
    ## Down          34        313          0
    ## NotSig     16795      16185      16881
    ## Up            52        383          0

``` r
#extract tables
Adult_Effect = topTable(fitDupCor, coef="Adult_Effect", number = Inf)
Preterm_Effect = topTable(fitDupCor, coef="Preterm_Effect", number = Inf)
Term_Effect = topTable(fitDupCor, coef="Term_Effect", number = Inf)
MOCK_AdTe = topTable(fitDupCor, coef="MOCK_AdTe", number = Inf)
MOCK_AdPre = topTable(fitDupCor, coef="MOCK_AdPre", number = Inf)
MOCK_PreTe = topTable(fitDupCor, coef="MOCK_PreTe", number = Inf)
GLUC_AdTe = topTable(fitDupCor, coef="GLUC_AdTe", number = Inf)
GLUC_AdPre = topTable(fitDupCor, coef="GLUC_AdPre", number = Inf)
GLUC_PreTe = topTable(fitDupCor, coef="GLUC_PreTe", number = Inf)

#repeat to make sure
Adult_Effect = subset(Adult_Effect, adj.P.Val < 0.05)
Preterm_Effect = subset(Preterm_Effect, adj.P.Val < 0.05)
Term_Effect = subset(Term_Effect, adj.P.Val < 0.05)
MOCK_AdTe = subset(MOCK_AdTe, adj.P.Val < 0.05)
MOCK_AdPre = subset(MOCK_AdPre, adj.P.Val < 0.05)
MOCK_PreTe = subset(MOCK_PreTe, adj.P.Val < 0.05)
GLUC_AdTe = subset(GLUC_AdTe, adj.P.Val < 0.05)
GLUC_AdPre= subset(GLUC_AdPre, adj.P.Val < 0.05)
GLUC_PreTe = subset(GLUC_PreTe, adj.P.Val < 0.05)

Adult_Effect$ENSEMBL <- rownames(Adult_Effect)
Adult_Effect$ENSEMBL <- str_replace(Adult_Effect$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
Adult_Effect <- Adult_Effect %>% left_join(fdata, by = "ENSEMBL") %>% drop_na() %>% select(ENSEMBL, SYMBOL, GENENAME, logFC, AveExpr, P.Value, adj.P.Val)
write.xlsx(Adult_Effect, here::here("export_RNA","Adult_Effect_supplement_genes.xlsx"))


Preterm_Effect$ENSEMBL <- rownames(Preterm_Effect)
Preterm_Effect$ENSEMBL <- str_replace(Preterm_Effect$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
Preterm_Effect <- Preterm_Effect %>% left_join(fdata, by = "ENSEMBL") %>% drop_na() %>% select(ENSEMBL, SYMBOL, GENENAME, logFC, AveExpr, P.Value, adj.P.Val)
write.xlsx(Preterm_Effect,  here::here("export_RNA","Preterm_Effect_supplement_genes.xlsx"))

Term_Effect$ENSEMBL <- rownames(Term_Effect)
Term_Effect$ENSEMBL <- str_replace(Term_Effect$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
Term_Effect <- Term_Effect %>% left_join(fdata, by = "ENSEMBL") %>% drop_na() %>% select(ENSEMBL, SYMBOL, GENENAME, logFC, AveExpr, P.Value, adj.P.Val)
write.xlsx(Term_Effect,  here::here("export_RNA","Term_Effect_supplement_genes.xlsx"))

MOCK_AdTe$ENSEMBL <- rownames(MOCK_AdTe)
MOCK_AdTe$ENSEMBL <- str_replace(MOCK_AdTe$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
MOCK_AdTe <- MOCK_AdTe %>% left_join(fdata, by = "ENSEMBL") %>% drop_na() %>% select(ENSEMBL, SYMBOL, GENENAME, logFC, AveExpr, P.Value, adj.P.Val)
write.xlsx(MOCK_AdTe,  here::here("export_RNA","MOCK_AdTe_supplement_genes.xlsx"))

MOCK_AdPre$ENSEMBL <- rownames(MOCK_AdPre)
MOCK_AdPre$ENSEMBL <- str_replace(MOCK_AdPre$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
MOCK_AdPre <- MOCK_AdPre %>% left_join(fdata, by = "ENSEMBL") %>% drop_na() %>% select(ENSEMBL, SYMBOL, GENENAME, logFC, AveExpr, P.Value, adj.P.Val)
write.xlsx(MOCK_AdPre, here::here("export_RNA","MOCK_AdPre_supplement_genes.xlsx"))

MOCK_PreTe$ENSEMBL <- rownames(MOCK_PreTe)
MOCK_PreTe$ENSEMBL <- str_replace(MOCK_PreTe$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
MOCK_PreTe <- MOCK_PreTe %>% left_join(fdata, by = "ENSEMBL") %>% drop_na() %>% select(ENSEMBL, SYMBOL, GENENAME, logFC, AveExpr, P.Value, adj.P.Val)
write.xlsx(MOCK_PreTe, here::here("export_RNA","MOCK_PreTe_supplement_genes.xlsx"))

GLUC_AdTe$ENSEMBL <- rownames(GLUC_AdTe)
GLUC_AdTe$ENSEMBL <- str_replace(GLUC_AdTe$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
GLUC_AdTe <- GLUC_AdTe %>% left_join(fdata, by = "ENSEMBL") %>% drop_na() %>% select(ENSEMBL, SYMBOL, GENENAME, logFC, AveExpr, P.Value, adj.P.Val)
write.xlsx(GLUC_AdTe, here::here("export_RNA","GLUC_AdTe_supplement_genes.xlsx"))

GLUC_AdPre$ENSEMBL <- rownames(GLUC_AdPre)
GLUC_AdPre$ENSEMBL <- str_replace(GLUC_AdPre$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
GLUC_AdPre <- GLUC_AdPre %>% left_join(fdata, by = "ENSEMBL") %>% drop_na() %>% select(ENSEMBL, SYMBOL, GENENAME, logFC, AveExpr, P.Value, adj.P.Val)
write.xlsx(GLUC_AdPre, here::here("export_RNA","GLUC_AdPre_supplement_genes.xlsx"))

GLUC_PreTe$ENSEMBL <- rownames(GLUC_PreTe)
GLUC_PreTe$ENSEMBL <- str_replace(GLUC_PreTe$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
GLUC_PreTe <- GLUC_PreTe %>% left_join(fdata, by = "ENSEMBL") %>% drop_na() %>% select(ENSEMBL, SYMBOL, GENENAME, logFC, AveExpr, P.Value, adj.P.Val)
write.xlsx(GLUC_PreTe, here::here("export_RNA","GLUC_PreTe_supplement_genes.xlsx"))


#set everything back to normal
#eBayes
fitDupCor <- eBayes( fit )

#summary
summary(decideTests(fitDupCor, p.value = 0.05))
```

    ##        Adult_Effect Preterm_Effect Term_Effect MOCK_AdTe MOCK_AdPre MOCK_PreTe
    ## Down            168              0           0        33        118         13
    ## NotSig        16484          16878       16881     16817      16622      16860
    ## Up              229              3           0        31        141          8
    ##        GLUC_AdTe GLUC_AdPre GLUC_PreTe
    ## Down          34        313          0
    ## NotSig     16795      16185      16881
    ## Up            52        383          0

``` r
#extract tables
Adult_Effect = topTable(fitDupCor, coef="Adult_Effect", number = Inf)
Preterm_Effect = topTable(fitDupCor, coef="Preterm_Effect", number = Inf)
Term_Effect = topTable(fitDupCor, coef="Term_Effect", number = Inf)
MOCK_AdTe = topTable(fitDupCor, coef="MOCK_AdTe", number = Inf)
MOCK_AdPre = topTable(fitDupCor, coef="MOCK_AdPre", number = Inf)
MOCK_PreTe = topTable(fitDupCor, coef="MOCK_PreTe", number = Inf)
GLUC_AdTe = topTable(fitDupCor, coef="GLUC_AdTe", number = Inf)
GLUC_AdPre = topTable(fitDupCor, coef="GLUC_AdPre", number = Inf)
GLUC_PreTe = topTable(fitDupCor, coef="GLUC_PreTe", number = Inf)

#use cutoffs
Adult_Effect = subset(Adult_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
Preterm_Effect = subset(Preterm_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
Term_Effect = subset(Term_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
MOCK_AdTe = subset(MOCK_AdTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
MOCK_AdPre = subset(MOCK_AdPre, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
MOCK_PreTe = subset(MOCK_PreTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
GLUC_AdTe = subset(GLUC_AdTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
GLUC_AdPre= subset(GLUC_AdPre, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
GLUC_PreTe = subset(GLUC_PreTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))

Adult_Effect_venn = Adult_Effect
Adult_Effect_venn$genes <- rownames(Adult_Effect_venn)
Adult_Effect_venn <- c(Adult_Effect_venn$genes)

Preterm_Effect_venn = Preterm_Effect
Preterm_Effect_venn$genes <- rownames(Preterm_Effect_venn)
Preterm_Effect_venn <- c(Preterm_Effect_venn$genes)

Term_Effect_venn = Term_Effect 
Term_Effect_venn$genes <- rownames(Term_Effect_venn)
Term_Effect_venn <- c(Term_Effect_venn$genes)

MOCK_AdTe_venn = MOCK_AdTe
MOCK_AdTe_venn$genes <- rownames(MOCK_AdTe_venn)
MOCK_AdTe_venn <- c(MOCK_AdTe_venn$genes)

MOCK_AdPre_venn = MOCK_AdPre
MOCK_AdPre_venn$genes <- rownames(MOCK_AdPre_venn)
MOCK_AdPre_venn <- c(MOCK_AdPre_venn$genes)

MOCK_PreTe_venn = MOCK_PreTe
MOCK_PreTe_venn$genes <- rownames(MOCK_PreTe_venn)
MOCK_PreTe_venn <- c(MOCK_PreTe_venn$genes)

GLUC_AdTe_venn = GLUC_AdTe
GLUC_AdTe_venn$genes <- rownames(GLUC_AdTe_venn)
GLUC_AdTe_venn <- c(GLUC_AdTe_venn$genes)

GLUC_AdPre_venn = GLUC_AdPre
GLUC_AdPre_venn$genes <- rownames(GLUC_AdPre_venn)
GLUC_AdPre_venn <- c(GLUC_AdPre_venn$genes)

GLUC_PreTe_venn = GLUC_PreTe
GLUC_PreTe_venn$genes <- rownames(GLUC_PreTe_venn)
GLUC_PreTe_venn <- c(GLUC_PreTe_venn$genes)
```

## Figure 3B: PCA of MOCK-samples

``` r
#get log_cpm values
x_MOCK <- x[,x$samples$Treatment == "mock"]
x_lcpm <- cpm(x_MOCK, log = TRUE)
genes <- x_lcpm %>% as.data.frame()

genes_PCA <- t(genes) %>% as.data.frame()
genes_PCA$SampleName <- c(rownames(genes_PCA))

genes_PCA <- genes_PCA %>% left_join(metadata, by = "SampleName")

genes_PCA <- genes_PCA %>% mutate_at(c(16882:16886), as.factor)
genes_PCA <- genes_PCA %>% dplyr::select(-SampleName, -Index_2)
genes_PCA <- genes_PCA %>% mutate(index_Treat = paste0(Group, Treatment, SEP = ""))
genes_PCA$index_Treat <- as.factor(genes_PCA$index_Treat)
X <- PCA(genes_PCA [1:16881], graph = FALSE)


#PCA visualization
PCA_Mock_samples <- fviz_pca_ind(X, geom.ind = c("point"),
                                col.ind = genes_PCA$Group,
                                fill.ind = genes_PCA$Group,palette = c("#b2c5b3","#ce8793", "#8dafd1"),
                                alpha.var ="contrib",
                                addEllipses = TRUE, # Concentration ellipses,
                                ellipse.alpha = 0.4, ellipse.type   = "confidence", ellipse.level   = 0.95,
                                legend.title = "Groups", mean.point = FALSE, axes.linetype = "blank") + theme_bw() + geom_point(aes(shape = factor(genes_PCA$Group), colour = factor(genes_PCA$Group), size = 2)) + theme(axis.text.y=element_text(size=8), 
                                                                                                                                                                                                                          axis.text.x=element_text(size=8), 
                                                                                                                                                                                                                          axis.title = element_text(size = 10, face = "bold"), 
                                                                                                                                                                                                                          legend.title = element_blank(),
                                                                                                                                                                                                                          legend.text = element_text(size = 10))
#print PCA
print(PCA_Mock_samples)
```

![](Trained_immunity_GIT_files/figure-gfm/setup30-1.png)<!-- -->

``` r
#save PCA as svg file
graph2svg(PCA_Mock_samples, file = here::here("plots","PCA_Mock_samples"), width = 3.5, height = 2.5)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/PCA_Mock_samples.svg

## Figure 3C: PCA of ß-glucan samples

``` r
#get log_cpm values
x_beta <- x[,x$samples$Treatment == "betaglucan"]
x_lcpm <- cpm(x_beta, log = TRUE)
genes <- x_lcpm %>% as.data.frame()

genes_PCA <- t(genes) %>% as.data.frame()
genes_PCA$SampleName <- c(rownames(genes_PCA))

genes_PCA <- genes_PCA %>% left_join(metadata, by = "SampleName")

genes_PCA <- genes_PCA %>% mutate_at(c(16882:16886), as.factor)
genes_PCA <- genes_PCA %>% dplyr::select(-SampleName, -Index_2)
genes_PCA <- genes_PCA %>% mutate(index_Treat = paste0(Group, Treatment, SEP = ""))
genes_PCA$index_Treat <- as.factor(genes_PCA$index_Treat)
X <- PCA(genes_PCA [1:16881], graph = FALSE)


#PCA visualization
PCA_beta_samples <- fviz_pca_ind(X, geom.ind = c("point"),
                                col.ind = genes_PCA$Group,
                                fill.ind = genes_PCA$Group,palette = c("#b2c5b3","#ce8793", "#8dafd1"),
                                alpha.var ="contrib",
                                addEllipses = TRUE, # Concentration ellipses
                                ellipse.alpha = 0.4, ellipse.type   = "confidence", ellipse.level   = 0.95,
                                legend.title = "Groups", mean.point = FALSE, axes.linetype = "blank") + theme_bw() + geom_point(aes(shape = factor(genes_PCA$Group), colour = factor(genes_PCA$Group), size = 3)) + theme(axis.text.y=element_text(size=8), 
                                                                                                                                                                                                                          axis.text.x=element_text(size=8), 
                                                                                                                                                                                                                          axis.title = element_text(size = 10, face = "bold"), 
                                                                                                                                                                                                                          legend.title = element_blank(),
                                                                                                                                                                                                                          legend.text = element_text(size = 10))

#print PCA
print(PCA_beta_samples)
```

![](Trained_immunity_GIT_files/figure-gfm/setup31-1.png)<!-- -->

``` r
#save PCA as svg file
graph2svg(PCA_beta_samples, file = here::here("plots","PCA_beta_samples"),  width = 3.5, height = 2.5)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/PCA_beta_samples.svg

## Additional: Overview PCA of all samples

``` r
#get log_cpm values
x_lcpm <- cpm(x, log = TRUE)
genes <- x_lcpm %>% as.data.frame()

genes_PCA <- t(genes) %>% as.data.frame()
genes_PCA$SampleName <- c(rownames(genes_PCA))

genes_PCA <- genes_PCA %>% left_join(metadata, by = "SampleName")

genes_PCA <- genes_PCA %>% mutate_at(c(16882:16886), as.factor)
genes_PCA <- genes_PCA %>% dplyr::select(-SampleName, -Index_2)
genes_PCA <- genes_PCA %>% mutate(index_Treat = paste0(Group, Treatment, SEP = ""))
genes_PCA$index_Treat <- as.factor(genes_PCA$index_Treat)
X <- PCA(genes_PCA [1:16881], graph = FALSE)


#PCA visualization
PCA_All_samples <- fviz_pca_ind(X, geom.ind = c("point"),
                          col.ind = genes_PCA$index_Treat,
                          fill.ind = genes_PCA$index_Treat,palette = c("#b2c5b3", "#829884","#ce8793","#A95463", "#8dafd1", "#5C85AD"),
                          alpha.var ="contrib",
                          addEllipses = TRUE, # Concentration ellipses
                          ellipse.alpha = 0.4, ellipse.type = "confidence", ellipse.level   = 0.95,
                          legend.title = "Groups", mean.point = FALSE, axes.linetype = "blank") + theme_bw() + geom_point(aes(shape = factor(genes_PCA$index_Treat), colour = factor(genes_PCA$index_Treat), size = 5)) + theme(axis.text.y=element_text(size=8), 
                                                                                                                                                                  axis.text.x=element_text(size=8), 
                                                                                                                                                                  axis.title = element_text(size = 10, face = "bold"), 
                                                                                                                                                                  legend.title = element_blank(),
                                                                                                                                                                  legend.text = element_text(size = 10))

#print PCA
print(PCA_All_samples)
```

![](Trained_immunity_GIT_files/figure-gfm/setup32-1.png)<!-- -->

``` r
#save PCA as svg file
graph2svg(PCA_All_samples, file = here::here("plots", "PCA_All_samples"), width = 4, height = 2.75)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/PCA_All_samples.svg

## Figure 3D: Heatmap of all DEG

``` r
#extract vector of relevant genes for heatmap
all_diff_genes <- c(Adult_Effect_venn, GLUC_AdPre_venn, GLUC_AdTe_venn, GLUC_PreTe_venn, MOCK_AdPre_venn, MOCK_AdTe_venn, MOCK_PreTe_venn, Preterm_Effect_venn, Term_Effect_venn)
all_diff_genes <- unique(all_diff_genes)

#log_cpm transformation + data wrangling
x_lcpm <- cpm(x, log = TRUE)
genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)
diff_exp_genes <- subset(genes, (genes$ENSEMBL %in% all_diff_genes))
diff_exp_genes <- diff_exp_genes %>% dplyr::select(-ENSEMBL)
diff_exp_genes <- t(diff_exp_genes)
diff_exp_genes_scaled <- scale(diff_exp_genes) %>% (t)
diff_exp_genes_metadaten <- x$samples %>% dplyr::select(Treatment, Group)

#annotation
col_an = HeatmapAnnotation(Stimulation = diff_exp_genes_metadaten$Treatment, Group = diff_exp_genes_metadaten$Group, col = list(Group = c("adult" = "#b2c5b3", "term" = "#8dafd1", "preterm" = "#ce8793"), Stimulation = c("betaglucan" = "#4D3C7E", "mock" = "#DEC08B")))

my_palette <-  colorRampPalette(c("#0765A0", "#FFF1E4","#C61700"))(100)

#create heatmap
htmp_all_diff_genes <- Heatmap(diff_exp_genes_scaled, show_row_names = FALSE, show_row_dend = TRUE, col = my_palette,
clustering_method_columns = "complete", clustering_method_rows = "complete",
 column_dend_side = "top", column_dend_height = unit(4, "cm"), column_km = 3, column_gap =unit(3, "mm"),   column_title_gp = gpar(fontsize = 10), top_annotation = col_an)

#print heatmap
print(htmp_all_diff_genes)
```

![](Trained_immunity_GIT_files/figure-gfm/setup33-1.png)<!-- -->

``` r
#save heatmap as svg file
graph2svg(htmp_all_diff_genes, file = here::here("plots","htmp_all_diff_genes"), width = 5.5, height = 7)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/htmp_all_diff_genes.svg

## Additionally: Heatmap of DEG only using RPMI samples

``` r
#extract vector of relevant genes for heatmap
all_diff_genes <- c(Adult_Effect_venn, GLUC_AdPre_venn, GLUC_AdTe_venn, GLUC_PreTe_venn, MOCK_AdPre_venn, MOCK_AdTe_venn, MOCK_PreTe_venn, Preterm_Effect_venn, Term_Effect_venn)
all_diff_genes <- unique(all_diff_genes)

#filter for RPMI = mock samples
x_MOCK <- x[,x$samples$Treatment == "mock"]

#log_cpm transformation + data wrangling
x_lcpm <- cpm(x_MOCK, log = TRUE)
genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)
diff_exp_genes <- subset(genes, (genes$ENSEMBL %in% all_diff_genes))
diff_exp_genes <- diff_exp_genes %>% dplyr::select(-ENSEMBL)
diff_exp_genes <- t(diff_exp_genes)
diff_exp_genes_scaled <- scale(diff_exp_genes) %>% (t)
diff_exp_genes_metadaten <- x_MOCK$samples %>% dplyr::select(Group)

#annotation
col_an = HeatmapAnnotation(Group = diff_exp_genes_metadaten$Group, col = list(Group = c("adult" = "#b2c5b3", "term" = "#8dafd1", "preterm" = "#ce8793")))

my_palette <-  colorRampPalette(c("#0765A0", "#FFF1E4","#C61700"))(100)

#create heatmap
htmp_mock_diff_genes <- Heatmap(diff_exp_genes_scaled, show_row_names = FALSE, show_row_dend = TRUE, col = my_palette,
                               clustering_method_columns = "complete", clustering_method_rows = "complete",
                               column_dend_side = "top", column_dend_height = unit(4, "cm"), column_km = 3, column_gap =unit(3, "mm"),   column_title_gp = gpar(fontsize = 10), top_annotation = col_an)

#print heatmap
print(htmp_mock_diff_genes)
```

![](Trained_immunity_GIT_files/figure-gfm/setup34-1.png)<!-- -->

``` r
#save heatmap as svg file
graph2svg(htmp_mock_diff_genes, file = here::here("plots","htmp_mock_diff_genes"), width = 7, height = 8)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/htmp_mock_diff_genes.svg

## Additionally: Heatmap of DEG only using ß-glucan samples

``` r
#extract vector of relevant genes for heatmap
all_diff_genes <- c(Adult_Effect_venn, GLUC_AdPre_venn, GLUC_AdTe_venn, GLUC_PreTe_venn, MOCK_AdPre_venn, MOCK_AdTe_venn, MOCK_PreTe_venn, Preterm_Effect_venn, Term_Effect_venn)
all_diff_genes <- unique(all_diff_genes)

#filter for RPMI = mock samples
x_beta <- x[,x$samples$Treatment == "betaglucan"]

#log_cpm transformation + data wrangling
x_lcpm <- cpm(x_beta, log = TRUE)
genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)
diff_exp_genes <- subset(genes, (genes$ENSEMBL %in% all_diff_genes))
diff_exp_genes <- diff_exp_genes %>% dplyr::select(-ENSEMBL)
diff_exp_genes <- t(diff_exp_genes)
diff_exp_genes_scaled <- scale(diff_exp_genes) %>% (t)
diff_exp_genes_metadaten <- x_beta$samples %>% dplyr::select(Group)

#annotation
col_an = HeatmapAnnotation(Group = diff_exp_genes_metadaten$Group, col = list(Group = c("adult" = "#b2c5b3", "term" = "#8dafd1", "preterm" = "#ce8793")))

my_palette <-  colorRampPalette(c("#0765A0", "#FFF1E4","#C61700"))(100)

#create heatmap
htmp_beta_diff_genes <- Heatmap(diff_exp_genes_scaled, show_row_names = FALSE, show_row_dend = TRUE, col = my_palette,
                                clustering_method_columns = "complete", clustering_method_rows = "complete",
                                column_dend_side = "top", column_dend_height = unit(4, "cm"), column_km = 3, column_gap =unit(3, "mm"),   column_title_gp = gpar(fontsize = 10), top_annotation = col_an)

#print heatmap
print(htmp_beta_diff_genes)
```

![](Trained_immunity_GIT_files/figure-gfm/setup35-1.png)<!-- -->

``` r
#save heatmap as svg file
graph2svg(htmp_beta_diff_genes, file = here::here("plots","htmp_beta_diff_genes"), width = 7, height = 8)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/htmp_beta_diff_genes.svg

## Figure 3E: Overrepresentation analysis of adult samples

``` r
#set back comparisons to baseline (in order to avoid wrong p-value filtering/logFC cutoffs!)
Adult_Effect = topTable(fitDupCor, coef="Adult_Effect", number = Inf)
Preterm_Effect = topTable(fitDupCor, coef="Preterm_Effect", number = Inf)
Term_Effect = topTable(fitDupCor, coef="Term_Effect", number = Inf)
MOCK_AdTe = topTable(fitDupCor, coef="MOCK_AdTe", number = Inf)
MOCK_AdPre = topTable(fitDupCor, coef="MOCK_AdPre", number = Inf)
MOCK_PreTe = topTable(fitDupCor, coef="MOCK_PreTe", number = Inf)
GLUC_AdTe = topTable(fitDupCor, coef="GLUC_AdTe", number = Inf)
GLUC_AdPre = topTable(fitDupCor, coef="GLUC_AdPre", number = Inf)
GLUC_PreTe = topTable(fitDupCor, coef="GLUC_PreTe", number = Inf)

#get ENSEMBL for annotation later
Adult_Effect$ENSEMBL = rownames(Adult_Effect)
Preterm_Effect$ENSEMBL = rownames(Preterm_Effect)
Term_Effect$ENSEMBL = rownames(Term_Effect)
MOCK_AdTe$ENSEMBL = rownames(MOCK_AdTe)
MOCK_AdPre$ENSEMBL = rownames(MOCK_AdPre)
MOCK_PreTe$ENSEMBL = rownames(MOCK_PreTe)
GLUC_AdTe$ENSEMBL = rownames(GLUC_AdTe)
GLUC_AdPre$ENSEMBL = rownames(GLUC_AdPre)
GLUC_PreTe$ENSEMBL = rownames(GLUC_PreTe)

#loosen p.value cutoff to 0.10 as ORA has an exploratory nature, while filtering for biological relevance using logFC
DE_Adult_Effect_pos = subset(Adult_Effect, adj.P.Val < 0.1 & (logFC > 0.58)) %>% mutate(Group = "Adult_Effect") %>% mutate(subgroup = "positive") 
DE_Adult_Effect_neg = subset(Adult_Effect, adj.P.Val < 0.1 & (logFC < -0.58)) %>% mutate(Group = "Adult_Effect") %>% mutate(subgroup = "negative")
DE_Preterm_Effect_pos = subset(Preterm_Effect, adj.P.Val < 0.1 & (logFC > 0.58)) %>% mutate(Group = "Preterm_Effect") %>% mutate(subgroup = "positive")
DE_Preterm_Effect_neg = subset(Preterm_Effect, adj.P.Val < 0.1 & (logFC < -0.58)) %>% mutate(Group = "Preterm_Effect") %>% mutate(subgroup = "negative")
DE_Term_Effect_pos = subset(Term_Effect, adj.P.Val < 0.1 & (logFC > 0.58)) %>% mutate(Group = "Term_Effect") %>% mutate(subgroup = "positive")
DE_Term_Effect_neg = subset(Term_Effect, adj.P.Val < 0.1 & (logFC < -0.58)) %>% mutate(Group = "Term_Effect") %>% mutate(subgroup = "negative")

all_DE_genes <- do.call("rbind", list(DE_Adult_Effect_pos, DE_Adult_Effect_neg, DE_Preterm_Effect_pos, DE_Preterm_Effect_neg, DE_Term_Effect_pos, DE_Term_Effect_neg)) %>% as.data.frame()
all_DE_genes$ENSEMBL <- str_replace(all_DE_genes$ENSEMBL,
                                    pattern = ".[0-9]+$",
                                    replacement = "")
#annotate ENSEMBL IDs
"entrez_id" = mapIds(
  # Replace with annotation package for the organism relevant to your data
  org.Hs.eg.db,
  keys =  all_DE_genes$ENSEMBL,
  # Replace with the type of gene identifiers in your data
  keytype = "ENSEMBL",
  # Replace with the type of gene identifiers you would like to map to
  column = "ENTREZID",
  # This will keep only the first mapped value for each Ensembl ID
  multiVals = "first"
)
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
#DEG wrangling
all_DE_genes$ENTREZ <- entrez_id
all_DE_genes <- all_DE_genes %>% dplyr::filter(!is.na(ENTREZ)) %>% dplyr::select(-ENSEMBL)

#set universe (background genes) as all genes detected in RNAseq
background_genes <- readRDS(here::here("data","se.rds"))
background_genes <- background_genes@rowRanges@partitioning@NAMES %>% as.data.frame()

background_genes$ENSEMBL <- str_replace(background_genes$.,
                                    pattern = ".[0-9]+$",
                                    replacement = "")

background_genes <- background_genes %>% dplyr::select(-.)

#map background genes
"entrez_id" = mapIds(
  # Replace with annotation package for the organism relevant to your data
  org.Hs.eg.db,
  keys =  background_genes$ENSEMBL,
  # Replace with the type of gene identifiers in your data
  keytype = "ENSEMBL",
  # Replace with the type of gene identifiers you would like to map to
  column = "ENTREZID",
  # This will keep only the first mapped value for each Ensembl ID
  multiVals = "first"
)
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
#background genes wrangling
background_genes$ENTREZ <- entrez_id
background_genes <- background_genes %>% dplyr::filter(!is.na(ENTREZ)) %>% dplyr::select(-ENSEMBL)
background_genes <- background_genes$ENTREZ %>% unique()

#perform ORA using the compareCluster function and GO::terms as reference
xx.formula.twogroups <- compareCluster(ENTREZ~Group+subgroup, data=all_DE_genes,
                                       fun='enrichGO', OrgDb='org.Hs.eg.db', universe = background_genes)

#filter CompareCluster output for adult effect (only there did we observe a group-specific relevant amount of sign. expressed genes for the comp. treatment vs. RPMI)
xx.formula.twogroups <- xx.formula.twogroups %>% dplyr::filter(xx.formula.twogroups@compareClusterResult$Group == "Adult_Effect")

#plot ORA
Figure_ORA_Terms_Preterms_Adults_Beta_Glucan_effect <- enrichplot::dotplot(xx.formula.twogroups, x="GeneRatio", size = "Count", by="GeneRatio", label_format = 30, font.size = 8, showCategory = 10) + theme(axis.text.y = element_text(angle=-35, vjust = 0.95))

#print ORA
print(Figure_ORA_Terms_Preterms_Adults_Beta_Glucan_effect)
```

![](Trained_immunity_GIT_files/figure-gfm/setup36-1.png)<!-- -->

``` r
#save ORA as svg file
graph2svg(Figure_ORA_Terms_Preterms_Adults_Beta_Glucan_effect, 
          file = here::here("plots","Figure_ORA_Terms_Preterms_Adults_Beta_Glucan_effect"), 
          height = 7.35, width = 4)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/Figure_ORA_Terms_Preterms_Adults_Beta_Glucan_effect.svg

## Figure 3F: Boxplots of interesting genes relevant to ORA-annotated functions

``` r
#get feature data
fdata <- AnnotationDbi::select(org.Hs.eg.db,
                               columns = c("ENTREZID",
                                           "ENSEMBL",
                                           "SYMBOL",
                                           "GENENAME"),
                               keys = keys(org.Hs.eg.db, keytype = "ENTREZID")) %>%
  group_by(ENTREZID) %>%
  summarize(ENSEMBL = paste(unique(ENSEMBL), collapse = ", "),
            SYMBOL = paste(unique(SYMBOL), collapse = ", "),
            GENENAME = paste(unique(GENENAME), collapse = ", ")) %>%
  mutate(rowname = ENTREZID) %>%
  column_to_rownames
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
#start pulling out genes from ORA pathways
adult_genes_ORA <- xx.formula.twogroups@compareClusterResult %>% dplyr::select(Description, geneID)
adult_genes_ORA <- adult_genes_ORA[c(1:4,21:24),]

adult_ORA_up <- adult_genes_ORA[c(7,8),]
adult_ORA_down <- adult_genes_ORA[c(1:6),]

transmembrane_transporter_respiration_ORA <- adult_ORA_down$geneID %>% c() %>% unlist()
collagen_binding <- adult_ORA_up[1,] %>% dplyr::select(geneID) %>% c() %>% unlist()
calmodulin_binding <- adult_ORA_up[2,] %>% dplyr::select(geneID) %>% c() %>% unlist()

transmembrane_transporter_respiration_ORA <- strsplit(transmembrane_transporter_respiration_ORA, "/")
transmembrane_transporter_respiration_ORA <- unlist(transmembrane_transporter_respiration_ORA) %>% unique()
transmembrane_transporter_respiration_ORA  <- unique(transmembrane_transporter_respiration_ORA) %>% as.data.frame()
colnames(transmembrane_transporter_respiration_ORA) <- c("ENTREZ")


collagen_binding <- strsplit(collagen_binding, "/")
collagen_binding <- unlist(collagen_binding) %>% unique()
collagen_binding  <- unique(collagen_binding) %>% as.data.frame()
colnames(collagen_binding) <- c("ENTREZ")

calmodulin_binding <- strsplit(calmodulin_binding, "/")
calmodulin_binding <- unlist(calmodulin_binding) %>% unique()
calmodulin_binding  <- unique(calmodulin_binding) %>% as.data.frame()
colnames(calmodulin_binding) <- c("ENTREZ")

#load background genes
background_genes <- readRDS(here::here("data","se.rds"))
background_genes <- background_genes@rowRanges@partitioning@NAMES %>% as.data.frame()

background_genes$ENSEMBL <- str_replace(background_genes$.,
                                        pattern = ".[0-9]+$",
                                        replacement = "")

background_genes <- background_genes %>% dplyr::select(-.)

"entrez_id" = mapIds(
  # Replace with annotation package for the organism relevant to your data
  org.Hs.eg.db,
  keys =  background_genes$ENSEMBL,
  # Replace with the type of gene identifiers in your data
  keytype = "ENSEMBL",
  # Replace with the type of gene identifiers you would like to map to
  column = "ENTREZID",
  # This will keep only the first mapped value for each Ensembl ID
  multiVals = "first"
)
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
background_genes$ENTREZ <- entrez_id
background_genes <- background_genes %>% dplyr::filter(!is.na(ENTREZ))

transmembrane_transporter_respiration_ORA <- transmembrane_transporter_respiration_ORA %>% left_join(background_genes, by = "ENTREZ")
transmembrane_transporter_vector <- transmembrane_transporter_respiration_ORA$ENSEMBL

collagen_binding <- collagen_binding %>% left_join(background_genes, by = "ENTREZ")
collagen_vector <- collagen_binding$ENSEMBL

calmodulin_binding <- calmodulin_binding %>% left_join(background_genes, by = "ENTREZ")
calmodulin_vector <- calmodulin_binding$ENSEMBL

#get cpm-transformed expression
x_lcpm <- cpm(x)

genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)
genes$ENSEMBL <- str_replace(genes$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
rownames(genes) <- genes$ENSEMBL
genes <- t(genes) %>% as.data.frame() %>% rownames_to_column("Sample")
genes <- genes[-27,]

#wrangling
meta <- x$samples
meta <- meta %>% rownames_to_column("Sample")
meta <- meta %>% dplyr::select(-lib.size, -group, -norm.factors, -SampleName)
genes <- genes %>% left_join(meta, by = "Sample")
rownames(genes) <- genes$Sample

genes <- genes %>% pivot_longer(cols = 2:16882, values_to = "expression", names_to = "ENSEMBL")
genes$expression <- as.numeric(genes$expression)
genes$Sample <- as.factor(genes$Sample)

ORA_adult_transporters <- subset(genes, (genes$ENSEMBL %in% transmembrane_transporter_vector))
ORA_adult_transporters <- ORA_adult_transporters %>% left_join(fdata, by = "ENSEMBL")
ORA_adult_transporters$GENENAME <- as.factor(ORA_adult_transporters$GENENAME)

ORA_adult_calmodulin <- subset(genes, (genes$ENSEMBL %in% calmodulin_vector))
ORA_adult_calmodulin <- ORA_adult_calmodulin %>% left_join(fdata, by = "ENSEMBL")
ORA_adult_calmodulin$GENENAME <- as.factor(ORA_adult_calmodulin$GENENAME)

ORA_adult_collagen<- subset(genes, (genes$ENSEMBL %in% collagen_vector))
ORA_adult_collagen <- ORA_adult_collagen %>% left_join(fdata, by = "ENSEMBL")
ORA_adult_collagen$GENENAME <- as.factor(ORA_adult_collagen$GENENAME)

#bind dataframes
ORA_adults <- rbind(ORA_adult_transporters, ORA_adult_calmodulin)
ORA_adults <- rbind(ORA_adults, ORA_adult_collagen)

#plot boxplots
ORA_adults <- ggplot(ORA_adults, aes(x = new_col, y = expression, fill = Group)) + 
  scale_x_discrete(limits = c("adultmock", "adultbetaglucan", "termmock", "termbetaglucan", "pretermmock", "pretermbetaglucan")) +
  geom_boxplot()  +
  geom_line(aes(group = interaction(Group, Index), col = Group),
            alpha = 0.4) + geom_point(aes(group = interaction(Group, Index), col = Group), alpha = 0.6) + 
  scale_fill_manual(values=c("#97C591", "#ce8793", "#649CD1")) +
  scale_color_manual(values=c("#97C591", "#ce8793", "#649CD1")) +
  facet_wrap(~GENENAME, scales ="free") + ylab("AUC") + xlab("Condition") + theme_bw() +
  theme(axis.text=element_text(size=8), 
        axis.title = element_text(size = 10,  face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 10))
#print boxplots
print(ORA_adults)
```

![](Trained_immunity_GIT_files/figure-gfm/setup37-1.png)<!-- -->

``` r
#save boxplots as svg file
graph2svg(ORA_adults, 
          file = here::here("plots","ORA_adults"), 
          height = 8, width = 9)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/ORA_adults.svg

### Back to baseline

``` r
#extract tables
Adult_Effect = topTable(fitDupCor, coef="Adult_Effect", number = Inf)
Preterm_Effect = topTable(fitDupCor, coef="Preterm_Effect", number = Inf)
Term_Effect = topTable(fitDupCor, coef="Term_Effect", number = Inf)
MOCK_AdTe = topTable(fitDupCor, coef="MOCK_AdTe", number = Inf)
MOCK_AdPre = topTable(fitDupCor, coef="MOCK_AdPre", number = Inf)
MOCK_PreTe = topTable(fitDupCor, coef="MOCK_PreTe", number = Inf)
GLUC_AdTe = topTable(fitDupCor, coef="GLUC_AdTe", number = Inf)
GLUC_AdPre = topTable(fitDupCor, coef="GLUC_AdPre", number = Inf)
GLUC_PreTe = topTable(fitDupCor, coef="GLUC_PreTe", number = Inf)

#use cutoffs
Adult_Effect = subset(Adult_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
Preterm_Effect = subset(Preterm_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
Term_Effect = subset(Term_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
MOCK_AdTe = subset(MOCK_AdTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
MOCK_AdPre = subset(MOCK_AdPre, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
MOCK_PreTe = subset(MOCK_PreTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
GLUC_AdTe = subset(GLUC_AdTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
GLUC_AdPre= subset(GLUC_AdPre, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
GLUC_PreTe = subset(GLUC_PreTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))
```

## Figure 4A: Venn Diagram

``` r
#wrangling
Adult_Effect_venn = Adult_Effect
Adult_Effect_venn$genes <- rownames(Adult_Effect_venn)
Adult_Effect_venn <- c(Adult_Effect_venn$genes)

Preterm_Effect_venn = Preterm_Effect
Preterm_Effect_venn$genes <- rownames(Preterm_Effect_venn)
Preterm_Effect_venn <- c(Preterm_Effect_venn$genes)

Term_Effect_venn = Term_Effect 
Term_Effect_venn$genes <- rownames(Term_Effect_venn)
Term_Effect_venn <- c(Term_Effect_venn$genes)

MOCK_AdTe_venn = MOCK_AdTe
MOCK_AdTe_venn$genes <- rownames(MOCK_AdTe_venn)
MOCK_AdTe_venn <- c(MOCK_AdTe_venn$genes)

MOCK_AdPre_venn = MOCK_AdPre
MOCK_AdPre_venn$genes <- rownames(MOCK_AdPre_venn)
MOCK_AdPre_venn <- c(MOCK_AdPre_venn$genes)

MOCK_PreTe_venn = MOCK_PreTe
MOCK_PreTe_venn$genes <- rownames(MOCK_PreTe_venn)
MOCK_PreTe_venn <- c(MOCK_PreTe_venn$genes)

GLUC_AdTe_venn = GLUC_AdTe
GLUC_AdTe_venn$genes <- rownames(GLUC_AdTe_venn)
GLUC_AdTe_venn <- c(GLUC_AdTe_venn$genes)

GLUC_AdPre_venn = GLUC_AdPre
GLUC_AdPre_venn$genes <- rownames(GLUC_AdPre_venn)
GLUC_AdPre_venn <- c(GLUC_AdPre_venn$genes)

GLUC_PreTe_venn = GLUC_PreTe
GLUC_PreTe_venn$genes <- rownames(GLUC_PreTe_venn)
GLUC_PreTe_venn <- c(GLUC_PreTe_venn$genes)


list <- list(MOCK_AdTe = MOCK_AdTe_venn,
             MOCK_AdPre = MOCK_AdPre_venn,
             GLUC_AdTe = GLUC_AdTe_venn,
             GLUC_AdPre = GLUC_AdPre_venn)

             
venn.diagram(list, filename = here::here("plots","venn_diagram_LogFc2.svg"), resolution = 1200, height = 2.4, width = 2.1, imagetype = "svg",
             category.names = c("M-adult vs. M-term" , "M-adult vs. M-preterm" , "BG-adult vs. BG-term", "BG-adult vs. BG-preterm"),
             # Circles
             lwd = 1.5, col = "#BFBFBF",fontfamily ="sans serif",
             #lty = 'blank',
             fill = c("#7dafd1", "#ce8793", "#8fadc1", "#be8793"), alpha = 0.4, 
             # Numbers
             cex = 0.9,
             #fontface = "bold",
             # Set names
             cat.cex = 0.7,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(0.22, 0.22, 0.11, 0.11))             
```

    ## [1] 1

## Analysis of relevant intersects

``` r
#get intersects
intersect <- calculate.overlap(list)

#manually annotate intersects (could also be looped!)

intersect_analysis_ALL_conserved <- intersect$a6 %>% as.data.frame()
intersect_analysis_ALL_conserved$ENSEMBL <- str_replace(intersect_analysis_ALL_conserved$., pattern = ".[0-9]+$", replacement = "")
colnames(intersect_analysis_ALL_conserved) <- c("x", "ENSEMBL")
intersect_analysis_ALL_conserved <- intersect_analysis_ALL_conserved %>% dplyr::select(-x)
intersect_analysis_ALL_conserved$ENTREZID <- mapIds(org.Hs.eg.db, keys = intersect_analysis_ALL_conserved$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
intersect_analysis_ALL_conserved <- intersect_analysis_ALL_conserved %>% drop_na(ENTREZID)
intersect_analysis_ALL_conserved <- intersect_analysis_ALL_conserved %>% dplyr::select(-ENSEMBL)

intersect_analysis_BETA_conserved <- intersect$a2 %>% as.data.frame()
intersect_analysis_BETA_conserved$ENSEMBL <- str_replace(intersect_analysis_BETA_conserved$., pattern = ".[0-9]+$", replacement = "")
colnames(intersect_analysis_BETA_conserved) <- c("x", "ENSEMBL")
intersect_analysis_BETA_conserved <- intersect_analysis_BETA_conserved %>% dplyr::select(-x)
intersect_analysis_BETA_conserved$ENTREZID <- mapIds(org.Hs.eg.db, keys = intersect_analysis_BETA_conserved$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
intersect_analysis_BETA_conserved <- intersect_analysis_BETA_conserved %>% drop_na(ENTREZID)
intersect_analysis_BETA_conserved <- intersect_analysis_BETA_conserved %>% dplyr::select(-ENSEMBL)

intersect_analysis_BETA_preterm <- intersect$a3 %>% as.data.frame()
intersect_analysis_BETA_preterm$ENSEMBL <- str_replace(intersect_analysis_BETA_preterm$., pattern = ".[0-9]+$", replacement = "")
colnames(intersect_analysis_BETA_preterm) <- c("x", "ENSEMBL")
intersect_analysis_BETA_preterm <- intersect_analysis_BETA_preterm %>% dplyr::select(-x)
intersect_analysis_BETA_preterm$ENTREZID <- mapIds(org.Hs.eg.db, keys = intersect_analysis_BETA_preterm$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
intersect_analysis_BETA_preterm <- intersect_analysis_BETA_preterm %>% drop_na(ENTREZID)
intersect_analysis_BETA_preterm <- intersect_analysis_BETA_preterm %>% dplyr::select(-ENSEMBL)

intersect_analysis_BETA_term <- intersect$a1 %>% as.data.frame()
intersect_analysis_BETA_term$ENSEMBL <- str_replace(intersect_analysis_BETA_term$., pattern = ".[0-9]+$", replacement = "")
colnames(intersect_analysis_BETA_term) <- c("x", "ENSEMBL")
intersect_analysis_BETA_term <- intersect_analysis_BETA_term %>% dplyr::select(-x)
intersect_analysis_BETA_term$ENTREZID <- mapIds(org.Hs.eg.db, keys = intersect_analysis_BETA_term$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
intersect_analysis_BETA_term <- intersect_analysis_BETA_term %>% drop_na(ENTREZID)
intersect_analysis_BETA_term <- intersect_analysis_BETA_term %>% dplyr::select(-ENSEMBL)

intersect_analysis_MOCK_term <- intersect$a9 %>% as.data.frame()
intersect_analysis_MOCK_term$ENSEMBL <- str_replace(intersect_analysis_MOCK_term$., pattern = ".[0-9]+$", replacement = "")
colnames(intersect_analysis_MOCK_term) <- c("x", "ENSEMBL")
intersect_analysis_MOCK_term <- intersect_analysis_MOCK_term %>% dplyr::select(-x)
intersect_analysis_MOCK_term$ENTREZID <- mapIds(org.Hs.eg.db, keys = intersect_analysis_MOCK_term$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
intersect_analysis_MOCK_term <- intersect_analysis_MOCK_term %>% drop_na(ENTREZID)
intersect_analysis_MOCK_term <- intersect_analysis_MOCK_term %>% dplyr::select(-ENSEMBL)

intersect_analysis_MOCK_preterm <- intersect$a14 %>% as.data.frame()
intersect_analysis_MOCK_preterm$ENSEMBL <- str_replace(intersect_analysis_MOCK_preterm$., pattern = ".[0-9]+$", replacement = "")
colnames(intersect_analysis_MOCK_preterm) <- c("x", "ENSEMBL")
intersect_analysis_MOCK_preterm <- intersect_analysis_MOCK_preterm %>% dplyr::select(-x)
intersect_analysis_MOCK_preterm$ENTREZID <- mapIds(org.Hs.eg.db, keys = intersect_analysis_MOCK_preterm$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
intersect_analysis_MOCK_preterm <- intersect_analysis_MOCK_preterm %>% drop_na(ENTREZID)
intersect_analysis_MOCK_preterm <- intersect_analysis_MOCK_preterm %>% dplyr::select(-ENSEMBL)

intersect_analysis_preterm_conserved <- intersect$a8 %>% as.data.frame()
intersect_analysis_preterm_conserved$ENSEMBL <- str_replace(intersect_analysis_preterm_conserved$., pattern = ".[0-9]+$", replacement = "")
colnames(intersect_analysis_preterm_conserved) <- c("x", "ENSEMBL")
intersect_analysis_preterm_conserved <- intersect_analysis_preterm_conserved %>% dplyr::select(-x)
intersect_analysis_preterm_conserved$ENTREZID <- mapIds(org.Hs.eg.db, keys = intersect_analysis_preterm_conserved$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
intersect_analysis_preterm_conserved <- intersect_analysis_preterm_conserved %>% drop_na(ENTREZID)
intersect_analysis_preterm_conserved <- intersect_analysis_preterm_conserved %>% dplyr::select(-ENSEMBL)

intersect_analysis_term_conserved <- intersect$a4 %>% as.data.frame()
intersect_analysis_term_conserved$ENSEMBL <- str_replace(intersect_analysis_term_conserved$., pattern = ".[0-9]+$", replacement = "")
colnames(intersect_analysis_term_conserved) <- c("x", "ENSEMBL")
intersect_analysis_term_conserved <- intersect_analysis_term_conserved %>% dplyr::select(-x)
intersect_analysis_term_conserved$ENTREZID <- mapIds(org.Hs.eg.db, keys = intersect_analysis_term_conserved$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
intersect_analysis_term_conserved <- intersect_analysis_term_conserved %>% drop_na(ENTREZID)
intersect_analysis_term_conserved <- intersect_analysis_term_conserved %>% dplyr::select(-ENSEMBL)
```

### Get SYMBOL and GENENAME and export intersects for manual annotation though literature search

``` r
#get annotations
map <- AnnotationDbi::select(org.Hs.eg.db,
                             columns = c("ENTREZID",
                                         "ENSEMBL"),
                             keys = keys(org.Hs.eg.db, keytype = "ENTREZID")) %>%
  drop_na
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
fdata <- AnnotationDbi::select(org.Hs.eg.db,
                               columns = c("ENTREZID",
                                           "ENSEMBL",
                                           "SYMBOL",
                                           "GENENAME"),
                               keys = keys(org.Hs.eg.db, keytype = "ENTREZID")) %>%
  group_by(ENTREZID) %>%
  summarize(ENSEMBL = paste(unique(ENSEMBL), collapse = ", "),
            SYMBOL = paste(unique(SYMBOL), collapse = ", "),
            GENENAME = paste(unique(GENENAME), collapse = ", ")) %>%
  mutate(rowname = ENTREZID) %>%
  column_to_rownames
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
#main ones
intersect_analysis_ALL_conserved <- intersect_analysis_ALL_conserved %>% left_join(fdata, by = "ENTREZID")
intersect_analysis_BETA_conserved <- intersect_analysis_BETA_conserved %>% left_join(fdata, by = "ENTREZID")
write.xlsx(intersect_analysis_ALL_conserved, here::here("export_RNA","intersect_analysis_ALL_conserved_old.xlsx")) #double-checked, correct
write.xlsx(intersect_analysis_BETA_conserved, here::here("export_RNA", "intersect_analysis_BETA_conserved_old.xlsx")) #double-checked, correct

#preterm
intersect_analysis_BETA_preterm <- intersect_analysis_BETA_preterm  %>% left_join(fdata, by = "ENTREZID") 
intersect_analysis_MOCK_preterm <- intersect_analysis_MOCK_preterm  %>% left_join(fdata, by = "ENTREZID")
intersect_analysis_preterm_conserved <- intersect_analysis_preterm_conserved  %>% left_join(fdata, by = "ENTREZID")
write.xlsx(intersect_analysis_BETA_preterm, here::here("export_RNA","intersect_analysis_BETA_preterm_old.xlsx")) #double-checked, correct
write.xlsx(intersect_analysis_MOCK_preterm, here::here("export_RNA","intersect_analysis_MOCK_preterm_old.xlsx")) #double-checked, correct
write.xlsx(intersect_analysis_preterm_conserved, here::here("export_RNA","intersect_analysis_preterm_conserved_old.xlsx")) #double-checked, correct

#term
intersect_analysis_BETA_term <- intersect_analysis_BETA_term  %>% left_join(fdata, by = "ENTREZID")
intersect_analysis_MOCK_term <- intersect_analysis_MOCK_term  %>% left_join(fdata, by = "ENTREZID")
intersect_analysis_term_conserved <- intersect_analysis_term_conserved  %>% left_join(fdata, by = "ENTREZID")
write.xlsx(intersect_analysis_BETA_term, here::here("export_RNA","intersect_analysis_BETA_term_old.xlsx")) #double-checked, correct
write.xlsx(intersect_analysis_MOCK_term, here::here("export_RNA","intersect_analysis_MOCK_term_old.xlsx")) #double-checked, correct
write.xlsx(intersect_analysis_term_conserved, here::here("export_RNA","intersect_analysis_term_conserved_old.xlsx")) #double-checked, correct

#relevant intersects were manually analysed through literature search
```

## Figure 5B: Conserved genes in neonates, only sign. upon treatment with ß-Glucan

``` r
x_lcpm <- cpm(x)

genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)
genes$ENSEMBL <- str_replace(genes$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
rownames(genes) <- genes$ENSEMBL

vector_beta_conserved <- rio::import(here::here("data", "ENTREZID_beta_conserved.xlsx"))
vector_beta_conserved <- vector_beta_conserved$ENTREZID

genes_new <- genes
genes_new$ENTREZID <- mapIds(org.Hs.eg.db, keys = genes_new$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
genes_new <- genes_new %>% drop_na(ENTREZID)

genes_of_interest <- genes_new %>% subset(ENTREZID %in% vector_beta_conserved)
genes_of_interest <- genes_of_interest %>% left_join(fdata, by = "ENTREZID")
genes_of_interest <- genes_of_interest %>% dplyr::select(-ENSEMBL.x, -ENSEMBL.y)
rownames(genes_of_interest) <- genes_of_interest$SYMBOL
genes_of_interest <- t(genes_of_interest) %>% as.data.frame()
genes_of_interest <- genes_of_interest %>% rownames_to_column("SampleName")
genes_of_interest <- genes_of_interest[-c(27:29),]
genes_of_interest <- genes_of_interest %>% left_join(metadata, by = "SampleName")
rownames(genes_of_interest) <- genes_of_interest$SampleName 
genes_of_interest <- genes_of_interest%>% dplyr::select(-SampleName)
genes_of_interest <- genes_of_interest %>% mutate_at(c(1:13), as.numeric)
genes_of_interest <- genes_of_interest %>% mutate_at(c(14:17), as.factor)

genes_of_interest <- genes_of_interest %>% pivot_longer(cols = 1:13, names_to = "gene", values_to = "expression")
genes_of_interest <- genes_of_interest %>% mutate(new_col = paste0(Group, Treatment, SEP = ""))
genes_of_interest$new_col <- as.factor(genes_of_interest$new_col)
genes_of_interest_betaglucan <- genes_of_interest %>% dplyr::filter(Treatment == "betaglucan")

genes_of_interest_beta_conserved <- ggplot(genes_of_interest_betaglucan, aes(x=Group, y=expression, fill = Group))+ geom_boxplot() + scale_x_discrete(limits = c("adult", "term", "preterm")) +
  geom_point(aes(x = Group, y = expression, col = Group), alpha = 0.6) + scale_fill_manual(values=c("#97C591", "#ce8793", "#649CD1")) + scale_color_manual(values=c("#97C591", "#ce8793", "#649CD1"))+
  facet_wrap(gene ~., scales='free') + ylab("Expression") + xlab("Condition") + theme_bw() +
  theme(axis.text=element_text(size=7), 
        axis.title = element_text(size = 8,  face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 8))

print(genes_of_interest_beta_conserved)
```

![](Trained_immunity_GIT_files/figure-gfm/setup42-1.png)<!-- -->

``` r
graph2svg(genes_of_interest_beta_conserved, file = here::here("plots","genes_of_interest_beta_conserved"), width = 5.45, height = 5.18)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/genes_of_interest_beta_conserved.svg

## Figure 5B: Conserved genes in neonates, significant irrespective of treatment

``` r
#most interesting four genes identified through manual analysis
vector_genes_conserved <- c("10643", "115265", "8528", "10642")

genes_new <- genes
genes_new$ENTREZID <- mapIds(org.Hs.eg.db, keys = genes_new$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
genes_new <- genes_new %>% drop_na(ENTREZID)

genes_of_interest <- genes_new %>% subset(ENTREZID %in% vector_genes_conserved)
genes_of_interest <- genes_of_interest %>% left_join(fdata, by = "ENTREZID")
genes_of_interest <- genes_of_interest %>% dplyr::select(-ENSEMBL.x, -ENSEMBL.y)
rownames(genes_of_interest) <- genes_of_interest$SYMBOL
genes_of_interest <- t(genes_of_interest) %>% as.data.frame()
genes_of_interest <- genes_of_interest %>% rownames_to_column("SampleName")
genes_of_interest <- genes_of_interest[-c(27:29),]
genes_of_interest <- genes_of_interest %>% left_join(metadata, by = "SampleName")
rownames(genes_of_interest) <- genes_of_interest$SampleName 
genes_of_interest <- genes_of_interest%>% dplyr::select(-SampleName)
genes_of_interest <- genes_of_interest %>% mutate_at(c(1:4), as.numeric)
genes_of_interest <- genes_of_interest %>% mutate_at(c(5:8), as.factor)

genes_of_interest <- genes_of_interest %>% pivot_longer(cols = 1:4, names_to = "gene", values_to = "expression")
genes_of_interest <- genes_of_interest %>% mutate(new_col = paste0(Group, Treatment, SEP = ""))
genes_of_interest$new_col <- as.factor(genes_of_interest$new_col)


genes_of_interest_conserved_ALL_no_ind_of_treatment <- ggplot(genes_of_interest, aes(x=new_col, y=expression, fill = Group)) + geom_boxplot()+ 
  scale_x_discrete(limits = c("adultmock", "adultbetaglucan", "termmock", "termbetaglucan", "pretermmock", "pretermbetaglucan")) + 
  geom_point(aes(x = new_col, y = expression, col = Group), alpha = 0.6) + 
  scale_fill_manual(values=c("#97C591", "#ce8793", "#649CD1")) + scale_color_manual(values=c("#97C591", "#ce8793", "#649CD1")) +
  facet_wrap(gene ~., scales='free') + theme_bw()

print(genes_of_interest_conserved_ALL_no_ind_of_treatment)
```

![](Trained_immunity_GIT_files/figure-gfm/setup43-1.png)<!-- -->

``` r
graph2svg(genes_of_interest_conserved_ALL_no_ind_of_treatment, file = here::here("plots","genes_of_interest_conserved_ALL_no_ind_of_treatment"), width = 3.38, height = 2.78) 
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/genes_of_interest_conserved_ALL_no_ind_of_treatment.svg

## Figure 5D: Genes sign. different between term neonates and adults upon treatment with ß-Glucan

``` r
x_lcpm <- cpm(x)

genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)
genes$ENSEMBL <- str_replace(genes$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
rownames(genes) <- genes$ENSEMBL

vector_beta_conserved_terms <- rio::import(here::here("data","ENTREZID_beta_term.xlsx"))
vector_beta_conserved_terms <- vector_beta_conserved_terms$ENTREZID

genes_new <- genes
genes_new$ENTREZID <- mapIds(org.Hs.eg.db, keys = genes_new$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
genes_new <- genes_new %>% drop_na(ENTREZID)

genes_of_interest <- genes_new %>% subset(ENTREZID %in% vector_beta_conserved_terms)
genes_of_interest <- genes_of_interest %>% left_join(fdata, by = "ENTREZID")
genes_of_interest <- genes_of_interest %>% dplyr::select(-ENSEMBL.x, -ENSEMBL.y)
rownames(genes_of_interest) <- genes_of_interest$SYMBOL
genes_of_interest <- t(genes_of_interest) %>% as.data.frame()
genes_of_interest <- genes_of_interest %>% rownames_to_column("SampleName")
genes_of_interest <- genes_of_interest[-c(27:29),]
genes_of_interest <- genes_of_interest %>% left_join(metadata, by = "SampleName")
rownames(genes_of_interest) <- genes_of_interest$SampleName 
genes_of_interest <- genes_of_interest %>% dplyr::select(-SampleName)
genes_of_interest <- genes_of_interest %>% mutate_at(c(1:7), as.numeric)
genes_of_interest <- genes_of_interest %>% mutate_at(c(8:11), as.factor)

genes_of_interest <- genes_of_interest %>% pivot_longer(cols = 1:7, names_to = "gene", values_to = "expression")
genes_of_interest <- genes_of_interest %>% mutate(new_col = paste0(Group, Treatment, SEP = ""))
genes_of_interest$new_col <- as.factor(genes_of_interest$new_col)
genes_of_interest_betaglucan_terms <- genes_of_interest %>% dplyr::filter(Treatment == "betaglucan")

genes_of_interest_betaglucan_terms <- ggplot(genes_of_interest_betaglucan_terms, aes(x=Group, y=expression, fill = Group))+ geom_boxplot() + scale_x_discrete(limits = c("adult", "term", "preterm")) +
  geom_point(aes(x = Group, y = expression, col = Group), alpha = 0.6) + scale_fill_manual(values=c("#97C591", "#ce8793", "#649CD1")) + scale_color_manual(values=c("#97C591", "#ce8793", "#649CD1"))+
  facet_wrap(gene ~., scales='free') + ylab("Expression") + xlab("Condition") + theme_bw() +
  theme(axis.text=element_text(size=7), 
        axis.title = element_text(size = 8,  face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 8))

print(genes_of_interest_betaglucan_terms)
```

![](Trained_immunity_GIT_files/figure-gfm/setup44-1.png)<!-- -->

``` r
graph2svg(genes_of_interest_betaglucan_terms, file = here::here("plots","genes_of_interest_betaglucan_terms"), width = 4.33, height = 3.95)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/genes_of_interest_betaglucan_terms.svg

## Figure 5E: Genes sign. different between preterm neonates and adults upon treatment with ß-Glucan

``` r
x_lcpm <- cpm(x)

genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)
genes$ENSEMBL <- str_replace(genes$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
rownames(genes) <- genes$ENSEMBL

vector_beta_conserved_preterms <- rio::import(here::here("data","ENTREZID_manually_curated_Beta_Preterm.xlsx"))
vector_beta_conserved_preterms <- vector_beta_conserved_preterms$ENTREZID

genes_new <- genes
genes_new$ENTREZID <- mapIds(org.Hs.eg.db, keys = genes_new$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
genes_new <- genes_new %>% drop_na(ENTREZID)

genes_of_interest <- genes_new %>% subset(ENTREZID %in% vector_beta_conserved_preterms)
genes_of_interest <- genes_of_interest %>% left_join(fdata, by = "ENTREZID")
genes_of_interest <- genes_of_interest %>% dplyr::select(-ENSEMBL.x, -ENSEMBL.y)
rownames(genes_of_interest) <- genes_of_interest$SYMBOL
genes_of_interest <- t(genes_of_interest) %>% as.data.frame()
genes_of_interest <- genes_of_interest %>% rownames_to_column("SampleName")
genes_of_interest <- genes_of_interest[-c(27:29),]
genes_of_interest <- genes_of_interest %>% left_join(metadata, by = "SampleName")
rownames(genes_of_interest) <- genes_of_interest$SampleName 
genes_of_interest <- genes_of_interest%>% dplyr::select(-SampleName)
genes_of_interest <- genes_of_interest %>% mutate_at(c(1:38), as.numeric)
genes_of_interest <- genes_of_interest %>% mutate_at(c(39:42), as.factor)

genes_of_interest <- genes_of_interest %>% pivot_longer(cols = 1:38, names_to = "gene", values_to = "expression")
genes_of_interest <- genes_of_interest %>% mutate(new_col = paste0(Group, Treatment, SEP = ""))
genes_of_interest$new_col <- as.factor(genes_of_interest$new_col)
genes_of_interest_betaglucan <- genes_of_interest %>% dplyr::filter(Treatment == "betaglucan")

genes_of_interest_2 <- genes_of_interest %>% dplyr::filter(gene == "LDHD" | gene == "KDM5B" | gene == "KDM6B" | gene == "HDAC11")


genes_of_interest_Treatment_Group_2 <- ggplot(genes_of_interest_2, aes(x=Group, y=expression, fill = Group))+ geom_boxplot() + scale_x_discrete(limits = c("adult", "term", "preterm")) +
  geom_point(aes(x = Group, y = expression, col = Group), alpha = 0.6) + scale_fill_manual(values=c("#97C591", "#ce8793", "#649CD1")) + scale_color_manual(values=c("#97C591", "#ce8793", "#649CD1"))+
  facet_wrap(gene ~., scales='free') + ylab("Expression") + xlab("Condition") + theme_bw() +
  theme(axis.text=element_text(size=7), 
        axis.title = element_text(size = 8,  face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 8))

print(genes_of_interest_Treatment_Group_2)
```

![](Trained_immunity_GIT_files/figure-gfm/setup45-1.png)<!-- -->

``` r
graph2svg(genes_of_interest_Treatment_Group_2, file = here::here("plots","genes_of_interest_Treatment_Group_2"), width = 3.353, height = 2.695)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/genes_of_interest_Treatment_Group_2.svg

## Figure 5G: Exploratory overview between term and preterm infants (nothing sign.)

``` r
x_lcpm <- cpm(x)

genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)
genes$ENSEMBL <- str_replace(genes$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
rownames(genes) <- genes$ENSEMBL

vector_expl_beta_term_preterm <- c("3111", "81793")

genes_new <- genes
genes_new$ENTREZID <- mapIds(org.Hs.eg.db, keys = genes_new$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
genes_new <- genes_new %>% drop_na(ENTREZID)

genes_of_interest <- genes_new %>% subset(ENTREZID %in% vector_expl_beta_term_preterm)
genes_of_interest <- genes_of_interest %>% left_join(fdata, by = "ENTREZID")
genes_of_interest <- genes_of_interest %>% dplyr::select(-ENSEMBL.x, -ENSEMBL.y)
rownames(genes_of_interest) <- genes_of_interest$SYMBOL
genes_of_interest <- t(genes_of_interest) %>% as.data.frame()
genes_of_interest <- genes_of_interest %>% rownames_to_column("SampleName")
genes_of_interest <- genes_of_interest[-c(27:29),]
genes_of_interest <- genes_of_interest %>% left_join(metadata, by = "SampleName")
rownames(genes_of_interest) <- genes_of_interest$SampleName 
genes_of_interest <- genes_of_interest%>% dplyr::select(-SampleName)
genes_of_interest <- genes_of_interest %>% mutate_at(c(1:2), as.numeric)
genes_of_interest <- genes_of_interest %>% mutate_at(c(3:6), as.factor)

genes_of_interest <- genes_of_interest %>% pivot_longer(cols = 1:2, names_to = "gene", values_to = "expression")
genes_of_interest <- genes_of_interest %>% mutate(new_col = paste0(Group, Treatment, SEP = ""))
genes_of_interest$new_col <- as.factor(genes_of_interest$new_col)
genes_of_interest_betaglucan <- genes_of_interest %>% dplyr::filter(Treatment == "betaglucan")

#remove one outlier
genes_of_interest_betaglucan <- genes_of_interest_betaglucan[-26,]

expl_beta_term_preterm <- ggplot(genes_of_interest_betaglucan, aes(x=Group, y=expression, fill = Group))+ geom_boxplot(outlier.shape = NA) + scale_x_discrete(limits = c("adult", "term", "preterm")) +
  geom_point(aes(x = Group, y = expression, col = Group), alpha = 0.6) + scale_fill_manual(values=c("#97C591", "#ce8793", "#649CD1")) + scale_color_manual(values=c("#97C591", "#ce8793", "#649CD1"))+
  facet_wrap(gene ~., scales='free') + ylab("Expression") + xlab("Condition") + theme_bw() +
  theme(axis.text=element_text(size=7), 
        axis.title = element_text(size = 8,  face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 8))

graph2svg(expl_beta_term_preterm, file = here::here("plots","expl_beta_term_preterm"), width = 3.3, height = 1.46)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/expl_beta_term_preterm.svg

## Figure 5 preparation: GSEA analysis loop

``` r
#return to initial results without filters to perform GSEA (needs ranked list of all genes)
Adult_Effect = topTable(fitDupCor, coef="Adult_Effect", number = Inf)
Preterm_Effect = topTable(fitDupCor, coef="Preterm_Effect", number = Inf)
Term_Effect = topTable(fitDupCor, coef="Term_Effect", number = Inf)
MOCK_AdTe = topTable(fitDupCor, coef="MOCK_AdTe", number = Inf)
MOCK_AdPre = topTable(fitDupCor, coef="MOCK_AdPre", number = Inf)
MOCK_PreTe = topTable(fitDupCor, coef="MOCK_PreTe", number = Inf)
GLUC_AdTe = topTable(fitDupCor, coef="GLUC_AdTe", number = Inf)
GLUC_AdPre = topTable(fitDupCor, coef="GLUC_AdPre", number = Inf)
GLUC_PreTe = topTable(fitDupCor, coef="GLUC_PreTe", number = Inf)


list_DF <- list(Adult_Effect = Adult_Effect, Preterm_Effect= Preterm_Effect, 
                Term_Effect = Term_Effect, MOCK_AdTe = MOCK_AdTe, MOCK_AdPre = MOCK_AdPre, 
                MOCK_PreTe = MOCK_PreTe, GLUC_AdTe = GLUC_AdTe, GLUC_AdPre = GLUC_AdPre, 
                GLUC_PreTe = GLUC_PreTe)


gene_list_MF <- list()
allgenes_GO_GSEA_MF <- list()
for (i in 1:9) {
  input.df <- list_DF[[i]]  #use LFC shrinkage df
  #extract the ids for GSEA - for ENSEMBL
  all_genes_ranked <- input.df$logFC
  #select LFC
  names(all_genes_ranked) <- row.names(input.df)
  all_genes_ranked <- sort(all_genes_ranked, decreasing = TRUE) #rank genes in descending order according to LFC
  
  all_genes_ranked <- data.frame(all_genes_ranked)
  all_genes_ranked$ENSEMBL <- rownames(all_genes_ranked)
  all_genes_ranked$ENSEMBL <- str_replace(all_genes_ranked$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
  rownames(all_genes_ranked) <- all_genes_ranked$ENSEMBL
  colnames(all_genes_ranked) <- c("logFC", "ENSEMBL")
  
  geneList = as.numeric(all_genes_ranked[,1])
  
  ## feature 2: named vector
  names(geneList) = as.character(all_genes_ranked[,2])
  
  ## feature 3: decreasing order
  geneList1 = sort(geneList, decreasing = TRUE)
  gene_list_MF[[i]] <- geneList1
  ##save to list
  GSEA_MF <- gseGO(geneList = geneList1, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENSEMBL", 
                      ont = "MF",
                      pvalueCutoff = 0.05,
                      minGSSize = 10,
                      maxGSSize = 500)
  allgenes_GO_GSEA_MF[[i]] <- GSEA_MF
}
```

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

``` r
for (i in c(1:9)) {
  dotGSEA <- enrichplot::dotplot(allgenes_GO_GSEA_MF[[i]], showCategory=20)
  print(dotGSEA)
}
```

![](Trained_immunity_GIT_files/figure-gfm/setup47-1.png)<!-- -->![](Trained_immunity_GIT_files/figure-gfm/setup47-2.png)<!-- -->![](Trained_immunity_GIT_files/figure-gfm/setup47-3.png)<!-- -->![](Trained_immunity_GIT_files/figure-gfm/setup47-4.png)<!-- -->![](Trained_immunity_GIT_files/figure-gfm/setup47-5.png)<!-- -->![](Trained_immunity_GIT_files/figure-gfm/setup47-6.png)<!-- -->![](Trained_immunity_GIT_files/figure-gfm/setup47-7.png)<!-- -->![](Trained_immunity_GIT_files/figure-gfm/setup47-8.png)<!-- -->![](Trained_immunity_GIT_files/figure-gfm/setup47-9.png)<!-- -->

## Figure 5C: Dotplot ß-Glucan term vs. adult

``` r
#Glucan term vs. adult
GSEA_MF_GO_Glucan_term_adult_x <- enrichplot::dotplot(allgenes_GO_GSEA_MF[[7]], showCategory=100, title = "Beta Term_Adult") + xlim(0,0.75)
GSEA_MF_GO_Glucan_term_adult_vec <- GSEA_MF_GO_Glucan_term_adult_x$data
GSEA_MF_GO_Glucan_term_adult_vec <- GSEA_MF_GO_Glucan_term_adult_vec[order(GSEA_MF_GO_Glucan_term_adult_vec$GeneRatio, decreasing = TRUE),]
GSEA_MF_GO_Glucan_term_adult_vec <- head(GSEA_MF_GO_Glucan_term_adult_vec, 5)
GSEA_MF_GO_Glucan_term_adult_vec <- c(GSEA_MF_GO_Glucan_term_adult_vec$Description)

GSEA_MF_GO_Glucan_term_adult_x <- enrichplot::dotplot(allgenes_GO_GSEA_MF[[7]], showCategory=GSEA_MF_GO_Glucan_term_adult_vec, title = "Beta Term_Adult") + 
  xlim(0,0.75)

print(GSEA_MF_GO_Glucan_term_adult_x)
```

![](Trained_immunity_GIT_files/figure-gfm/setup48-1.png)<!-- -->

``` r
graph2svg(GSEA_MF_GO_Glucan_term_adult_x, 
          file = here::here("plots","GSEA Beta Term_Adult MF"), 
          height = 1.72, width = 4.3)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/GSEA Beta Term_Adult MF.svg

## Figure 5E: Dotplot ß-Glucan preterm vs. adult

``` r
GSEA_MF_GO_Glucan_preterm_adult_x <- enrichplot::dotplot(allgenes_GO_GSEA_MF[[8]], showCategory=100, title = "Beta Preterm_Adult") + xlim(0,0.75)
GSEA_MF_GO_Glucan_preterm_adult_vec <- GSEA_MF_GO_Glucan_preterm_adult_x$data
GSEA_MF_GO_Glucan_preterm_adult_vec <- GSEA_MF_GO_Glucan_preterm_adult_vec[order(GSEA_MF_GO_Glucan_preterm_adult_vec$GeneRatio, decreasing = TRUE),]
GSEA_MF_GO_Glucan_preterm_adult_vec <- head(GSEA_MF_GO_Glucan_preterm_adult_vec, 5)
GSEA_MF_GO_Glucan_preterm_adult_vec <- c(GSEA_MF_GO_Glucan_preterm_adult_vec$Description)

GSEA_MF_GO_Glucan_preterm_adult_x <- enrichplot::dotplot(allgenes_GO_GSEA_MF[[8]], showCategory=GSEA_MF_GO_Glucan_preterm_adult_vec, title = "Beta Preterm_Adult") + 
  xlim(0,0.75)

print(GSEA_MF_GO_Glucan_preterm_adult_x)
```

![](Trained_immunity_GIT_files/figure-gfm/setup49-1.png)<!-- -->

``` r
graph2svg(GSEA_MF_GO_Glucan_preterm_adult_x, 
          file = here::here("plots","GSEA Beta Preterm_Adult MF"), 
          height = 1.7, width = 4.08)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/GSEA Beta Preterm_Adult MF.svg

## Figure 5G: Dotplot ß-Glucan preterm vs. term

``` r
GSEA_MF_GO_Glucan_term_preterm_x <- enrichplot::dotplot(allgenes_GO_GSEA_MF[[9]], showCategory=100, title = "Beta Term_Preterm") + xlim(0,0.75)
GSEA_MF_GO_Glucan_term_preterm_vec <- GSEA_MF_GO_Glucan_term_preterm_x$data
GSEA_MF_GO_Glucan_term_preterm_vec <- GSEA_MF_GO_Glucan_term_preterm_vec[order(GSEA_MF_GO_Glucan_term_preterm_vec$GeneRatio, decreasing = TRUE),]
GSEA_MF_GO_Glucan_term_preterm_vec <- head(GSEA_MF_GO_Glucan_term_preterm_vec, 5)
GSEA_MF_GO_Glucan_term_preterm_vec <- c(GSEA_MF_GO_Glucan_term_preterm_vec$Description)

GSEA_MF_GO_Glucan_term_preterm_x <- enrichplot::dotplot(allgenes_GO_GSEA_MF[[9]], showCategory=GSEA_MF_GO_Glucan_term_preterm_vec, title = "Beta Term_Preterm") + 
  xlim(0,0.75)

print(GSEA_MF_GO_Glucan_term_preterm_x)
```

![](Trained_immunity_GIT_files/figure-gfm/setup50-1.png)<!-- -->

``` r
graph2svg(GSEA_MF_GO_Glucan_term_preterm_x, 
          file = here::here("plots","GSEA Beta Term_Preterm MF"), 
          height = 1.70, width = 4.123)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/GSEA Beta Term_Preterm MF.svg

### Back to baseline

``` r
Adult_Effect = topTable(fitDupCor, coef="Adult_Effect", number = Inf)
Preterm_Effect = topTable(fitDupCor, coef="Preterm_Effect", number = Inf)
Term_Effect = topTable(fitDupCor, coef="Term_Effect", number = Inf)
MOCK_AdTe = topTable(fitDupCor, coef="MOCK_AdTe", number = Inf)
MOCK_AdPre = topTable(fitDupCor, coef="MOCK_AdPre", number = Inf)
MOCK_PreTe = topTable(fitDupCor, coef="MOCK_PreTe", number = Inf)
GLUC_AdTe = topTable(fitDupCor, coef="GLUC_AdTe", number = Inf)
GLUC_AdPre = topTable(fitDupCor, coef="GLUC_AdPre", number = Inf)
GLUC_PreTe = topTable(fitDupCor, coef="GLUC_PreTe", number = Inf)

Adult_Effect = subset(Adult_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) %>% rownames_to_column("ENSEMBL") #37
Adult_Effect$ENSEMBL <- str_replace(Adult_Effect$ENSEMBL, pattern = ".[0-9]+$", replacement = "")

Preterm_Effect = subset(Preterm_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) %>% rownames_to_column("ENSEMBL") #3
Preterm_Effect$ENSEMBL <- str_replace(Preterm_Effect$ENSEMBL, pattern = ".[0-9]+$", replacement = "")

Term_Effect = subset(Term_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) %>% rownames_to_column("ENSEMBL") #1
Term_Effect$ENSEMBL <- str_replace(Term_Effect$ENSEMBL, pattern = ".[0-9]+$", replacement = "")

MOCK_AdTe = subset(MOCK_AdTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) %>% rownames_to_column("ENSEMBL")#59
MOCK_AdTe$ENSEMBL <- str_replace(MOCK_AdTe$ENSEMBL, pattern = ".[0-9]+$", replacement = "")

MOCK_AdPre = subset(MOCK_AdPre, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) %>% rownames_to_column("ENSEMBL")#86
MOCK_AdPre$ENSEMBL <- str_replace(MOCK_AdPre$ENSEMBL, pattern = ".[0-9]+$", replacement = "")

MOCK_PreTe = subset(MOCK_PreTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) %>% rownames_to_column("ENSEMBL")#0
MOCK_PreTe$ENSEMBL <- str_replace(MOCK_PreTe$ENSEMBL, pattern = ".[0-9]+$", replacement = "")

GLUC_AdTe = subset(GLUC_AdTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) %>% rownames_to_column("ENSEMBL")#83
GLUC_AdTe$ENSEMBL <- str_replace(GLUC_AdTe$ENSEMBL, pattern = ".[0-9]+$", replacement = "")

GLUC_AdPre= subset(GLUC_AdPre, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) %>% rownames_to_column("ENSEMBL")#354
GLUC_AdPre$ENSEMBL <- str_replace(GLUC_AdPre$ENSEMBL, pattern = ".[0-9]+$", replacement = "")

GLUC_PreTe = subset(GLUC_PreTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) %>% rownames_to_column("ENSEMBL")#0
GLUC_PreTe$ENSEMBL <- str_replace(GLUC_PreTe$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
```

## Figure 5C: Heatplot ß-Glucan term vs. adult

``` r
GSEA_MF_GO_Glucan_term_adult <- enrichplot::dotplot(allgenes_GO_GSEA_MF[[7]], showCategory=100, title = "Beta Term_Adult")
graph2svg(GSEA_MF_GO_Glucan_term_adult, 
          file = here::here("plots","GSEA Beta Term_Adult MF"), 
          height = 4, width = 4.35)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/GSEA Beta Term_Adult MF.svg

``` r
GSEA_MF_GO_Glucan_term_adult_data <- GSEA_MF_GO_Glucan_term_adult[["data"]] %>% as.data.frame()
GSEA_MF_GO_Glucan_term_adult_data <- GSEA_MF_GO_Glucan_term_adult_data %>% dplyr::filter(Description =="organic acid binding" | 
                                                                                         Description == "molecular carrier activity" |
                                                                                         Description == "N6-methyladenosine-containing RNA binding" | 
                                                                                         Description == "mRNA 5'-UTR binding"|
                                                                                         Description == "oxygen binding")

vector_Glucan_term_adult <- GSEA_MF_GO_Glucan_term_adult_data$core_enrichment
vector_Glucan_term_adult <- strsplit(vector_Glucan_term_adult, "/")[[1]]
vector_Glucan_term_adult <- vector_Glucan_term_adult %>% unique()

Glucan_AdTe_GSEA <- subset(GLUC_AdTe,(GLUC_AdTe$ENSEMBL %in% vector_Glucan_term_adult))
#0
edox <- setReadable(allgenes_GO_GSEA_MF[[7]], 'org.Hs.eg.db', 'ENSEMBL')
edox <- edox %>% dplyr::filter(edox@result$Description == "organic acid binding" |
                          edox@result$Description == "molecular carrier activity" | 
                          edox@result$Description == "N6-methyladenosine-containing RNA binding" | 
                          edox@result$Description == "mRNA 5'-UTR binding" | 
                          edox@result$Description == "oxygen binding")

gene_list_logFC <-  gene_list_MF[[7]]
gene_list_filtered_pos <- gene_list_logFC[gene_list_logFC >= 1] 
gene_list_filtered_neg <- gene_list_logFC[gene_list_logFC <= -1] 
gene_list_logFC <- c(gene_list_filtered_pos, gene_list_filtered_neg)

Heatplot_GSEA_GLUCAN_AdTe <- heatplot(edox,  showCategory=5, foldChange = gene_list_logFC) + scale_y_discrete(limits = c("oxygen binding", 
                                                                                                                         "mRNA 5'-UTR binding",
                                                                                                                         "N6-methyladenosine-containing RNA binding",
                                                                                                                         "molecular carrier activity",
                                                                                                                         "organic acid binding")) + scale_fill_gradientn(name = "fold change",
                                                                                                                                                                                          colours = c("#0765A0","#337eac","#5d95b8","#86acc3","#aec3ce","#d6dad9","#FFF1E4",
                                                                                                                                                                                                      "#FFF1E4","#f7d2c4","#efb3a3","#e69382","#de725f","#d44c37","#c61700"), 
                                                                                                                                                                                          values = rescale(c(-4,-3.5,-3,-2.5,-2,-1.5,-1,1,1.5,2,2.5,3,3.5, 4)),
                                                                                                                                                                                          limits=c(-4, 4),
                                                                                                                                                                                          breaks=c(-4 ,-2, 0,2, 4)) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

``` r
graph2svg(Heatplot_GSEA_GLUCAN_AdTe, 
          file = here::here("plots","Heatplot_GSEA_GLUCAN_AdTe"), 
          height = 1.70, width = 8.7)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/Heatplot_GSEA_GLUCAN_AdTe.svg

## Figure 5E: Heatplot ß-Glucan preterm vs. adult

``` r
GSEA_MF_GO_Glucan_preterm_adult <- enrichplot::dotplot(allgenes_GO_GSEA_MF[[8]], showCategory=100, title = "Beta Preterm_Adult")
graph2svg(GSEA_MF_GO_Glucan_preterm_adult, 
          file = here::here("plots","GSEA Beta Preterm_Adult MF"), 
          height = 4, width = 4.40)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/GSEA Beta Preterm_Adult MF.svg

``` r
GSEA_MF_GO_Glucan_preterm_adult_data <- GSEA_MF_GO_Glucan_preterm_adult[["data"]] %>% as.data.frame()
GSEA_MF_GO_Glucan_preterm_adult_data <- GSEA_MF_GO_Glucan_preterm_adult_data %>% dplyr::filter(Description == "single-stranded DNA helicase activity" | 
                                                                                           Description =="MHC class II protein complex binding" | 
                                                                                           Description == "peptide antigen binding" | 
                                                                                           Description == "structural constituent of ribosome"| 
                                                                                           Description == "MHC protein complex binding")

vector_Glucan_preterm_adult <- GSEA_MF_GO_Glucan_preterm_adult_data$core_enrichment
vector_Glucan_preterm_adult <- strsplit(vector_Glucan_preterm_adult, "/")[[1]]
vector_Glucan_preterm_adult <- vector_Glucan_preterm_adult %>% unique()

Glucan_AdPre_GSEA <- subset(GLUC_AdPre,(GLUC_AdPre$ENSEMBL %in% vector_Glucan_preterm_adult))
#0
edox <- setReadable(allgenes_GO_GSEA_MF[[8]], 'org.Hs.eg.db', 'ENSEMBL')
edox <- edox %>% dplyr::filter(edox@result$Description == "single-stranded DNA helicase activity" |
                        edox@result$Description =="MHC class II protein complex binding" | 
                        edox@result$Description == "peptide antigen binding"  | 
                        edox@result$Description == "structural constituent of ribosome"| 
                        edox@result$Description == "MHC protein complex binding")

gene_list_logFC <-  gene_list_MF[[8]]
gene_list_filtered_pos <- gene_list_logFC[gene_list_logFC >= 1] 
gene_list_filtered_neg <- gene_list_logFC[gene_list_logFC <= -1] 
gene_list_logFC <- c(gene_list_filtered_pos, gene_list_filtered_neg)

Heatplot_GSEA_GLUCAN_AdPre <- heatplot(edox,  showCategory=5, foldChange =gene_list_logFC) + scale_y_discrete(limits = c("MHC protein complex binding", 
                                                                                                                          "structural constituent of ribosome",
                                                                                                                          "peptide antigen binding",
                                                                                                                          "MHC class II protein complex binding",
                                                                                                                          "single-stranded DNA helicase activity")) + scale_fill_gradientn(name = "fold change",
                                                                                    colours = c("#0765A0","#337eac","#5d95b8","#86acc3","#aec3ce","#d6dad9","#FFF1E4",
                                                                                                "#FFF1E4","#f7d2c4","#efb3a3","#e69382","#de725f","#d44c37","#c61700"), 
                                                                                    values = rescale(c(-4,-3.5,-3,-2.5,-2,-1.5,-1,1,1.5,2,2.5,3,3.5, 4)),
                                                                                    limits=c(-4, 4),
                                                                                    breaks=c(-4 ,-2, 0,2, 4)) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

``` r
graph2svg(Heatplot_GSEA_GLUCAN_AdPre, 
          file = here::here("plots","Heatplot_GSEA_GLUCAN_AdPre"), 
          height = 1.56, width = 17)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/Heatplot_GSEA_GLUCAN_AdPre.svg

## Figure 5G: Heatplot ß-Glucan preterm vs. term

``` r
GSEA_MF_GO_Glucan_term_preterm <- enrichplot::dotplot(allgenes_GO_GSEA_MF[[9]], showCategory=20, title = "Beta Term_Preterm")
graph2svg(GSEA_MF_GO_Glucan_term_preterm, 
          file = here::here("plots","GSEA Beta Term_Preterm MF"), 
          height = 8.5, width = 6.3)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/GSEA Beta Term_Preterm MF.svg

``` r
GSEA_MF_GO_Glucan_term_preterm_data <- GSEA_MF_GO_Glucan_term_preterm[["data"]] %>% as.data.frame()
GSEA_MF_GO_Glucan_term_preterm_data <- GSEA_MF_GO_Glucan_term_preterm_data %>% dplyr::filter(Description == "MHC class II protein complex binding" |
                                                                                                 Description == "MHC protein complex binding" |
                                                                                                 Description =="immune receptor activity" | 
                                                                                                 Description == "peptide antigen binding" |
                                                                                                 Description == "hydrolase activity, acting on glycosyl bonds")

vector_Glucan_term_preterm <- GSEA_MF_GO_Glucan_term_preterm_data$core_enrichment
vector_Glucan_term_preterm <- strsplit(vector_Glucan_term_preterm, "/")[[1]]
vector_Glucan_term_preterm <- vector_Glucan_term_preterm %>% unique()

Glucan_PreTe_GSEA <- subset(GLUC_PreTe,(GLUC_PreTe$ENSEMBL %in% vector_Glucan_term_preterm))
#0
edox <- setReadable(allgenes_GO_GSEA_MF[[9]], 'org.Hs.eg.db', 'ENSEMBL')
edox <- edox %>% dplyr::filter(edox@result$Description == "MHC class II protein complex binding" |
                          edox@result$Description =="MHC protein complex binding" | 
                          edox@result$Description == "immune receptor activity" | 
                          edox@result$Description == "peptide antigen binding"| 
                          edox@result$Description == "hydrolase activity, acting on glycosyl bonds")

gene_list_logFC <-  gene_list_MF[[9]]
gene_list_filtered_pos <- gene_list_logFC[gene_list_logFC >= 1] 
gene_list_filtered_neg <- gene_list_logFC[gene_list_logFC <= -1] 
gene_list_logFC <- c(gene_list_filtered_pos, gene_list_filtered_neg)

Heatplot_GSEA_GLUCAN_PreTe <- heatplot(edox,  showCategory=5, foldChange =gene_list_logFC) +  scale_y_discrete(limits = c("hydrolase activity, acting on glycosyl bonds", 
                                                                                                                          "peptide antigen binding","immune receptor activity",
                                                                                                                          "MHC protein complex binding",
                                                                                                                          "MHC class II protein complex binding")) + scale_fill_gradientn(name = "fold change",
                                                                                    colours = c("#0765A0","#337eac","#5d95b8","#86acc3","#aec3ce","#d6dad9","#FFF1E4",
                                                                                                "#FFF1E4","#f7d2c4","#efb3a3","#e69382","#de725f","#d44c37","#c61700"), 
                                                                                    values = rescale(c(-4,-3.5,-3,-2.5,-2,-1.5,-1,1,1.5,2,2.5,3,3.5, 4)),
                                                                                    limits=c(-4, 4),
                                                                                    breaks=c(-4 ,-2, 0,2, 4)) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

``` r
graph2svg(Heatplot_GSEA_GLUCAN_PreTe, 
          file = here::here("plots","Heatplot_GSEA_GLUCAN_PreTe"), 
          height = 1.5, width = 13.8)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/Heatplot_GSEA_GLUCAN_PreTe.svg

## lnc-RNA analysis

``` r
#get ensemble annotations
hub <- AnnotationHub()
query(hub, c("homo sapiens","ensdb"))
```

    ## AnnotationHub with 26 records
    ## # snapshotDate(): 2023-10-23
    ## # $dataprovider: Ensembl
    ## # $species: Homo sapiens
    ## # $rdataclass: EnsDb
    ## # additional mcols(): taxonomyid, genome, description,
    ## #   coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
    ## #   rdatapath, sourceurl, sourcetype 
    ## # retrieve records with, e.g., 'object[["AH53211"]]' 
    ## 
    ##              title                             
    ##   AH53211  | Ensembl 87 EnsDb for Homo Sapiens 
    ##   AH53715  | Ensembl 88 EnsDb for Homo Sapiens 
    ##   AH56681  | Ensembl 89 EnsDb for Homo Sapiens 
    ##   AH57757  | Ensembl 90 EnsDb for Homo Sapiens 
    ##   AH60773  | Ensembl 91 EnsDb for Homo Sapiens 
    ##   ...        ...                               
    ##   AH104864 | Ensembl 107 EnsDb for Homo sapiens
    ##   AH109336 | Ensembl 108 EnsDb for Homo sapiens
    ##   AH109606 | Ensembl 109 EnsDb for Homo sapiens
    ##   AH113665 | Ensembl 110 EnsDb for Homo sapiens
    ##   AH116291 | Ensembl 111 EnsDb for Homo sapiens

``` r
ensdb <- hub[["AH113665"]]
```

    ## loading from cache

``` r
require("ensembldb")
gns <- genes(ensdb)

#retrieve differentially expressed genes from all comparisons
all_diff_genes <- c(Adult_Effect_venn, GLUC_AdPre_venn, GLUC_AdTe_venn, GLUC_PreTe_venn, MOCK_AdPre_venn, MOCK_AdTe_venn, MOCK_PreTe_venn, Preterm_Effect_venn, Term_Effect_venn)
all_diff_genes <- unique(all_diff_genes)

#wrangling
x_lcpm <- cpm(x, log = TRUE)
genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)
diff_exp_genes <- subset(genes, (genes$ENSEMBL %in% all_diff_genes))
diff_exp_genes$ENSEMBL <- str_replace(diff_exp_genes$ENSEMBL,
                             pattern = ".[0-9]+$",
                             replacement = "")
#get lncRNAs
haha <- c(gns$gene_biotype)
haha <- unique(haha)
lncs <- gns[gns$gene_biotype %in% "lncRNA"]
xy <- c(lncs@ranges@NAMES)
xz <- unique(xy)

#subset genes
lncs_genes <- subset(diff_exp_genes, (diff_exp_genes$ENSEMBL %in% xz))
lncs_genes <- lncs_genes %>% dplyr::select(-ENSEMBL)
lncs_genes_t <- t(lncs_genes)
```

## PCA Figure 5: Overview on all lncRNAs

``` r
x_lcpm <- cpm(x, log = TRUE)
genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)

genes$ENSEMBL <- str_replace(genes$ENSEMBL,
                        pattern = ".[0-9]+$",
                        replacement = "")

haha <- c(gns$gene_biotype)
haha <- unique(haha)
lncs <- gns[gns$gene_biotype %in% "lncRNA"]
xy <- c(lncs@ranges@NAMES)
xz <- unique(xy)

lncs_genes <- subset(genes, (genes$ENSEMBL %in% xz))
lncs_genes <- lncs_genes %>% dplyr::select(-ENSEMBL)
lncs_genes <- t(lncs_genes) %>% as.data.frame()

lncs_genes$SampleName <- c(rownames(lncs_genes))
lncs_genes <- lncs_genes %>% left_join(metadata, by = "SampleName")

lncs_genes <- lncs_genes %>% mutate_at(c(2284:2288), as.factor)
lncs_genes <- lncs_genes %>% dplyr::select(-SampleName, -Index_2)
lncs_genes <- lncs_genes %>% mutate(index_Treat = paste0(Group, Treatment, SEP = ""))
lncs_genes$index_Treat <- as.factor(lncs_genes$index_Treat)
X <- PCA(lncs_genes [1:2283], graph = FALSE)



PCA_All_lnc_RNA_samples <- fviz_pca_ind(X, geom.ind = c("point"),
                                col.ind = lncs_genes$index_Treat,
                                fill.ind = lncs_genes$index_Treat,palette = c("#b2c5b3", "#829884","#ce8793","#A95463", "#8dafd1", "#5C85AD"),
                                alpha.var ="contrib",
                                addEllipses = TRUE, # Concentration ellipses
                                ellipse.alpha = 0.4, ellipse.type   = "confidence", ellipse.level   = 0.95,
                                legend.title = "Groups", mean.point = FALSE, axes.linetype = "blank") + theme_bw() + geom_point(aes(shape = factor(lncs_genes$index_Treat), colour = factor(lncs_genes$index_Treat), size = 5)) +
                                theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=8), axis.title = element_text(size = 10, face = "bold"), legend.title = element_blank(), legend.text = element_text(size = 10))

print(PCA_All_lnc_RNA_samples)
```

![](Trained_immunity_GIT_files/figure-gfm/setup56-1.png)<!-- -->

``` r
graph2svg(PCA_All_lnc_RNA_samples, file = here::here("plots","PCA_All_lnc_RNA_samples"), width = 3.5, height = 2.5)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/PCA_All_lnc_RNA_samples.svg

## Figure 5I: PCA MOCK

``` r
#PCA MOCK lncRNA
x_MOCK <- x[,x$samples$Treatment == "mock"]
x_lcpm <- cpm(x_MOCK, log = TRUE)
genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)

genes$ENSEMBL <- str_replace(genes$ENSEMBL,
                             pattern = ".[0-9]+$",
                             replacement = "")

haha <- c(gns$gene_biotype)
haha <- unique(haha)
lncs <- gns[gns$gene_biotype %in% "lncRNA"]
xy <- c(lncs@ranges@NAMES)
xz <- unique(xy)

lncs_genes <- subset(genes, (genes$ENSEMBL %in% xz))
lncs_genes <- lncs_genes %>% dplyr::select(-ENSEMBL)
lncs_genes <- t(lncs_genes) %>% as.data.frame()

lncs_genes$SampleName <- c(rownames(lncs_genes))
lncs_genes <- lncs_genes %>% left_join(metadata, by = "SampleName")

lncs_genes <- lncs_genes %>% mutate_at(c(2284:2288), as.factor)
lncs_genes <- lncs_genes %>% dplyr::select(-SampleName, -Index_2)
lncs_genes <- lncs_genes %>% mutate(index_Treat = paste0(Group, Treatment, SEP = ""))
lncs_genes$index_Treat <- as.factor(lncs_genes$index_Treat)

#remove one sample that has been shown to be an lnc_RNA outlier in MOCK PCA that has previously been performed (not in here)
lncs_genes <- lncs_genes[-13,]
X <- PCA(lncs_genes [1:2283], graph = FALSE)


PCA_Mock_lncRNA_samples <- fviz_pca_ind(X, geom.ind = c("point"),
                                        col.ind = lncs_genes$Group,
                                        fill.ind = lncs_genes$Group,palette = c("#b2c5b3","#ce8793", "#8dafd1"),
                                        alpha.var ="contrib",
                                        select.ind = lncs_genes$SampleName,
                                        select.car = "name",
                                        addEllipses = TRUE, # Concentration ellipses,
                                        ellipse.alpha = 0.4, ellipse.type   = "confidence", ellipse.level   = 0.95,
                                        legend.title = "Groups", mean.point = FALSE, axes.linetype = "blank") + theme_bw() + 
                                        geom_point(aes(shape = factor(lncs_genes$Group), colour = factor(lncs_genes$Group), size = 2)) + theme(axis.text.y=element_text(size=8),
                                                                                                                                               axis.text.x=element_text(size=8),
                                                                                                                                               axis.title = element_text(size = 10, face = "bold"),
                                                                                                                                               legend.title = element_blank(),
                                                                                                                                               legend.text = element_text(size = 10))

print(PCA_Mock_lncRNA_samples)
```

![](Trained_immunity_GIT_files/figure-gfm/setup57-1.png)<!-- -->

``` r
graph2svg(PCA_Mock_lncRNA_samples, file = here::here("plots","PCA_Mock_lncRNA_samples"), width = 3, height = 1.95)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/PCA_Mock_lncRNA_samples.svg

## Figure 5I: PCA ß-Glucan

``` r
x_beta <- x[,x$samples$Treatment == "betaglucan"]
x_lcpm <- cpm(x_beta, log = TRUE)
genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)

genes$ENSEMBL <- str_replace(genes$ENSEMBL,
                             pattern = ".[0-9]+$",
                             replacement = "")

haha <- c(gns$gene_biotype)
haha <- unique(haha)
lncs <- gns[gns$gene_biotype %in% "lncRNA"]
xy <- c(lncs@ranges@NAMES)
xz <- unique(xy)

lncs_genes <- subset(genes, (genes$ENSEMBL %in% xz))
lncs_genes <- lncs_genes %>% dplyr::select(-ENSEMBL)
lncs_genes <- t(lncs_genes) %>% as.data.frame()

lncs_genes$SampleName <- c(rownames(lncs_genes))
lncs_genes <- lncs_genes %>% left_join(metadata, by = "SampleName")

lncs_genes <- lncs_genes %>% mutate_at(c(2284:2288), as.factor)
lncs_genes <- lncs_genes %>% dplyr::select(-SampleName, -Index_2)
lncs_genes <- lncs_genes %>% mutate(index_Treat = paste0(Group, Treatment, SEP = ""))
lncs_genes$index_Treat <- as.factor(lncs_genes$index_Treat)

#remove one sample that has been shown to be an lnc_RNA outlier in MOCK PCA
lncs_genes <- lncs_genes[-13,]
X <- PCA(lncs_genes [1:2283], graph = FALSE)


PCA_Beta_lncRNA_samples <- fviz_pca_ind(X, geom.ind = c("point"),
                                 col.ind = lncs_genes$Group,
                                 fill.ind = lncs_genes$Group,palette = c("#b2c5b3","#ce8793", "#8dafd1"),
                                 alpha.var ="contrib",
                                 addEllipses = TRUE, # Concentration ellipses,
                                 ellipse.alpha = 0.4, ellipse.type  = "confidence", ellipse.level   = 0.95,
                                 legend.title = "Groups", mean.point = FALSE, axes.linetype = "blank") + theme_bw() + geom_point(aes(shape = factor(lncs_genes$Group), colour = factor(lncs_genes$Group), size = 2)) +
                                                                                                                                                                              theme(axis.text.y=element_text(size=8),
                                                                                                                                                                                    axis.text.x=element_text(size=8),
                                                                                                                                                                                    axis.title = element_text(size = 10, face = "bold"),
                                                                                                                                                                                    legend.title = element_blank(),
                                                                                                                                                                                    legend.text = element_text(size = 10))
print(PCA_Beta_lncRNA_samples)
```

![](Trained_immunity_GIT_files/figure-gfm/setup58-1.png)<!-- -->

``` r
graph2svg(PCA_Beta_lncRNA_samples, file = here::here("plots","PCA_Beta_lncRNA_samples"), width = 3, height = 1.95)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/PCA_Beta_lncRNA_samples.svg

## Figure 5J: lncRNA heatmap

``` r
lncs_genes_scaled <- scale(lncs_genes_t) %>% (t)
lncs_genes_metadaten <- x$samples %>% dplyr::select(Treatment, Group)

col_an = HeatmapAnnotation(Stimulation = lncs_genes_metadaten$Treatment, Group = lncs_genes_metadaten$Group, col = list(Group = c(adult = "#b2c5b3", "term" = "#8dafd1", "preterm" = "#ce8793"), Stimulation = c("betaglucan" = "#4D3C7E", "mock" = "#DEC08B")))

my_palette <-  colorRampPalette(c("#0765A0", "#FFF1E4","#C61700"))(100)

htmp_long_coding <- Heatmap(lncs_genes_scaled, show_row_names = FALSE, show_row_dend = TRUE, col = my_palette, show_column_names = FALSE, 
                               clustering_method_columns = "complete", clustering_method_rows = "complete",
                               column_dend_side = "top", column_dend_height = unit(0.75, "cm"), column_km = 3, column_gap =unit(3, "mm"),   column_title_gp = gpar(fontsize = 10), top_annotation = col_an)

print(htmp_long_coding)
```

![](Trained_immunity_GIT_files/figure-gfm/setup59-1.png)<!-- -->

``` r
graph2svg(htmp_long_coding, file = here::here("plots","htmp_long_coding"), width = 4.65, height = 3.85)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/htmp_long_coding.svg

### Data preparation for pathway inferences using lncPath

``` r
all_diff_genes <- c(Adult_Effect_venn, GLUC_AdPre_venn, GLUC_AdTe_venn, GLUC_PreTe_venn, MOCK_AdPre_venn, MOCK_AdTe_venn, MOCK_PreTe_venn, Preterm_Effect_venn, Term_Effect_venn)
all_diff_genes <- unique(all_diff_genes)

x_lcpm <- cpm(x, log = TRUE)

genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)

diff_exp_genes <- subset(genes, (genes$ENSEMBL %in% all_diff_genes))
diff_exp_genes$ENSEMBL <- str_replace(diff_exp_genes$ENSEMBL,
                                      pattern = ".[0-9]+$",
                                      replacement = "")
haha <- c(gns$gene_biotype)
haha <- unique(haha)
lncs <- gns[gns$gene_biotype %in% "lncRNA"]
xy <- c(lncs@ranges@NAMES)
xz <- unique(xy)


lncs_genes <- subset(diff_exp_genes, (diff_exp_genes$ENSEMBL %in% xz))

#back to basic
Adult_Effect = topTable(fitDupCor, coef="Adult_Effect", number = Inf)
Preterm_Effect = topTable(fitDupCor, coef="Preterm_Effect", number = Inf)
Term_Effect = topTable(fitDupCor, coef="Term_Effect", number = Inf)
MOCK_AdTe = topTable(fitDupCor, coef="MOCK_AdTe", number = Inf)
MOCK_AdPre = topTable(fitDupCor, coef="MOCK_AdPre", number = Inf)
MOCK_PreTe = topTable(fitDupCor, coef="MOCK_PreTe", number = Inf)
GLUC_AdTe = topTable(fitDupCor, coef="GLUC_AdTe", number = Inf)
GLUC_AdPre = topTable(fitDupCor, coef="GLUC_AdPre", number = Inf)
GLUC_PreTe = topTable(fitDupCor, coef="GLUC_PreTe", number = Inf)

Adult_Effect = subset(Adult_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #37
Preterm_Effect = subset(Preterm_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #3
Term_Effect = subset(Term_Effect, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #1
MOCK_AdTe = subset(MOCK_AdTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #59
MOCK_AdPre = subset(MOCK_AdPre, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #86
MOCK_PreTe = subset(MOCK_PreTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #0
GLUC_AdTe = subset(GLUC_AdTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #83
GLUC_AdPre= subset(GLUC_AdPre, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #354
GLUC_PreTe = subset(GLUC_PreTe, adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) #0
```

### Retrieval of differentially expressed lncRNAs for respective comparisons of interest

``` r
#BETA adults vs. neonates lncRNA
GLUC_AdTe_vector <- rownames(GLUC_AdTe)
GLUC_AdPre_vector <- rownames(GLUC_AdPre)
GLUC_neo_vector <- c(GLUC_AdTe_vector, GLUC_AdPre_vector) %>% unique()

GLUC_neo_vector <- str_replace(GLUC_neo_vector,
                               pattern = ".[0-9]+$", 
                               replacement = "")

lncs_genes_diff_adult_neonates_GLUC <- rownames(lncs_genes)
lncs_genes_diff_adult_neonates_GLUC <- str_replace(lncs_genes_diff_adult_neonates_GLUC,
                                               pattern = ".[0-9]+$", 
                                               replacement = "")
lncs_genes_diff_adult_neonates_GLUC <- subset(lncs_genes_diff_adult_neonates_GLUC, (lncs_genes_diff_adult_neonates_GLUC %in% GLUC_neo_vector)) %>% as.list()



#MOCK adults vs. neonates lncRNA
MOCK_AdTe_vector <- rownames(MOCK_AdTe)
MOCK_AdPre_vector <- rownames(MOCK_AdPre)

MOCK_neo_vector <- c(MOCK_AdTe_vector, MOCK_AdPre_vector) %>% unique()
MOCK_neo_vector <- str_replace(MOCK_neo_vector,
                              pattern = ".[0-9]+$", 
                              replacement = "")

lncs_genes_diff_adult_neonates_MOCK <- rownames(lncs_genes)
lncs_genes_diff_adult_neonates_MOCK <- str_replace(lncs_genes_diff_adult_neonates_MOCK,
                                                   pattern = ".[0-9]+$", 
                                                   replacement = "")
lncs_genes_diff_adult_neonates_MOCK <- subset(lncs_genes_diff_adult_neonates_MOCK, (lncs_genes_diff_adult_neonates_MOCK %in% MOCK_neo_vector)) %>% as.list()




#Adult treatment diff. exp. lncRNAs
Adult_Effect_vector <- rownames(Adult_Effect)

Adult_Effect_vector <- Adult_Effect_vector %>% unique()
Adult_Effect_vector <- str_replace(Adult_Effect_vector,
                               pattern = ".[0-9]+$", 
                               replacement = "")

lncs_genes_diff_adult_effect <- rownames(lncs_genes)
lncs_genes_diff_adult_effect <- str_replace(lncs_genes_diff_adult_effect,
                                                   pattern = ".[0-9]+$", 
                                                   replacement = "")
lncs_genes_diff_adult_effect <- subset(lncs_genes_diff_adult_effect, (lncs_genes_diff_adult_effect %in% Adult_Effect_vector)) %>% as.list()
```

### Import of lncPath data

``` r
#get lncRNA-mRNA interaction network
NetLncPath <- getNet();
dim(NetLncPath)

#Result_neonates_betaglucan_KEGG <- lncPath(lncs_genes_diff_adult_neonates_GLUC, NetLncPath, Weighted = TRUE, PathwayDataSet = "KEGG", nperm = 1000,
                  #minPathSize = 15, maxPathSize = 500)
#saveRDS(Result_neonates_betaglucan_KEGG, "/home/michael/Trained immunity/RDS/Result_neonates_betaglucan_KEGG.rds")

Result_neonates_betaglucan_KEGG <-rio::import(here::here("data", "Result_neonates_betaglucan_KEGG.rds"))

#Result_neonates_MOCK_KEGG <- lncPath(lncs_genes_diff_adult_neonates_MOCK, NetLncPath, Weighted = TRUE, PathwayDataSet = "KEGG", nperm = 1000,
                                      #minPathSize = 15, maxPathSize = 500)
#saveRDS(Result_neonates_MOCK_KEGG, "C:/Users/Michi/Documents/R-Trained_Immunity/data/RDS_objects/Result_neonates_MOCK_KEGG.rds")

Result_neonates_MOCK_KEGG <-rio::import(here::here("data","Result_neonates_MOCK_KEGG.rds"))

#Result_adults_effect_KEGG <- lncPath(lncs_genes_diff_adult_effect, NetLncPath, Weighted = TRUE, PathwayDataSet = "KEGG", nperm = 1000,
                                      #minPathSize = 15, maxPathSize = 500)
#saveRDS(Result_adults_effect_KEGG, "C:/Users/Michi/Documents/R-Trained_Immunity/data/RDS_objects/Result_adults_effect_KEGG.rds")

Result_adults_effect_KEGG <-rio::import(here::here("data","Result_adults_effect_KEGG.rds"))
```

### Retrieval of pathways

``` r
#get Pathways neonates beta
Neo_beta_KEGG <- lncPath2Table(Result_neonates_betaglucan_KEGG)
Neo_beta_KEGG <- Neo_beta_KEGG %>% dplyr::filter(`False Discovery Rate` <= 0.05)
print(head(Neo_beta_KEGG), row.names = FALSE)
```

    ##                                                  Gene Set Name Gene Set Size
    ##                                 KEGG_OXIDATIVE_PHOSPHORYLATION           114
    ##  KEGG_GLYCOSPHINGOLIPID_BIOSYNTHESIS_LACTO_AND_NEOLACTO_SERIES            23
    ##                                               KEGG_SPLICEOSOME           124
    ##                               KEGG_CHEMOKINE_SIGNALING_PATHWAY           186
    ##                                                KEGG_CELL_CYCLE           122
    ##                      KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION           113
    ##  Enrichment Scores Normalized Enrichment Scores P Value False Discovery Rate
    ##             0.8324             1.39255377167857       0                    0
    ##            0.89414             1.49584097718486       0                    0
    ##            0.76785             1.28456560978303       0                    0
    ##            0.68477             1.14557790272986       0                    0
    ##            0.74068             1.23911187843211       0                    0
    ##            0.74491             1.24618840708924       0                    0

``` r
#get Pathways neonates mock
Neo_mock_KEGG <- lncPath2Table(Result_neonates_MOCK_KEGG)
Neo_mock_KEGG <- Neo_mock_KEGG %>% dplyr::filter(`False Discovery Rate` <= 0.05)
print(head(Neo_mock_KEGG), row.names = FALSE)
```

    ##                     Gene Set Name Gene Set Size Enrichment Scores
    ##    KEGG_OXIDATIVE_PHOSPHORYLATION           114           0.88745
    ##        KEGG_PYRIMIDINE_METABOLISM            94           0.76852
    ##              KEGG_DNA_REPLICATION            36           0.89504
    ##   KEGG_NUCLEOTIDE_EXCISION_REPAIR            44           0.90852
    ##              KEGG_MISMATCH_REPAIR            23           0.89066
    ##  KEGG_CHEMOKINE_SIGNALING_PATHWAY           186           0.66463
    ##  Normalized Enrichment Scores P Value False Discovery Rate
    ##              1.50235460110479       0                    0
    ##              1.30101927775205       0                    0
    ##              1.51520363082183       0                    0
    ##              1.53802377846158       0                    0
    ##               1.5077887757282       0                    0
    ##              1.12514500933267       0                    0

``` r
#get Pathways adults effect
Adults_effect_KEGG <- lncPath2Table(Result_adults_effect_KEGG)
Adults_effect_KEGG <- Adults_effect_KEGG %>% dplyr::filter(`False Discovery Rate` <= 0.05)
print(head(Adults_effect_KEGG), row.names = FALSE)
```

    ##                                     Gene Set Name Gene Set Size
    ##  KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM            43
    ##                               KEGG_FOCAL_ADHESION           198
    ##            KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY           107
    ##            KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY            74
    ##                    KEGG_OXIDATIVE_PHOSPHORYLATION           114
    ##                      KEGG_NOTCH_SIGNALING_PATHWAY            46
    ##  Enrichment Scores Normalized Enrichment Scores P Value False Discovery Rate
    ##             0.8252             1.43683440734685       0                    0
    ##            0.66035             1.14979835299502       0                    0
    ##            0.76814             1.33748180036283       0                    0
    ##             0.7837             1.36457480009419       0                    0
    ##            0.75954             1.32250752030565   0.001               0.0173
    ##            0.81903              1.4260912320035   0.001               0.0173

``` r
#use KEGG as proposed in workflow and depict selected pathways < 0.05 FDR
#top 5 most biologically relevant per comparison

#ß-glucan neonates
colnames(Neo_beta_KEGG) <- c("Geneset", "SetSize", "EnrichmentScore", "nEnrichmentScore", "p.value", "FDR")
Neo_beta_KEGG <- Neo_beta_KEGG %>% dplyr::filter(Geneset == "KEGG_OXIDATIVE_PHOSPHORYLATION" |
                                                 Geneset == "KEGG_GLYCOSPHINGOLIPID_BIOSYNTHESIS_LACTO_AND_NEOLACTO_SERIES" |
                                                 Geneset == "KEGG_SPLICEOSOME" |
                                                 Geneset == "KEGG_CHEMOKINE_SIGNALING_PATHWAY" |
                                                 Geneset == "KEGG_PENTOSE_PHOSPHATE_PATHWAY")
Neo_beta_KEGG <- Neo_beta_KEGG %>% mutate_at(c(2:6), as.numeric)
Neo_beta_KEGG <- Neo_beta_KEGG %>% mutate_at(c(1), as.character)

genes_beta_neo_KEGG_OXIDATIVE_PHOSPHORYLATION<- geneSetDetail(Result_neonates_betaglucan_KEGG, Name = "KEGG_OXIDATIVE_PHOSPHORYLATION")
head(genes_beta_neo_KEGG_OXIDATIVE_PHOSPHORYLATION)
```

    ##   #   GENE LIST LOC     S2N    RES CORE_ENRICHMENT
    ## 1 1 NDUFB1       49  0.0314 0.0717             YES
    ## 2 2 NDUFB3       50  0.0314  0.146             YES
    ## 3 3  ATP5O       51  0.0314   0.22             YES
    ## 4 4 NDUFC1       56  0.0297   0.29             YES
    ## 5 5  UQCRB      212 0.00864  0.302             YES
    ## 6 6  ATP5J      313 0.00634  0.312             YES

``` r
genes_beta_neo_KEGG_GLYCOSPHINGOLIPID_BIOSYNTHESIS_LACTO_AND_NEOLACTO_SERIES <- geneSetDetail(Result_neonates_betaglucan_KEGG, Name = "KEGG_GLYCOSPHINGOLIPID_BIOSYNTHESIS_LACTO_AND_NEOLACTO_SERIES")
head(genes_beta_neo_KEGG_GLYCOSPHINGOLIPID_BIOSYNTHESIS_LACTO_AND_NEOLACTO_SERIES)
```

    ##   #    GENE LIST LOC     S2N   RES CORE_ENRICHMENT
    ## 1 1 ST3GAL4      678 0.00425 0.063             YES
    ## 2 2 B4GALT1      730 0.00405 0.155             YES
    ## 3 3 ST3GAL3      738 0.00402 0.248             YES
    ## 4 4  B3GNT2      747 0.00402 0.342             YES
    ## 5 5 ST3GAL6      765 0.00398 0.434             YES
    ## 6 6 B4GALT4      774 0.00397 0.526             YES

``` r
genes_beta_neo_KEGG_SPLICEOSOME <- geneSetDetail(Result_neonates_betaglucan_KEGG, Name = "KEGG_SPLICEOSOME")
head(genes_beta_neo_KEGG_SPLICEOSOME)
```

    ##   #    GENE LIST LOC     S2N   RES CORE_ENRICHMENT
    ## 1 1 PRPF38A        2  0.0885 0.216             YES
    ## 2 2    PPIE       18  0.0444 0.324             YES
    ## 3 3   U2AF1      139    0.01 0.343             YES
    ## 4 4   CHERP      147 0.00973 0.366             YES
    ## 5 5    LSM2      151 0.00967 0.389             YES
    ## 6 6   SRSF4      158 0.00965 0.413             YES

``` r
genes_beta_neo_KEGG_CHEMOKINE_SIGNALING_PATHWAY <- geneSetDetail(Result_neonates_betaglucan_KEGG, Name = "KEGG_CHEMOKINE_SIGNALING_PATHWAY")
head(genes_beta_neo_KEGG_CHEMOKINE_SIGNALING_PATHWAY)
```

    ##   #   GENE LIST LOC     S2N    RES CORE_ENRICHMENT
    ## 1 1  GNG12       22  0.0426 0.0985             YES
    ## 2 2   GNB5       98  0.0131  0.125             YES
    ## 3 3  CXCR6      114  0.0122  0.153             YES
    ## 4 4 PRKACA      169 0.00956  0.172             YES
    ## 5 5 PIK3R1      174 0.00949  0.194             YES
    ## 6 6   RAF1      199 0.00913  0.214             YES

``` r
genes_beta_neo_KEGG_PENTOSE_PHOSPHATE_PATHWAY <- geneSetDetail(Result_neonates_betaglucan_KEGG, Name = "KEGG_PENTOSE_PHOSPHATE_PATHWAY")
head(genes_beta_neo_KEGG_PENTOSE_PHOSPHATE_PATHWAY)
```

    ##   #    GENE LIST LOC     S2N   RES CORE_ENRICHMENT
    ## 1 1   ALDOB       36  0.0339 0.263             YES
    ## 2 2    PGM1       59  0.0296 0.493             YES
    ## 3 3    RPIA      105  0.0131 0.592             YES
    ## 4 4 PRPS1L1      177 0.00948 0.662             YES
    ## 5 5    PFKL      759   0.004 0.663             YES
    ## 6 6    PFKM      937 0.00358 0.681             YES

``` r
#Mock neonates
Neo_mock_KEGG_2 <- Neo_mock_KEGG

colnames(Neo_mock_KEGG) <- c("Geneset", "SetSize", "EnrichmentScore", "nEnrichmentScore", "p.value", "FDR")
Neo_mock_KEGG <- Neo_mock_KEGG %>% dplyr::filter(Geneset == "KEGG_OXIDATIVE_PHOSPHORYLATION" |
                                          Geneset == "KEGG_PYRIMIDINE_METABOLISM" |
                                          Geneset == "KEGG_SPLICEOSOME" |
                                          Geneset == "KEGG_CHEMOKINE_SIGNALING_PATHWAY" |
                                          Geneset == "KEGG_PEROXISOME")
Neo_mock_KEGG <- Neo_mock_KEGG %>% mutate_at(c(2:6), as.numeric)
Neo_mock_KEGG <- Neo_mock_KEGG %>% mutate_at(c(1), as.factor)

genes_mock_neo_KEGG_OXIDATIVE_PHOSPHORYLATION <- geneSetDetail(Result_neonates_MOCK_KEGG, Name = "KEGG_OXIDATIVE_PHOSPHORYLATION")
head(genes_mock_neo_KEGG_OXIDATIVE_PHOSPHORYLATION)
```

    ##   #    GENE LIST LOC     S2N    RES CORE_ENRICHMENT
    ## 1 1  NDUFB1       23   0.047 0.0937             YES
    ## 2 2  NDUFB3       24   0.047  0.189             YES
    ## 3 3   ATP5O       25   0.047  0.284             YES
    ## 4 4   ATP5J      102 0.00913  0.298             YES
    ## 5 5  UQCR11      271 0.00641  0.302             YES
    ## 6 6 UQCRFS1      312 0.00579  0.311             YES

``` r
genes_mock_neo_KEGG_PYRIMIDINE_METABOLISM <- geneSetDetail(Result_neonates_MOCK_KEGG, Name = "KEGG_PYRIMIDINE_METABOLISM")
head(genes_mock_neo_KEGG_PYRIMIDINE_METABOLISM)
```

    ##   #   GENE LIST LOC     S2N   RES CORE_ENRICHMENT
    ## 1 1   TYMS       30  0.0403 0.173             YES
    ## 2 2   NME1      138 0.00848 0.204             YES
    ## 3 3   NT5E      175 0.00825 0.238             YES
    ## 4 4 POLR2D      238 0.00687 0.264             YES
    ## 5 5 POLR2K      245 0.00665 0.292             YES
    ## 6 6 POLR2F      246 0.00663 0.321             YES

``` r
genes_mock_neo_KEGG_SPLICEOSOME <- geneSetDetail(Result_neonates_MOCK_KEGG, Name = "KEGG_SPLICEOSOME")
head(genes_mock_neo_KEGG_SPLICEOSOME)
```

    ##   #   GENE LIST LOC     S2N   RES CORE_ENRICHMENT
    ## 1 1   PPIE        4  0.0666 0.253             YES
    ## 2 2   XAB2      279  0.0063 0.263             YES
    ## 3 3  U2AF2      375 0.00523 0.278             YES
    ## 4 4 SRSF10      420  0.0049 0.294             YES
    ## 5 5 PRPF19      515  0.0047 0.307             YES
    ## 6 6   ISY1      555 0.00459 0.322             YES

``` r
genes_mock_neo_KEGG_CHEMOKINE_SIGNALING_PATHWAY <- geneSetDetail(Result_neonates_MOCK_KEGG, Name = "KEGG_CHEMOKINE_SIGNALING_PATHWAY")
head(genes_mock_neo_KEGG_CHEMOKINE_SIGNALING_PATHWAY)
```

    ##   #  GENE LIST LOC     S2N    RES CORE_ENRICHMENT
    ## 1 1 GNAI2      192 0.00818 0.0158             YES
    ## 2 2  GNB1      193 0.00818 0.0419             YES
    ## 3 3 GNAI1      195 0.00805 0.0676             YES
    ## 4 4 STAT3      239 0.00679 0.0869             YES
    ## 5 5   LYN      301 0.00599  0.103             YES
    ## 6 6 CXCR4      322 0.00567   0.12             YES

``` r
genes_mock_neo_KEGG_PEROXISOME <- geneSetDetail(Result_neonates_MOCK_KEGG, Name = "KEGG_PEROXISOME")
head(genes_mock_neo_KEGG_PEROXISOME)
```

    ##   #  GENE LIST LOC     S2N   RES CORE_ENRICHMENT
    ## 1 1  PEX6        1  0.0942 0.438             YES
    ## 2 2 HMGCL       43  0.0244  0.55             YES
    ## 3 3 PEX26       57  0.0163 0.625             YES
    ## 4 4  PEX1       58  0.0163 0.701             YES
    ## 5 5   CAT      576 0.00454 0.695             YES
    ## 6 6  SOD1      620 0.00442 0.713             YES

``` r
#ß-Glucan adults
Adults_effect_KEGG_2 <- Adults_effect_KEGG

colnames(Adults_effect_KEGG) <- c("Geneset", "SetSize", "EnrichmentScore", "nEnrichmentScore", "p.value", "FDR")
Adults_effect_KEGG <- Adults_effect_KEGG %>% dplyr::filter(Geneset == "KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM" |
                                                    Geneset == "KEGG_SPLICEOSOME" |
                                                    Geneset == "KEGG_OXIDATIVE_PHOSPHORYLATION" |
                                                    Geneset== "KEGG_CHEMOKINE_SIGNALING_PATHWAY" |
                                                    Geneset == "KEGG_REGULATION_OF_ACTIN_CYTOSKELETON")
Adults_effect_KEGG <- Adults_effect_KEGG %>% mutate_at(c(2:6), as.numeric)
Adults_effect_KEGG <- Adults_effect_KEGG %>% mutate_at(c(1), as.factor)

genes_adults_effect_KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM <- geneSetDetail(Result_adults_effect_KEGG, Name = "KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM")
head(genes_adults_effect_KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM)
```

    ##   #  GENE LIST LOC     S2N   RES CORE_ENRICHMENT
    ## 1 1 CHIT1        6  0.0695 0.463             YES
    ## 2 2  PGM1       23  0.0397 0.727             YES
    ## 3 3   HK2      701 0.00373 0.715             YES
    ## 4 4  UGP2      796 0.00356 0.734             YES
    ## 5 5   GCK      815  0.0035 0.756             YES
    ## 6 6   HK3      819 0.00349 0.779             YES

``` r
genes_adults_effect_KEGG_SPLICEOSOME <- geneSetDetail(Result_adults_effect_KEGG, Name = "KEGG_SPLICEOSOME")
head(genes_adults_effect_KEGG_SPLICEOSOME)
```

    ##   #    GENE LIST LOC     S2N    RES CORE_ENRICHMENT
    ## 1 1    SNW1      109  0.0107 0.0449             YES
    ## 2 2  HSPA1L      112  0.0106 0.0951             YES
    ## 3 3 PRPF40A      316 0.00627  0.114             YES
    ## 4 4    LSM3      417 0.00491  0.132             YES
    ## 5 5   RBM8A      436 0.00472  0.153             YES
    ## 6 6  SNRPB2      447 0.00469  0.175             YES

``` r
genes_adults_effect_KEGG_OXIDATIVE_PHOSPHORYLATION <- geneSetDetail(Result_adults_effect_KEGG, Name = "KEGG_OXIDATIVE_PHOSPHORYLATION")
head(genes_adults_effect_KEGG_OXIDATIVE_PHOSPHORYLATION)
```

    ##   #    GENE LIST LOC     S2N   RES CORE_ENRICHMENT
    ## 1 1  NDUFC1       27  0.0396 0.149             YES
    ## 2 2   UQCRB      107  0.0108 0.185             YES
    ## 3 3  ATP5G2      388 0.00536  0.19             YES
    ## 4 4  NDUFS4      405 0.00511 0.209             YES
    ## 5 5  NDUFA8      570 0.00422 0.216             YES
    ## 6 6 NDUFA11      572 0.00421 0.232             YES

``` r
genes_adults_effect_KEGG_CHEMOKINE_SIGNALING_PATHWAY <- geneSetDetail(Result_adults_effect_KEGG, Name = "KEGG_CHEMOKINE_SIGNALING_PATHWAY")
head(genes_adults_effect_KEGG_CHEMOKINE_SIGNALING_PATHWAY)
```

    ##   #  GENE LIST LOC     S2N     RES CORE_ENRICHMENT
    ## 1 1 MAPK3      284 0.00705 0.00561             YES
    ## 2 2 CDC42      287 0.00679  0.0256             YES
    ## 3 3   LYN      289 0.00661  0.0451             YES
    ## 4 4 CXCL5      306 0.00635   0.063             YES
    ## 5 5 ADCY7      314 0.00628  0.0812             YES
    ## 6 6  GRB2      408 0.00506  0.0912             YES

``` r
genes_adults_effect_KEGG_REGULATION_OF_ACTIN_CYTOSKELETON <- geneSetDetail(Result_adults_effect_KEGG, Name = "KEGG_REGULATION_OF_ACTIN_CYTOSKELETON")
head(genes_adults_effect_KEGG_REGULATION_OF_ACTIN_CYTOSKELETON)
```

    ##   #    GENE LIST LOC     S2N     RES CORE_ENRICHMENT
    ## 1 1      F2      243 0.00749 0.00659             YES
    ## 2 2   ACTN4      254 0.00739  0.0254             YES
    ## 3 3   MAPK3      284 0.00705  0.0424             YES
    ## 4 4   CDC42      287 0.00679  0.0601             YES
    ## 5 5 NCKAP1L      290 0.00655  0.0772             YES
    ## 6 6 ARHGEF6      292 0.00653  0.0943             YES

## Figure 5K: Visualization of pathways

``` r
#ß-glucan neonates vs. adults
lncpath_neo_beta_KEGG <- ggplot(Neo_beta_KEGG, aes(x = nEnrichmentScore, y = fct_reorder(Geneset, nEnrichmentScore))) +
  geom_point(aes(size = SetSize)) +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  ggtitle("KEGG-lnc-neo_beta") + xlim(1,1.6) +
  theme(axis.text.y=element_text(size=7), 
        axis.text.x=element_text(size=7), 
        axis.title = element_text(size = 10, face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 10))

print(lncpath_neo_beta_KEGG)
```

![](Trained_immunity_GIT_files/figure-gfm/setup64-1.png)<!-- -->

``` r
graph2svg(lncpath_neo_beta_KEGG, file = here::here("plots","lncpath_neo_beta_KEGG"), width = 5.85, height = 1.7)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/lncpath_neo_beta_KEGG.svg

``` r
#Mock neonates vs. adults
lncpath_neo_mock_KEGG <- ggplot(Neo_mock_KEGG, aes(x = nEnrichmentScore, y = fct_reorder(Geneset, nEnrichmentScore))) +
  geom_point(aes(size = SetSize)) +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  ggtitle("KEGG-lnc-neo_mock") + xlim(1,1.6) +
  theme(axis.text.y=element_text(size=7), 
        axis.text.x=element_text(size=7), 
        axis.title = element_text(size = 10, face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 10))

print(lncpath_neo_mock_KEGG)
```

![](Trained_immunity_GIT_files/figure-gfm/setup64-2.png)<!-- -->

``` r
graph2svg(lncpath_neo_mock_KEGG, file = here::here("plots","lncpath_neo_beta_KEGG"), width = 4.13, height = 1.7)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/lncpath_neo_beta_KEGG.svg

``` r
#ß-glucan adults vs. mock adults
lncpath_adult_effect_KEGG <- ggplot(Adults_effect_KEGG, aes(x = nEnrichmentScore, y = fct_reorder(Geneset, nEnrichmentScore))) +
  geom_point(aes(size = SetSize)) +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  ggtitle("KEGG-lnc-adult_effect") + xlim(1,1.6) +
  theme(axis.text.y=element_text(size=7), 
        axis.text.x=element_text(size=7), 
        axis.title = element_text(size = 10, face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 10))

print(lncpath_adult_effect_KEGG)
```

![](Trained_immunity_GIT_files/figure-gfm/setup64-3.png)<!-- -->

``` r
graph2svg(lncpath_adult_effect_KEGG, file = here::here("plots","lncpath_neo_beta_KEGG"), width = 5.16, height = 1.7)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/lncpath_neo_beta_KEGG.svg

# WGCNA analysis

### Load necessary packages

``` r
pacman::p_load(WGCNA)
pacman::p_load(gplots)
pacman::p_load(Biobase)
pacman::p_load(mixtools)
pacman::p_load(genefilter)
pacman::p_load(DESeq2)
pacman::p_load(org.Hs.eg.db)
pacman::p_load(clusterProfiler)
pacman::p_load(msigdbr)
pacman::p_load(ReactomePA)
```

## Data import, wrangling and initiation

``` r
#import data
metadata <- readRDS(here::here("data","metadata.rds"))
se <- readRDS(here::here("data","se.rds"))

#short wrangling
metadata <- metadata %>% mutate_at(c(1:5), as.factor) %>% dplyr::select(-Index)

#add metadata
colData(se) <- DataFrame(metadata)
colData(se)
```

    ## DataFrame with 28 rows and 4 columns
    ##     SampleName    Group  Treatment  Index_2
    ##       <factor> <factor>   <factor> <factor>
    ## 1      LW01001    adult mock              1
    ## 2      LW01003    adult mock              2
    ## 3      LW01004    adult betaglucan        2
    ## 4      LW01005    adult mock              3
    ## 5      LW01007    adult mock              4
    ## ...        ...      ...        ...      ...
    ## 24     LW01026  preterm betaglucan        3
    ## 25     LW01027  adult   betaglucan        1
    ## 26     LW01028  adult   betaglucan        3
    ## 27     LW01029  adult   mock              6
    ## 28     LW01030  adult   betaglucan        6

``` r
#rownames explicitly assigned
rownames(colData(se)) <- colData(se)$SampleName

#initiate DESeq-object
dds <- DESeqDataSet(se, design = ~ 1)

#ENSEMBL - get rid of versions
rownames(dds) <- str_sub(rownames(dds), 1, 15)

#extract counts
counts <- assay(dds) %>%
  data.frame

#extract pdata
pdata <- colData(dds) %>%
  data.frame

#create map between ENSEMBLE and ENTREZID
map <- AnnotationDbi::select(org.Hs.eg.db,
                             columns = c("ENTREZID",
                                         "ENSEMBL"),
                             keys = keys(org.Hs.eg.db, keytype = "ENTREZID")) %>%
  drop_na
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
#create annotation for features
fdata <- AnnotationDbi::select(org.Hs.eg.db,
                               columns = c("ENTREZID",
                                           "ENSEMBL",
                                           "SYMBOL",
                                           "GENENAME"),
                               keys = keys(org.Hs.eg.db, keytype = "ENTREZID")) %>%
  group_by(ENTREZID) %>%
  summarize(ENSEMBL = paste(unique(ENSEMBL), collapse = ", "),
            SYMBOL = paste(unique(SYMBOL), collapse = ", "),
            GENENAME = paste(unique(GENENAME), collapse = ", ")) %>%
  mutate(rowname = ENTREZID) %>%
  column_to_rownames
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
#convert the ENSEMBLE based count matrix into ENTREZ based counts 
counts_entrezid <- counts %>%
  rownames_to_column("ENSEMBL") %>%
  inner_join(fdata) %>%
  dplyr::select(ENTREZID, starts_with ("LW")) %>%
  dplyr::group_by(ENTREZID) %>%
  summarize_all(sum) %>%
  column_to_rownames("ENTREZID")
```

    ## Joining with `by = join_by(ENSEMBL)`

``` r
#initiate a new DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_entrezid,
                              colData = pdata[colnames(counts_entrezid),],
                              design = ~ 1)

featureData <- fdata[rownames(dds),]
mcols(dds) <- DataFrame(mcols(dds), featureData)

#filter out lowly expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#apply variance stabilizing transformation
dds_vst <- dds %>%
  vst

#IQR filtering
IQR_vds <- assay(dds_vst)

#calculate IQR across rows
extracted.iqr <- apply(IQR_vds, 1, iqr) %>% as.data.frame()

#arrange and filter for most variable genes
extracted.iqr <- extracted.iqr %>% arrange(desc(.))
extracted.iqr$genes <- rownames(extracted.iqr)
extracted.iqr <- extracted.iqr[1:10000,]

filter_vector <- c(extracted.iqr$genes)

#sanity check
filter_vector <- filter_vector %>% unique()
length(filter_vector)
```

    ## [1] 10000

``` r
#filter
dds_vst_IQR <- dds_vst %>% subset(dds_vst@rowRanges@elementMetadata@listData[["ENTREZID"]] %in% filter_vector)

#create an ExpressionSet for WGCNA
eset <- ExpressionSet(assay(dds_vst_IQR))
pData(eset) <- colData(dds_vst_IQR) %>% data.frame
fData(eset) <- rowData(dds_vst_IQR) %>% data.frame
```

### Additional preparations

``` r
#allow multithreading
enableWGCNAThreads(10)
```

    ## Allowing parallel execution with up to 10 working processes.

``` r
options(stringsAsFactors = FALSE)

#additional filtering-step disabled here - that´s why code seems redundant
selected_samples <- sampleNames(eset)
selected_features <- featureNames(eset)

# filter based on ENTREZID and remove duplicates
eset_wgcna <- eset[selected_features,selected_samples]

# extract the expression matrix
expression <- exprs(eset_wgcna)

networkConcepts <- WGCNA::networkConcepts(expression,
                                          power = 2,
                                          networkType = "signed")

connectivity <- networkConcepts$Connectivity
clustercoefficient <- networkConcepts$ClusterCoef
zK <- (connectivity - mean(connectivity))/sqrt(sd(connectivity))
zC <- (clustercoefficient - mean(clustercoefficient))/sqrt(sd(clustercoefficient))

outliers <- names(zK)[abs(zK) > 2]

expresssion_outlier_removed <- expression[,!(colnames(expression) %in% outliers)]
```

## Parameters for the blockwisemodules function + network creation

``` r
params <- list(networkType = "signed",
               corType = "bicor",
               maxBlockSize = 16000,
               # TOMType = "signed",
               minModuleSize = 80,
               detectCutHeight = 0.988,
               reassignThreshold = 1e-6,
               mergeCutHeight = 0.282,#0.25
               deepSplit = 2,#3
               numericLabels = TRUE,
               pamStage = TRUE,
               pamRespectsDendro = TRUE,
               verbose = 6,
               saveTOMs = TRUE,
               saveTOMFileBase = "trainedTOM",
               datExpr = t(expresssion_outlier_removed))

powers = c(seq(1,10,by=1), seq(12,20, by=2));

sft = pickSoftThreshold(params$datExpr,
                        corFnc = bicor,
                        RsquaredCut = 0.8,
                        powerVector=powers,
                        networkType = params$networkType,
                        verbose = 6)
```

    ## pickSoftThreshold: will use block size 4473.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 4473 of 10000
    ##    ..working on genes 4474 through 8946 of 10000
    ##    ..working on genes 8947 through 10000 of 10000
    ##    Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1  0.01050  7.22          0.985 5020.00   5020.00 5220.0
    ## 2      2  0.00294  1.53          0.966 2690.00   2690.00 2970.0
    ## 3      3  0.01240 -1.43          0.904 1530.00   1520.00 1880.0
    ## 4      4  0.11000 -2.51          0.889  911.00    892.00 1290.0
    ## 5      5  0.29200 -2.87          0.911  565.00    546.00  923.0
    ## 6      6  0.45000 -2.82          0.935  363.00    345.00  686.0
    ## 7      7  0.56500 -2.70          0.947  241.00    224.00  525.0
    ## 8      8  0.64000 -2.58          0.953  164.00    149.00  411.0
    ## 9      9  0.70500 -2.44          0.964  115.00    102.00  328.0
    ## 10    10  0.73900 -2.34          0.965   82.20     70.40  265.0
    ## 11    12  0.80000 -2.16          0.973   44.60     35.60  180.0
    ## 12    14  0.85900 -2.00          0.989   25.90     18.90  127.0
    ## 13    16  0.86000 -2.00          0.978   15.90     10.60   95.1
    ## 14    18  0.87800 -2.01          0.984   10.20      6.12   75.0
    ## 15    20  0.88600 -2.03          0.987    6.77      3.68   60.6

``` r
collectGarbage();

#pick threshold
params$power <- 12

net <- do.call(blockwiseModules, c(params))
```

    ##  Calculating module eigengenes block-wise from all genes
    ##    Flagging genes and samples with too many missing values...
    ##     ..step 1
    ##  ..Working on block 1 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 1 into file trainedTOM-block.1.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##      ..Going through the merge tree
    ##  
    ##      ..Going through detected branches and marking clusters..
    ##      ..Assigning Tree Cut stage labels..
    ##      ..Assigning PAM stage labels..
    ##      ....assigned 2779 objects to existing clusters.
    ##      ..done.
    ##  ....calculating module eigengenes..
    ##      moduleEigengenes : Working on ME for module 1
    ##       ... 1010 genes
    ##      moduleEigengenes : Working on ME for module 2
    ##       ... 544 genes
    ##      moduleEigengenes : Working on ME for module 3
    ##       ... 511 genes
    ##      moduleEigengenes : Working on ME for module 4
    ##       ... 417 genes
    ##      moduleEigengenes : Working on ME for module 5
    ##       ... 400 genes
    ##      moduleEigengenes : Working on ME for module 6
    ##       ... 380 genes
    ##      moduleEigengenes : Working on ME for module 7
    ##       ... 350 genes
    ##      moduleEigengenes : Working on ME for module 8
    ##       ... 342 genes
    ##      moduleEigengenes : Working on ME for module 9
    ##       ... 339 genes
    ##      moduleEigengenes : Working on ME for module 10
    ##       ... 301 genes
    ##      moduleEigengenes : Working on ME for module 11
    ##       ... 297 genes
    ##      moduleEigengenes : Working on ME for module 12
    ##       ... 261 genes
    ##      moduleEigengenes : Working on ME for module 13
    ##       ... 261 genes
    ##      moduleEigengenes : Working on ME for module 14
    ##       ... 229 genes
    ##      moduleEigengenes : Working on ME for module 15
    ##       ... 216 genes
    ##      moduleEigengenes : Working on ME for module 16
    ##       ... 207 genes
    ##      moduleEigengenes : Working on ME for module 17
    ##       ... 198 genes
    ##      moduleEigengenes : Working on ME for module 18
    ##       ... 192 genes
    ##      moduleEigengenes : Working on ME for module 19
    ##       ... 152 genes
    ##      moduleEigengenes : Working on ME for module 20
    ##       ... 139 genes
    ##      moduleEigengenes : Working on ME for module 21
    ##       ... 126 genes
    ##      moduleEigengenes : Working on ME for module 22
    ##       ... 105 genes
    ##      moduleEigengenes : Working on ME for module 23
    ##       ... 102 genes
    ##      moduleEigengenes : Working on ME for module 24
    ##       ... 92 genes
    ##      moduleEigengenes : Working on ME for module 25
    ##       ... 88 genes
    ##  ....checking kME in modules..
    ##      ..removing 28 genes from module 1 because their KME is too low.
    ##      ..removing 28 genes from module 3 because their KME is too low.
    ##      ..removing 1 genes from module 4 because their KME is too low.
    ##      ..removing 20 genes from module 5 because their KME is too low.
    ##      ..removing 9 genes from module 7 because their KME is too low.
    ##      ..removing 1 genes from module 8 because their KME is too low.
    ##      ..removing 5 genes from module 9 because their KME is too low.
    ##      ..removing 11 genes from module 10 because their KME is too low.
    ##      ..removing 2 genes from module 12 because their KME is too low.
    ##      ..removing 1 genes from module 17 because their KME is too low.
    ##   ..reassigning 1 genes from module 1 to modules with higher KME.
    ##   ..reassigning 2 genes from module 2 to modules with higher KME.
    ##   ..reassigning 1 genes from module 12 to modules with higher KME.
    ##   ..reassigning 4 genes from module 13 to modules with higher KME.
    ##  ..merging modules that are too close..
    ##      mergeCloseModules: Merging modules whose distance is less than 0.282
    ##        .. will look for grey label ME0
    ##        multiSetMEs: Calculating module MEs.
    ##          Working on set 1 ...
    ##          moduleEigengenes : Working on ME for module 1
    ##          moduleEigengenes : Working on ME for module 2
    ##          moduleEigengenes : Working on ME for module 3
    ##          moduleEigengenes : Working on ME for module 4
    ##          moduleEigengenes : Working on ME for module 5
    ##          moduleEigengenes : Working on ME for module 6
    ##          moduleEigengenes : Working on ME for module 7
    ##          moduleEigengenes : Working on ME for module 8
    ##          moduleEigengenes : Working on ME for module 9
    ##          moduleEigengenes : Working on ME for module 10
    ##          moduleEigengenes : Working on ME for module 11
    ##          moduleEigengenes : Working on ME for module 12
    ##          moduleEigengenes : Working on ME for module 13
    ##          moduleEigengenes : Working on ME for module 14
    ##          moduleEigengenes : Working on ME for module 15
    ##          moduleEigengenes : Working on ME for module 16
    ##          moduleEigengenes : Working on ME for module 17
    ##          moduleEigengenes : Working on ME for module 18
    ##          moduleEigengenes : Working on ME for module 19
    ##          moduleEigengenes : Working on ME for module 20
    ##          moduleEigengenes : Working on ME for module 21
    ##          moduleEigengenes : Working on ME for module 22
    ##          moduleEigengenes : Working on ME for module 23
    ##          moduleEigengenes : Working on ME for module 24
    ##          moduleEigengenes : Working on ME for module 25
    ##        Merging original colors 3, 6
    ##        Merging original colors 2, 9
    ##        Merging original colors 8, 17
    ##        Merging original colors 4, 11
    ##        Merging original colors 1, 5
    ##        Merging original colors 18, 19
    ##        multiSetMEs: Calculating module MEs.
    ##          Working on set 1 ...
    ##          moduleEigengenes : Working on ME for module 1
    ##          moduleEigengenes : Working on ME for module 2
    ##          moduleEigengenes : Working on ME for module 3
    ##          moduleEigengenes : Working on ME for module 4
    ##          moduleEigengenes : Working on ME for module 7
    ##          moduleEigengenes : Working on ME for module 8
    ##          moduleEigengenes : Working on ME for module 10
    ##          moduleEigengenes : Working on ME for module 12
    ##          moduleEigengenes : Working on ME for module 13
    ##          moduleEigengenes : Working on ME for module 14
    ##          moduleEigengenes : Working on ME for module 15
    ##          moduleEigengenes : Working on ME for module 16
    ##          moduleEigengenes : Working on ME for module 18
    ##          moduleEigengenes : Working on ME for module 20
    ##          moduleEigengenes : Working on ME for module 21
    ##          moduleEigengenes : Working on ME for module 22
    ##          moduleEigengenes : Working on ME for module 23
    ##          moduleEigengenes : Working on ME for module 24
    ##          moduleEigengenes : Working on ME for module 25
    ##         Changing original colors:
    ##             1 to  1
    ##             2 to  2
    ##             3 to  3
    ##             4 to  4
    ##             8 to  5
    ##             18 to  6
    ##             7 to  7
    ##             10 to  8
    ##             12 to  9
    ##             13 to  10
    ##             14 to  11
    ##             15 to  12
    ##             16 to  13
    ##             20 to  14
    ##             21 to  15
    ##             22 to  16
    ##             23 to  17
    ##             24 to  18
    ##             25 to  19
    ##        Calculating new MEs...
    ##        multiSetMEs: Calculating module MEs.
    ##          Working on set 1 ...
    ##          moduleEigengenes : Working on ME for module 0
    ##          moduleEigengenes : Working on ME for module 1
    ##          moduleEigengenes : Working on ME for module 2
    ##          moduleEigengenes : Working on ME for module 3
    ##          moduleEigengenes : Working on ME for module 4
    ##          moduleEigengenes : Working on ME for module 5
    ##          moduleEigengenes : Working on ME for module 6
    ##          moduleEigengenes : Working on ME for module 7
    ##          moduleEigengenes : Working on ME for module 8
    ##          moduleEigengenes : Working on ME for module 9
    ##          moduleEigengenes : Working on ME for module 10
    ##          moduleEigengenes : Working on ME for module 11
    ##          moduleEigengenes : Working on ME for module 12
    ##          moduleEigengenes : Working on ME for module 13
    ##          moduleEigengenes : Working on ME for module 14
    ##          moduleEigengenes : Working on ME for module 15
    ##          moduleEigengenes : Working on ME for module 16
    ##          moduleEigengenes : Working on ME for module 17
    ##          moduleEigengenes : Working on ME for module 18
    ##          moduleEigengenes : Working on ME for module 19

``` r
# attach underlying parameters and data
net$params <- as.list(args(blockwiseModules))
net$params[names(params)] <- params
net$eset <- eset_wgcna

#save network
saveRDS(net, here::here("export_WGCNA", "WGCNAnet_trained_immunity_110123.rds"))
```

## Visualisation of WGCNA-results

### Dendrogram structure and modules

``` r
plotDendroAndColors(dendro = net$dendrograms[[1]], 
                    colors = cbind(net$unmergedColors, net$colors),
                    groupLabels = c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05, 
                    main = "Consensus gene dendrogram and module colors")
```

![](Trained_immunity_GIT_files/figure-gfm/setup69-1.png)<!-- -->

## Figure 6A: TOMplot - not done here, computationally very heavy

``` r
#import network
net <- rio::import(here::here("export_WGCNA", "WGCNAnet_trained_immunity_110123.rds"))

#TOMplot
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

#calculate topological overlap anew
nGenes = ncol(t(expresssion_outlier_removed))
nSamples = nrow(t(expresssion_outlier_removed))

#all genes
nSelect = 10000

#select + calculate dissTOM
select = sample(nGenes, size = nSelect)
dissTOM = 1-TOMsimilarityFromExpr(t(expresssion_outlier_removed), networkType = "signed", corType = "bicor", power = 12, nThreads = 12)
```

    ## TOM calculation: adjacency..
    ## ..will not use multithreading.
    ##  Fraction of slow calculations: 0.000000
    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.

``` r
#transform dissTOM to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7

#set diagonal to NA for a nicer plot
diag(plotTOM) = NA

#define parameters
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
geneTree = net$dendrograms[[1]];

#call the plot function
sizeGrWindow(9,9)

#adapt coloring
myheatcol = colorpanel(250,'red',"orange","lemonchiffon")

#perform vizualisation (CAVE: takes a lot of time)
#TOMplot_self <- TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes", col=myheatcol)

#print(TOMplot_self)

#save
#graph2svg(TOMplot_self, file = "TOMplot_self", width = 15, height = 15)
```

### Module Trait association preparation

``` r
#wrangling/preparation
pdata <- pData(net$eset) %>%
  mutate(sample_type = paste0(Group, Treatment, SEP = ""))
pdata$sample_type <- as.factor(pdata$sample_type)
traits <- data.frame(pdata) %>% 
  dplyr::select(sample_type)

#binarize traits
allTraits <- WGCNA::binarizeCategoricalColumns(data = traits$sample_type,
                                               includePairwise = FALSE,
                                               includeLevelVsAll = TRUE,
                                               dropFirstLevelVsAll = FALSE,
                                               minCount = 1)
#wrangling/preparation
moduleColors = paste("M", net$colors, sep = "")
nGenes = ncol(net[["params"]][["datExpr"]])
nSamples <- nrow(net$params$datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(net$params$datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
MEs_plot <- MEs %>% as.data.frame()

#Trait Correlations
moduleTraitCor = cor(MEs, allTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)


par(mar= c(6, 8.5, 3, 3))
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(allTraits),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.35,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))
```

![](Trained_immunity_GIT_files/figure-gfm/setup71-1.png)<!-- -->

``` r
# -> a little too many comparisons to adequately interpret. Merge bigger groups and do it again
```

## Figure 6B: Module Trait associations

``` r
#wrangling + preparation
traits <- data.frame(pdata) %>% dplyr::select(Group, Treatment)

traits_group <- WGCNA::binarizeCategoricalColumns(data = traits$Group,
                                               includePairwise = FALSE,
                                               includeLevelVsAll = TRUE,
                                               dropFirstLevelVsAll = FALSE,
                                               minCount = 1)
traits_Treatment <- WGCNA::binarizeCategoricalColumns(data = traits$Treatment,
                                                  includePairwise = FALSE,
                                                  includeLevelVsAll = TRUE,
                                                  dropFirstLevelVsAll = FALSE,
                                                  minCount = 1)
allTraits <-  cbind(traits_Treatment, traits_group)

nGenes = ncol(net[["params"]][["datExpr"]])
nSamples <- nrow(net$params$datExpr)

MEs0 = moduleEigengenes(net$params$datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, allTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(6, 8.5, 3, 3))

my_palette2 <-  colorRampPalette(c("#0765A0","#FFF1E4","#C61700"))(100)

module_traits_overall <- labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(allTraits),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= my_palette2,
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.5,
               cex.lab.x = 0.8,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))

labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(allTraits),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= my_palette2,
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.5,
               cex.lab.x = 0.8,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))
```

![](Trained_immunity_GIT_files/figure-gfm/setup72-1.png)<!-- -->

### Vizualisation of all modules using violin-plots

``` r
#define the groups
pdata <- pData(net$eset) %>%
  mutate(sample_type = paste0(Group, Treatment, SEP = ""))

pdata$sample_type <- as.factor(pdata$sample_type)

groups <- pdata[rownames(net$params$datExpr), "sample_type"]
#define comparisons
comparisons <- list(c("betaglucan", "mock"))

#extract module colors
moduleColors = paste("M", net$colors, sep = "")

#calculate MEs with color labels
data <- WGCNA::moduleEigengenes(net$params$datExpr, moduleColors)$eigengenes %>%
  
  #remove grey
  dplyr::select(-MEM0) %>%
  
  #rename the color names from ME <color> to <color>
  rename_all(~ str_replace(.,"ME", "")) %>%
  
  #add sampleID
  rownames_to_column(var = "SAMPLE") %>%
  
  #add groups
  mutate(GROUP = groups) %>%
  
  #pivot longer
  pivot_longer(!(matches("SAMPLE") | matches("GROUP")),
               names_to = "MODULE",
               values_to = "VALUE") %>%
  
  #add color
  mutate(COLOR = WGCNA::labels2colors(as.numeric(str_replace(MODULE, "M", ""))))

eigengene.plots.violin <- data %>%
  tidyr::nest(GROUP = -"MODULE") %>%
  deframe() %>%
  purrr::map2(names(.),
              function(x,y){
                ggplot(x, aes(x = GROUP,
                              y = VALUE))+ 
                  geom_violin(stat="ydensity",
                              fill = unique(x$COLOR),
                              size = 0.25) +
                  geom_boxplot(width=0.1,
                               size = 0.25,
                               outlier.size=0.1) +
                  ylim(-1, 1) +
                  annotate(geom = "text",
                           x = 1,
                           y = 1,
                           label = paste("Module:", y),
                           size = 2,
                           hjust = 0) +
                  geom_hline(yintercept = 0,
                             linewidth = 0.25) +
                  labs(x = "",
                       y = "eigengene expression") +
                theme_bw(base_size = 10) +
                  theme(axis.text.x = element_text(angle=90))}
)

print(eigengene.plots.violin)
```

    ## $M1

![](Trained_immunity_GIT_files/figure-gfm/setup73-1.png)<!-- -->

    ## 
    ## $M10

![](Trained_immunity_GIT_files/figure-gfm/setup73-2.png)<!-- -->

    ## 
    ## $M11

![](Trained_immunity_GIT_files/figure-gfm/setup73-3.png)<!-- -->

    ## 
    ## $M12

![](Trained_immunity_GIT_files/figure-gfm/setup73-4.png)<!-- -->

    ## 
    ## $M13

![](Trained_immunity_GIT_files/figure-gfm/setup73-5.png)<!-- -->

    ## 
    ## $M14

![](Trained_immunity_GIT_files/figure-gfm/setup73-6.png)<!-- -->

    ## 
    ## $M15

![](Trained_immunity_GIT_files/figure-gfm/setup73-7.png)<!-- -->

    ## 
    ## $M16

![](Trained_immunity_GIT_files/figure-gfm/setup73-8.png)<!-- -->

    ## 
    ## $M17

![](Trained_immunity_GIT_files/figure-gfm/setup73-9.png)<!-- -->

    ## 
    ## $M18

![](Trained_immunity_GIT_files/figure-gfm/setup73-10.png)<!-- -->

    ## 
    ## $M19

![](Trained_immunity_GIT_files/figure-gfm/setup73-11.png)<!-- -->

    ## 
    ## $M2

![](Trained_immunity_GIT_files/figure-gfm/setup73-12.png)<!-- -->

    ## 
    ## $M3

![](Trained_immunity_GIT_files/figure-gfm/setup73-13.png)<!-- -->

    ## 
    ## $M4

![](Trained_immunity_GIT_files/figure-gfm/setup73-14.png)<!-- -->

    ## 
    ## $M5

![](Trained_immunity_GIT_files/figure-gfm/setup73-15.png)<!-- -->

    ## 
    ## $M6

![](Trained_immunity_GIT_files/figure-gfm/setup73-16.png)<!-- -->

    ## 
    ## $M7

![](Trained_immunity_GIT_files/figure-gfm/setup73-17.png)<!-- -->

    ## 
    ## $M8

![](Trained_immunity_GIT_files/figure-gfm/setup73-18.png)<!-- -->

    ## 
    ## $M9

![](Trained_immunity_GIT_files/figure-gfm/setup73-19.png)<!-- -->

## Figure 6C: Analysis of trait assocations over time (only those that also give biol. meaningful results in ORA)

``` r
#prepare data
MEs_plot$SampleName <- rownames(MEs_plot)
ME_plot <- MEs_plot %>% left_join(metadata, by = "SampleName")

ME_plot <- ME_plot %>% mutate(time = case_when(Group == "adult" ~ 6, 
                                               Group == "term" ~ 2.5, 
                                               Group == "preterm" ~ 1, 
                                               TRUE ~ NA))

ME_plot <- ME_plot %>% mutate_at(c(21:24), as.factor)
ME_plot <- ME_plot %>% mutate_at(c(1:20,25), as.numeric)

ME_plot <- ME_plot %>% pivot_longer(cols = 1:20, values_to = "eigengene_value", names_to = "module")
ME_plot <- ME_plot %>% group_by(Group, Treatment, time, module) %>% mutate(mean = mean(eigengene_value), sd = sd(eigengene_value)) %>% ungroup()

ME_plot_mock <- ME_plot %>% dplyr::filter(Treatment == "mock")
ME_plot_beta <- ME_plot %>% dplyr::filter(Treatment == "betaglucan")

#MOCK "over-time"
module_eigengenes_mock <- ggplot(ME_plot_mock, aes(x=time, y=eigengene_value)) + 
  geom_point(aes(color=Group)) + 
  geom_point(aes(x= time, y = mean)) + 
  geom_line(aes(x= time, y = mean)) + 
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), fill = "grey", alpha = 0.35) +
  facet_wrap(~module) +
  theme_bw() + ylim(-0.55,0.55) +  theme(axis.text=element_text(size=8), 
                                       axis.title = element_text(size = 10, face = "bold"), 
                                       legend.title = element_blank(),
                                       legend.text = element_text(size = 10))

#Beta-Glucan "over-time"
module_eigengenes_beta<-  ggplot(ME_plot_beta, aes(x=time, y=eigengene_value)) + 
  geom_point(aes(color=Group)) + 
  geom_point(aes(x= time, y = mean)) + 
  geom_line(aes(x= time, y = mean)) + 
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), fill = "grey", alpha = 0.35) +
  facet_wrap(~module) +
  theme_bw() + ylim(-0.55,0.55) +  theme(axis.text=element_text(size=8), 
                                       axis.title = element_text(size = 10, face = "bold"), 
                                       legend.title = element_blank(),
                                       legend.text = element_text(size = 10))

#save MOCK
graph2svg(module_eigengenes_mock, file = here::here("plots", "module_eigengenes_mock"), width = 6, height = 4.7)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/module_eigengenes_mock.svg

``` r
#save Beta
graph2svg(module_eigengenes_beta, file = here::here("plots", "module_eigengenes_beta"), width = 6, height = 4.7)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/module_eigengenes_beta.svg

## Figure 6D: Overrepresentation Analysis of WGCNA-modules

``` r
#define genes + background
gene_list <- data.frame(gene = names(net$colors),
                        module = net$colors) %>%
  dplyr::mutate(module = paste0("M", module)) %>%
  tidyr::nest(gg = -"module") %>%
  deframe %>%
  purrr::map(deframe)

universe <- names(net$colors)

#loop GO:BP analysis
Modules_WGCNA_GOBP<- list()
for (i in 1:20)  {
  a <- enrichGO(gene = gene_list[[i]],
                     universe = universe,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.1,
                     qvalueCutoff = 0.1)
  Modules_WGCNA_GOBP[[i]] <- a
}

#redefine stuff (unclearly assigned by loop, could ofc be computed easier - but already done manually at the time)
M3 <- Modules_WGCNA_GOBP[[2]]
M12 <- Modules_WGCNA_GOBP[[3]]
M9 <- Modules_WGCNA_GOBP[[4]]
M7 <- Modules_WGCNA_GOBP[[5]]
M1 <- Modules_WGCNA_GOBP[[6]]
M6 <- Modules_WGCNA_GOBP[[7]]
M17 <- Modules_WGCNA_GOBP[[8]]
M2 <- Modules_WGCNA_GOBP[[9]]
M18 <- Modules_WGCNA_GOBP[[10]]
M5 <- Modules_WGCNA_GOBP[[11]]
M4 <- Modules_WGCNA_GOBP[[12]]
M16 <- Modules_WGCNA_GOBP[[13]]
M8 <- Modules_WGCNA_GOBP[[14]]
M14 <- Modules_WGCNA_GOBP[[15]]
M10 <- Modules_WGCNA_GOBP[[16]]
M19 <- Modules_WGCNA_GOBP[[17]]
M15 <- Modules_WGCNA_GOBP[[18]]
M13 <- Modules_WGCNA_GOBP[[19]]
M11 <- Modules_WGCNA_GOBP[[20]]

WGCNA_res_GO_BP <- list("M2" = M2, "M3" = M3, "M4" = M4, "M5"= M5,
                                                   "M6" = M6, "M7" = M7, "M8" = M8, "M9" = M9, 
                                                   "M10" = M10,"M11" = M11, "M12" = M12, "M13" = M13, 
                                                   "M14" = M14, "M15" = M15, "M16" = M16, "M17" = M17,
                                                   "M18" = M18,"M19" = M19)

#make dotplots
BP_M2 <- enrichplot::dotplot(WGCNA_res_GO_BP$M2, showCategory=5) + 
  xlim(0,0.3) + theme(axis.text.y = element_text(size = 6), axis.text.x= element_text(size = 6))
print(BP_M2)
```

![](Trained_immunity_GIT_files/figure-gfm/setup75-1.png)<!-- -->

``` r
graph2svg(BP_M2, 
          file = here::here("plots", "BP_M2"), 
          height = 1.4, width = 3.75)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/BP_M2.svg

``` r
BP_M8 <- enrichplot::dotplot(WGCNA_res_GO_BP$M8, showCategory=9) + 
  xlim(0,0.3) + theme(axis.text.y = element_text(size = 6), axis.text.x= element_text(size = 6))
print(BP_M8)
```

![](Trained_immunity_GIT_files/figure-gfm/setup75-2.png)<!-- -->

``` r
#not all make biol. sense here -> thats why cat. = 9 and those that make sense were considered for the figure
graph2svg(BP_M8, 
          file = here::here("plots", "BP_M8"), 
          height = 1.9, width = 3.24)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/BP_M8.svg

``` r
BP_M10 <- enrichplot::dotplot(WGCNA_res_GO_BP$M10, showCategory=5) + 
  xlim(0,0.3) + theme(axis.text.y = element_text(size = 6), axis.text.x= element_text(size = 6))
print(BP_M10)
```

![](Trained_immunity_GIT_files/figure-gfm/setup75-3.png)<!-- -->

``` r
graph2svg(BP_M10, 
          file = here::here("plots", "BP_M10"), 
          height = 1.6, width = 3.05)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/BP_M10.svg

``` r
BP_M12 <- enrichplot::dotplot(WGCNA_res_GO_BP$M12, showCategory=5) + 
  xlim(0,0.3) + theme(axis.text.y = element_text(size = 6), axis.text.x= element_text(size = 30))
print(BP_M12)
```

![](Trained_immunity_GIT_files/figure-gfm/setup75-4.png)<!-- -->

``` r
graph2svg(BP_M12, 
          file = here::here("plots", "BP_M12"), 
          height = 1.7, width =3.56)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/BP_M12.svg

### Identification of hub-genes from biologically relevant modules

``` r
#calculate connectivity, identify module hub-genes
moduleColors <- net$colors %>%
  paste0("M", .)
adjacency <- WGCNA::adjacency(net$params$datExpr,
                              power = net$params$power,
                              type = net$params$networkType,
                              corFnc = net$params$corType)
konnectivity <- WGCNA::intramodularConnectivity(adjacency, moduleColors) %>%
  mutate(module = moduleColors) %>%
  rownames_to_column("ENTREZID")
konnectivity_norm <- konnectivity %>%
  tidyr::nest(gg = -"module") %>%
  deframe %>%
  purrr::map(function(df){
    df %>%
      mutate(kWithin.norm = kWithin/max(kWithin))
  })

#gene annotations
map <- AnnotationDbi::select(org.Hs.eg.db,
                             columns = c("ENTREZID",
                                         "ENSEMBL"),
                             keys = keys(org.Hs.eg.db, keytype = "ENTREZID")) %>%
  drop_na
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
#create  annotation for features
fdata <- AnnotationDbi::select(org.Hs.eg.db,
                               columns = c("ENTREZID",
                                           "ENSEMBL",
                                           "SYMBOL",
                                           "GENENAME"),
                               keys = keys(org.Hs.eg.db, keytype = "ENTREZID")) %>%
  group_by(ENTREZID) %>%
  summarize(ENSEMBL = paste(unique(ENSEMBL), collapse = ", "),
            SYMBOL = paste(unique(SYMBOL), collapse = ", "),
            GENENAME = paste(unique(GENENAME), collapse = ", ")) %>%
  mutate(rowname = ENTREZID) %>%
  column_to_rownames
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
#signed KME
eigengenes <- WGCNA::moduleEigengenes(net$params$datExpr, moduleColors)$eigengenes
sKME <- WGCNA::signedKME(net$params$datExpr, eigengenes)
sKME$ENTREZID <- rownames(sKME)

#module membership

#M2
M2 <- konnectivity_norm$M2
M2 <- M2 %>% left_join(fdata, by = "ENTREZID")
sKME_2 <- sKME %>% dplyr::select(kMEM2, ENTREZID)
M2 <- M2 %>% left_join(sKME_2, by = "ENTREZID")
M2 <- M2 %>% arrange(desc(kMEM2))
#M2 <- M2 %>% arrange(desc(kWithin.norm))
print(c(M2$ENTREZID))
```

    ##   [1] "51514"     "7083"      "55388"     "6241"      "22974"     "3161"     
    ##   [7] "9768"      "157570"    "991"       "983"       "890"       "8318"     
    ##  [13] "54443"     "332"       "374393"    "9156"      "7153"      "64151"    
    ##  [19] "11004"     "195828"    "5888"      "990"       "4605"      "6790"     
    ##  [25] "55165"     "9133"      "24137"     "51659"     "9787"      "55635"    
    ##  [31] "891"       "221150"    "29128"     "2305"      "81624"     "79019"    
    ##  [37] "5347"      "10403"     "55355"     "79733"     "9833"      "81620"    
    ##  [43] "10112"     "91687"     "4085"      "79022"     "54821"     "3832"     
    ##  [49] "7298"      "993"       "55789"     "4751"      "1719"      "54478"    
    ##  [55] "1666"      "11065"     "150468"    "79801"     "5832"      "7465"     
    ##  [61] "9212"      "51512"     "113130"    "29028"     "23594"     "1869"     
    ##  [67] "10615"     "7390"      "995"       "2597"      "79843"     "4171"     
    ##  [73] "83879"     "3070"      "1164"      "83461"     "11130"     "2491"     
    ##  [79] "9928"      "9319"      "9088"      "2791"      "4173"      "2237"     
    ##  [85] "79968"     "10733"     "6240"      "147841"    "55005"     "6491"     
    ##  [91] "63967"     "4175"      "699"       "90488"     "151648"    "55143"    
    ##  [97] "79172"     "5427"      "157313"    "83990"     "64105"     "81610"    
    ## [103] "55872"     "8208"      "57650"     "9582"      "7272"      "29899"    
    ## [109] "91057"     "580"       "9232"      "4288"      "259266"    "79980"    
    ## [115] "10721"     "701"       "84791"     "10024"     "9134"      "5352"     
    ## [121] "5983"      "2882"      "162681"    "1031"      "57082"     "90381"    
    ## [127] "10036"     "83540"     "898"       "2634"      "4176"      "94081"    
    ## [133] "29127"     "150771"    "146909"    "55010"     "8438"      "8309"     
    ## [139] "728833"    "101928132" "63875"     "51611"     "63979"     "4998"     
    ## [145] "81930"     "554282"    "51661"     "7516"      "729533"    "81563"    
    ## [151] "1063"      "1033"      "26147"     "79682"     "11113"     "90843"    
    ## [157] "1861"      "27306"     "116159"    "159371"    "57545"     "79000"    
    ## [163] "8801"      "7169"      "55055"     "1376"      "129642"    "64757"    
    ## [169] "55247"     "378708"    "4174"      "115106"    "56947"     "151827"   
    ## [175] "8099"      "653820"    "29940"     "339929"    "1058"      "51056"    
    ## [181] "5985"      "55277"     "3014"      "23743"     "152137"    "84296"    
    ## [187] "81853"     "401397"    "672"       "5901"      "200035"    "64756"    
    ## [193] "646324"    "51371"     "10026"     "8323"      "5984"      "51071"    
    ## [199] "11034"     "55013"     "5557"      "29789"     "112399"    "79723"    
    ## [205] "5608"      "79075"     "11274"     "105370333" "8914"      "51390"    
    ## [211] "11077"     "137392"    "5358"      "29980"     "81691"     "3094"     
    ## [217] "57124"     "105379362" "117854"    "60526"     "974"       "202052"   
    ## [223] "51699"     "9837"      "1350"      "119710"    "54841"     "54892"    
    ## [229] "963"       "5723"      "89858"     "283487"    "101926907" "91368"    
    ## [235] "84817"     "5230"      "65008"     "23658"     "5558"      "64080"    
    ## [241] "25902"     "104310351" "114883"    "219865"    "92667"     "185"      
    ## [247] "727"       "6119"      "135154"    "28985"     "55215"     "5255"     
    ## [253] "541471"    "5111"      "81611"     "4258"      "51248"     "5033"     
    ## [259] "497258"    "56952"     "2335"      "10982"     "130589"    "53918"    
    ## [265] "92014"     "55732"     "144455"    "7867"      "5756"      "134147"   
    ## [271] "102724661" "5148"      "8794"      "153642"    "92092"     "4172"     
    ## [277] "51768"     "401258"    "284403"    "80777"     "10478"     "100506881"
    ## [283] "7388"      "26150"     "5464"      "5019"      "6746"      "101928687"
    ## [289] "105378305" "1111"      "55163"     "27067"     "169834"    "9532"     
    ## [295] "80071"     "340481"    "101927978" "26586"     "151056"    "6126"     
    ## [301] "84259"     "58527"     "26499"     "641"       "5950"      "10739"    
    ## [307] "27075"     "159091"    "5743"      "122769"    "131076"    "79989"    
    ## [313] "2177"      "23741"     "79807"     "55825"     "51053"     "3015"     
    ## [319] "7512"      "55775"     "93627"     "51115"     "28976"     "132320"   
    ## [325] "5889"      "6772"      "285343"    "643155"    "54908"     "56704"    
    ## [331] "106479038" "4430"      "1163"      "60"        "390940"    "55706"    
    ## [337] "79600"     "10276"     "9700"      "5158"      "1368"      "113115"   
    ## [343] "116832"    "22824"     "1910"      "84650"     "57697"     "389119"   
    ## [349] "471"       "9055"      "55081"     "10484"     "11142"     "2673"     
    ## [355] "51527"     "10247"     "8872"      "27284"     "11147"     "3835"     
    ## [361] "4130"      "87769"     "4038"      "664"       "11257"     "501"      
    ## [367] "1776"      "126820"    "51075"     "3071"      "2958"      "284252"   
    ## [373] "41"        "101928136" "9553"      "58505"     "10381"     "149837"   
    ## [379] "348235"    "397"       "79858"     "57447"     "50804"     "144608"   
    ## [385] "8228"      "92340"     "56954"     "80312"     "51255"     "64078"    
    ## [391] "28977"     "5686"      "5163"      "140688"    "54414"     "25886"    
    ## [397] "10203"     "148741"    "3329"      "10566"     "8787"      "23516"    
    ## [403] "9455"      "11164"     "7381"      "647946"    "22917"     "56253"    
    ## [409] "11107"     "83463"     "9459"      "2944"      "4833"      "26271"    
    ## [415] "574016"    "116028"    "254948"    "192683"    "644100"    "9601"     
    ## [421] "132299"    "60592"     "100188953" "105376159" "64801"     "4001"     
    ## [427] "346389"    "892"       "51338"     "101929541" "64374"     "6354"     
    ## [433] "128854"    "2001"      "5887"      "79077"     "79814"     "2670"     
    ## [439] "9401"      "83937"     "125144"    "1854"      "79823"     "7167"     
    ## [445] "85437"     "81618"     "1603"      "1282"      "65062"     "105369391"
    ## [451] "57181"     "51657"     "57414"     "2946"      "648791"    "83856"    
    ## [457] "4925"      "2996"      "154664"    "8611"      "8395"      "8029"     
    ## [463] "91614"     "51019"     "644656"    "115827"    "205327"    "319138"   
    ## [469] "4353"      "54969"     "55711"     "4600"      "2187"      "729178"   
    ## [475] "23557"     "6450"      "6342"      "122622"    "102724163" "24145"    
    ## [481] "2203"      "58511"     "10592"     "55144"     "58499"     "100130967"
    ## [487] "3006"      "26521"     "128346"    "3081"      "3606"      "4717"     
    ## [493] "79582"     "100506119" "50808"     "286151"    "55615"     "93550"    
    ## [499] "10570"     "89958"     "23530"     "23464"     "85824"     "2706"     
    ## [505] "30001"     "55766"     "666"       "220042"    "25874"     "51315"    
    ## [511] "3251"      "80221"     "79736"     "105378425" "64839"     "81839"    
    ## [517] "6136"      "550643"    "10745"     "522"       "7837"      "5939"     
    ## [523] "140564"    "84693"     "10953"     "79140"     "219285"    "203328"   
    ## [529] "4437"      "54492"     "26228"     "160365"    "100506098" "80008"    
    ## [535] "388389"    "55786"     "120224"    "79414"     "643617"    "100131211"
    ## [541] "139324"    "5127"      "92960"     "634"       "140462"    "23643"    
    ## [547] "9585"      "6915"      "340120"    "10669"     "348751"    "55166"    
    ## [553] "8364"      "22883"     "8819"      "54107"     "25788"     "5213"     
    ## [559] "5368"      "3796"      "79031"     "121355"    "84274"     "84833"    
    ## [565] "149628"    "94031"     "79053"     "4608"      "101927111" "10923"    
    ## [571] "51501"     "9562"      "7108"      "84937"     "4724"      "5627"     
    ## [577] "167555"    "84515"     "38"        "57573"     "50636"     "1846"     
    ## [583] "26548"     "23600"     "101928101" "1996"      "8975"      "55843"    
    ## [589] "3217"      "3667"      "9141"      "83695"     "2842"      "115123"   
    ## [595] "90557"     "401152"    "101928227" "101927888" "4072"      "55771"    
    ## [601] "5366"      "441549"    "6285"      "401494"    "143471"    "2272"     
    ## [607] "23708"     "132321"    "266655"    "57576"     "105378049" "5479"     
    ## [613] "9918"      "1175"      "79816"     "909"       "5372"      "10440"    
    ## [619] "154467"    "7643"      "83903"     "83594"     "79865"     "219736"   
    ## [625] "1277"      "55052"     "338099"    "10412"     "538"       "92521"    
    ## [631] "5357"      "160897"    "54976"     "89894"     "101954278" "340351"   
    ## [637] "102723663" "441051"    "3396"      "55520"     "55300"     "1291"     
    ## [643] "26472"     "136332"    "80023"     "140465"    "441024"    "105371814"
    ## [649] "2995"      "1843"      "3601"      "1902"      "84217"     "6839"     
    ## [655] "83734"     "2004"      "101928104" "9702"      "252839"    "27068"    
    ## [661] "25804"     "100507392" "230"       "758"       "29058"     "1044"     
    ## [667] "145957"    "160428"    "91227"     "7161"      "221895"    "3772"     
    ## [673] "7266"      "5801"      "1894"      "8372"      "100128191" "440176"   
    ## [679] "101929530" "101928438" "728554"    "647415"    "642938"    "11162"    
    ## [685] "57379"     "134549"    "5881"      "29893"     "57001"     "285605"   
    ## [691] "29086"     "8076"      "631"       "57037"     "140686"    "1124"     
    ## [697] "26085"     "6427"      "4697"      "677802"    "25941"     "26276"    
    ## [703] "85417"     "110806273" "90355"     "83448"     "10907"     "55204"    
    ## [709] "6647"      "23169"     "51110"     "54951"     "93183"     "360132"   
    ## [715] "81931"     "55335"     "644820"    "57116"     "222235"    "9465"     
    ## [721] "10189"     "27101"     "79583"     "57537"     "79694"     "539"      
    ## [727] "7691"      "27338"     "284716"    "400750"    "84899"     "84733"    
    ## [733] "56913"     "10051"     "9958"      "10098"     "7504"      "388685"   
    ## [739] "105375650" "22998"     "92689"     "100287284" "5682"      "128308"   
    ## [745] "9519"      "5820"      "653857"    "3507"      "11315"     "54558"    
    ## [751] "811"       "4739"      "64116"     "79641"     "857"       "65999"    
    ## [757] "1292"      "2528"      "3336"      "51142"     "51643"     "100996404"
    ## [763] "57482"     "145226"    "126375"    "5168"      "57094"     "440145"   
    ## [769] "7184"      "644464"    "83871"     "28978"     "80114"     "100507650"
    ## [775] "100129543" "79730"     "1062"      "116444"    "55851"     "347734"   
    ## [781] "57167"     "1646"      "6415"      "100289373" "2632"      "4942"     
    ## [787] "1573"      "285613"    "63926"     "2258"      "400451"    "51014"    
    ## [793] "2923"      "112703"    "8876"      "162239"    "100422885" "51760"    
    ## [799] "9407"      "2109"      "55321"     "116150"    "9823"      "57121"    
    ## [805] "7697"      "729920"    "272"       "133688"    "84191"     "6917"     
    ## [811] "440957"    "9194"      "100271332" "2132"      "1158"      "51026"    
    ## [817] "112267986" "400512"    "100507633" "23639"     "58473"     "101928063"
    ## [823] "24"        "7485"      "7517"      "8974"      "342035"    "667"      
    ## [829] "84814"     "2987"      "6584"      "10644"     "57489"     "644873"   
    ## [835] "101927379" "1645"      "57336"     "10998"     "6426"      "51186"    
    ## [841] "206412"    "10231"     "26225"     "83699"     "51375"     "59286"    
    ## [847] "55687"     "91147"     "113452"    "100422872" "126526"    "205428"   
    ## [853] "9949"      "84690"     "101929746" "10715"     "283726"    "57205"    
    ## [859] "64928"     "51226"     "51218"     "5244"      "2770"      "113612"   
    ## [865] "5300"      "1410"      "54885"     "100420505" "284184"    "1515"     
    ## [871] "105371941" "100874282" "101927189" "90523"     "55969"     "79953"

``` r
write.xlsx(M2, here::here("export_WGCNA", "M2.xlsx"))

#M8
M8 <- konnectivity_norm$M8
M8 <- M8 %>% left_join(fdata, by = "ENTREZID")
sKME_8 <- sKME %>% dplyr::select(kMEM8, ENTREZID)
M8 <- M8 %>% left_join(sKME_8, by = "ENTREZID")
M8 <- M8 %>% arrange(desc(kMEM8))
#M8 <- M8 %>% arrange(desc(kWithin.norm))
print(c(M8$ENTREZID))
```

    ##   [1] "2113"      "916"       "3932"      "9806"      "919"       "915"      
    ##   [7] "81606"     "10225"     "28639"     "9750"      "201633"    "3702"     
    ##  [13] "28638"     "9402"      "22806"     "28755"     "940"       "29909"    
    ##  [19] "917"       "11278"     "4753"      "27334"     "10663"     "3820"     
    ##  [25] "10365"     "64919"     "9235"      "50852"     "926"       "3560"     
    ##  [31] "51176"     "925"       "6095"      "28589"     "23224"     "8320"     
    ##  [37] "28692"     "6932"      "4773"      "814"       "4068"      "79652"    
    ##  [43] "28558"     "28567"     "26051"     "6967"      "128611"    "53637"    
    ##  [49] "1901"      "197358"    "81539"     "3004"      "1233"      "1236"     
    ##  [55] "5243"      "2833"      "959"       "6402"      "6504"      "84636"    
    ##  [61] "27065"     "5588"      "3001"      "100130231" "23348"     "54900"    
    ##  [67] "53347"     "154075"    "340547"    "9840"      "225"       "60468"    
    ##  [73] "3738"      "2625"      "28517"     "399"       "28559"     "9047"     
    ##  [79] "28606"     "2043"      "6983"      "3983"      "126259"    "8807"     
    ##  [85] "914"       "28576"     "3003"      "10578"     "6097"      "28673"    
    ##  [91] "28671"     "101928100" "8707"      "924"       "927"       "28568"    
    ##  [97] "51440"     "50861"     "9452"      "28684"     "256380"    "6785"     
    ## [103] "55001"     "356"       "83888"     "129804"    "8911"      "168537"   
    ## [109] "84174"     "23305"     "921"       "6297"      "9026"      "3824"     
    ## [115] "5144"      "51301"     "28674"     "54855"     "2838"      "3669"     
    ## [121] "360"       "128553"    "27086"     "102724246" "100996286" "115362"   
    ## [127] "28660"     "83988"     "2775"      "6966"      "219790"    "939"      
    ## [133] "7704"      "115352"    "7089"      "28682"     "4818"      "26053"    
    ## [139] "53405"     "28614"     "3575"      "388228"    "115650"    "6489"     
    ## [145] "284486"    "64061"     "28609"     "8821"      "145864"    "3902"     
    ## [151] "92922"     "23178"     "56967"     "8519"      "387882"    "28670"    
    ## [157] "151888"    "8313"      "28573"     "146206"    "28603"     "55450"    
    ## [163] "28611"     "3131"      "59271"     "9241"      "101929241" "4538"     
    ## [169] "100506169" "8728"      "54039"     "116442"    "29121"     "28815"    
    ## [175] "4512"      "4539"      "107075141" "2149"      "7220"      "28662"    
    ## [181] "64425"     "6304"      "22997"     "28596"     "84525"     "1380"     
    ## [187] "10316"     "167681"    "647121"    "4540"      "124909475" "3502"     
    ## [193] "4541"      "8743"      "283149"    "4514"      "2841"      "3706"     
    ## [199] "105370525" "4513"      "100996455" "28586"     "3710"      "254559"   
    ## [205] "83700"     "106480796" "2039"      "4508"      "107075270" "26157"    
    ## [211] "2533"      "4603"      "4535"      "8745"      "644353"    "10016"    
    ## [217] "106478943" "101927837" "28669"     "54149"     "5174"      "3620"     
    ## [223] "54674"     "100887749" "3848"      "53344"     "107075310" "3212"     
    ## [229] "3500"      "27177"     "5473"      "11168"     "101241892" "23314"    
    ## [235] "54762"     "4892"      "6653"      "55303"     "4537"      "100271358"
    ## [241] "56114"     "10020"     "101929469" "84984"     "51305"     "101927480"
    ## [247] "643911"    "131450"    "375593"    "158471"    "106480310" "9254"     
    ## [253] "101927233" "149297"    "387036"    "79673"     "29902"     "84689"    
    ## [259] "102724368" "387644"    "390598"    "5327"      "89890"     "25845"    
    ## [265] "1102"      "79006"     "91653"     "131566"    "5874"      "9685"     
    ## [271] "2729"      "100506697" "6894"      "101059953" "79838"     "102723167"
    ## [277] "158297"    "26278"     "100420899" "54958"     "56987"     "10207"    
    ## [283] "3493"      "4199"      "340267"    "255061"    "7490"      "692093"   
    ## [289] "83999"     "60673"

``` r
write.xlsx(M8, here::here("export_WGCNA","M8.xlsx"))

#M10
M10 <- konnectivity_norm$M10
M10 <- M10 %>% left_join(fdata, by = "ENTREZID")
sKME_10 <- sKME %>% dplyr::select(kMEM10, ENTREZID)
M10 <- M10 %>% left_join(sKME_10, by = "ENTREZID")
M10 <- M10 %>% arrange(desc(kMEM10))
#M10 <- M10 %>% arrange(desc(kWithin.norm))
print(c(M10$ENTREZID))
```

    ##   [1] "23764"     "23221"     "101928707" "5603"      "100131244" "126321"   
    ##   [7] "467"       "9754"      "1605"      "84440"     "729359"    "4145"     
    ##  [13] "1822"      "10362"     "116984"    "152519"    "9372"      "7552"     
    ##  [19] "151"       "2026"      "1844"      "85019"     "3156"      "1959"     
    ##  [25] "2268"      "2539"      "219833"    "54509"     "387590"    "23457"    
    ##  [31] "1795"      "1263"      "28231"     "197187"    "27242"     "92345"    
    ##  [37] "64411"     "8525"      "91683"     "1847"      "84444"     "28984"    
    ##  [43] "3744"      "100874214" "311"       "79778"     "10758"     "7040"     
    ##  [49] "8491"      "10985"     "27"        "5990"      "101927034" "649"      
    ##  [55] "163702"    "26119"     "80176"     "2669"      "3101"      "388939"   
    ##  [61] "91252"     "8650"      "51564"     "5091"      "8740"      "400793"   
    ##  [67] "2534"      "1889"      "57184"     "253982"    "4861"      "167465"   
    ##  [73] "150290"    "6483"      "6510"      "5792"      "9170"      "56926"    
    ##  [79] "23365"     "374378"    "57572"     "5777"      "139735"    "9242"     
    ##  [85] "51606"     "11131"     "641649"    "56920"     "388121"    "100507599"
    ##  [91] "7016"      "10908"     "8913"      "3202"      "22853"     "3652"     
    ##  [97] "84536"     "29109"     "390627"    "23141"     "9531"      "9816"     
    ## [103] "3737"      "5330"      "10019"     "22876"     "57553"     "642311"   
    ## [109] "5606"      "431705"    "64773"     "64924"     "9722"      "2194"     
    ## [115] "23529"     "415116"    "49854"     "339983"    "23095"     "115804232"
    ## [121] "23235"     "400960"    "100289187" "54145"     "6197"      "154807"   
    ## [127] "56965"     "375057"    "339665"    "6713"      "29920"     "55604"    
    ## [133] "84798"     "65985"     "1832"      "3157"      "126661"    "10116"    
    ## [139] "84612"     "122616"    "27042"     "23550"     "2036"      "170850"   
    ## [145] "4628"      "165829"    "203286"    "80205"     "22993"     "54663"    
    ## [151] "9912"      "26025"     "10411"     "55703"     "84101"     "10079"    
    ## [157] "401074"    "113655"    "11103"     "57719"     "23368"     "728121"   
    ## [163] "79751"     "256364"    "79906"     "8412"      "65997"     "5017"     
    ## [169] "4598"      "10144"     "5467"      "249"       "85236"     "4953"     
    ## [175] "7450"      "57546"     "160518"    "54903"     "1911"      "10076"    
    ## [181] "7378"      "401491"    "79864"     "6307"      "29951"     "727914"   
    ## [187] "116412"    "4248"      "387066"    "11161"     "286256"    "54872"    
    ## [193] "23005"     "9181"      "112577465" "219539"    "114327"    "7832"     
    ## [199] "55582"     "6441"      "55679"     "56165"     "5287"      "29775"    
    ## [205] "283209"    "10257"     "113828"    "5833"      "79639"     "5335"     
    ## [211] "101928303" "8884"      "10531"     "6309"      "149175"    "55509"    
    ## [217] "60490"     "56180"     "79759"     "5786"      "100506207" "60676"    
    ## [223] "79029"     "114792"    "29953"     "3422"      "103344932" "1802"     
    ## [229] "3949"      "220136"    "126917"    "92999"     "64799"     "51330"    
    ## [235] "9508"      "51144"     "55283"     "84803"     "2063"      "128153"   
    ## [241] "399774"    "9051"      "492"       "8987"      "10693"     "51059"    
    ## [247] "2557"      "55262"     "51172"     "729288"    "400224"    "6576"     
    ## [253] "387496"    "3177"      "338785"    "9540"      "100996732"

``` r
write.xlsx(M10, here::here("export_WGCNA", "M10.xlsx"))

#M12
M12 <- konnectivity_norm$M12
M12 <- M12 %>% left_join(fdata, by = "ENTREZID")
sKME_12 <- sKME %>% dplyr::select(kMEM12, ENTREZID)
M12 <- M12 %>% left_join(sKME_12, by = "ENTREZID")
M12 <- M12 %>% arrange(desc(kMEM12))
#M12 <- M12 %>% arrange(desc(kWithin.norm))
print(c(M12$ENTREZID))
```

    ##   [1] "55422"     "5376"      "972"       "10855"     "29760"     "10123"    
    ##   [7] "23092"     "54704"     "4261"      "2885"      "254065"    "7454"     
    ##  [13] "9445"      "8905"      "64127"     "19"        "6310"      "2123"     
    ##  [19] "2212"      "1794"      "943"       "9060"      "57162"     "1393"     
    ##  [25] "1439"      "6001"      "51646"     "752"       "9586"      "7100"     
    ##  [31] "79891"     "25840"     "29116"     "10924"     "7132"      "4810"     
    ##  [37] "55799"     "2535"      "217"       "3188"      "57095"     "153222"   
    ##  [43] "259230"    "8495"      "283537"    "11213"     "55092"     "25900"    
    ##  [49] "653653"    "5641"      "7538"      "1152"      "359948"    "9802"     
    ##  [55] "26301"     "79974"     "5621"      "54434"     "941"       "29062"    
    ##  [61] "374655"    "127544"    "144110"    "128077"    "55812"     "102724307"
    ##  [67] "54838"     "4811"      "27239"     "140738"    "1230"      "433"      
    ##  [73] "5205"      "719"       "9936"      "5738"      "2581"      "123"      
    ##  [79] "3028"      "2207"      "2829"      "27185"     "9873"      "10608"    
    ##  [85] "9902"      "51284"     "3082"      "10766"     "10509"     "1389"     
    ##  [91] "2621"      "150275"    "27235"     "2867"      "1522"      "170575"   
    ##  [97] "286144"    "50619"     "55313"     "7127"      "168544"    "64581"    
    ## [103] "6038"      "100506388" "57486"     "119032"    "6430"      "54796"    
    ## [109] "27071"     "84967"     "84939"     "653583"    "23670"     "51313"    
    ## [115] "256586"    "3587"      "200728"    "836"       "80162"     "340152"   
    ## [121] "5166"      "144402"    "10957"     "151195"    "27202"     "51744"    
    ## [127] "2139"      "54165"     "10577"     "122953"    "55106"     "29103"    
    ## [133] "124540"    "1107"      "100505881" "100"       "81029"     "54777"    
    ## [139] "7539"      "64386"     "64798"     "79734"     "2876"      "7305"     
    ## [145] "79772"     "4678"      "9770"      "84448"     "283284"    "51444"    
    ## [151] "83478"     "10500"     "3418"      "53373"     "338339"    "162387"   
    ## [157] "4835"      "100507639" "10807"     "101927143" "2357"      "216"      
    ## [163] "57026"     "155368"    "120892"    "1536"      "684"       "54557"    
    ## [169] "3127"      "220032"    "55207"     "51601"     "53834"     "28232"    
    ## [175] "374872"    "100379592" "6275"      "148709"    "7006"      "100131465"
    ## [181] "10628"     "51421"     "344887"    "25847"     "6668"      "1486"     
    ## [187] "692205"    "326276"    "23142"     "105374758" "51259"     "475"      
    ## [193] "84224"     "23432"     "55884"     "90390"     "9398"      "9901"     
    ## [199] "646532"    "7570"      "91947"     "57136"     "55213"     "10335"    
    ## [205] "8992"      "642946"    "27345"     "56994"     "9159"      "1793"     
    ## [211] "10354"     "100506810" "692225"    "3662"      "203522"    "1356"

``` r
write.xlsx(M12, here::here("export_WGCNA","M12.xlsx"))
```

# Transcriptions factor network inference (<https://bioconductor.org/packages/devel/bioc/vignettes/decoupleR/inst/doc/tf_bk.html#1_Loading_packages>)

### Wrangle the data and preparation for input

``` r
metadata <- readRDS(here::here("data","metadata.rds"))

#provide log transformed count matrix with SYMBOL as rownames and samples as columns
x_lcpm <- cpm(x, log = TRUE)

genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)
genes$ENSEMBL <- str_replace(genes$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
rownames(genes) <- genes$ENSEMBL

genes_new <- genes
genes_new$SYMBOL <- mapIds(org.Hs.eg.db, keys = genes_new$ENSEMBL, keytype="ENSEMBL", column = "SYMBOL")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
genes_new <- genes_new %>% drop_na(SYMBOL)

#random sampling of duplicate entries
genes_new <- genes_new %>%
  group_by(SYMBOL) %>%
  sample_n(1) %>% as.data.frame()


rownames(genes_new) <- genes_new$SYMBOL
genes_new <- genes_new %>% dplyr::select(-SYMBOL, -ENSEMBL)

metadata_tf <- metadata
colnames(metadata_tf) <- c("sample" ,"Group", "Treatment", "Index", "Index_2")

genes_new_all <- genes_new
genes_new_mock <- genes_new
genes_new_beta <- genes_new

#set DEG in limma to all with p.adj_value < 0.05 to get holistic TF-view
Adult_Effect = topTable(fitDupCor, coef="Adult_Effect", number = Inf)
Preterm_Effect = topTable(fitDupCor, coef="Preterm_Effect", number = Inf)
Term_Effect = topTable(fitDupCor, coef="Term_Effect", number = Inf)
MOCK_AdTe = topTable(fitDupCor, coef="MOCK_AdTe", number = Inf)
MOCK_AdPre = topTable(fitDupCor, coef="MOCK_AdPre", number = Inf)
MOCK_PreTe = topTable(fitDupCor, coef="MOCK_PreTe", number = Inf)
GLUC_AdTe = topTable(fitDupCor, coef="GLUC_AdTe", number = Inf)
GLUC_AdPre = topTable(fitDupCor, coef="GLUC_AdPre", number = Inf)
GLUC_PreTe = topTable(fitDupCor, coef="GLUC_PreTe", number = Inf)


Adult_Effect_noFC = subset(Adult_Effect, adj.P.Val < 0.05)
Preterm_Effect_noFC = subset(Preterm_Effect, adj.P.Val < 0.05)
Term_Effect_noFC = subset(Term_Effect, adj.P.Val < 0.05)
MOCK_AdTe_noFC = subset(MOCK_AdTe, adj.P.Val < 0.05)
MOCK_AdPre_noFC = subset(MOCK_AdPre, adj.P.Val < 0.05)
MOCK_PreTe_noFC = subset(MOCK_PreTe, adj.P.Val < 0.05)
GLUC_AdTe_noFC = subset(GLUC_AdTe, adj.P.Val < 0.05)
GLUC_AdPre_noFC = subset(GLUC_AdPre, adj.P.Val < 0.05)
GLUC_PreTe_noFC = subset(GLUC_PreTe, adj.P.Val < 0.05)

#prepare DEGs for TF

#adults
deg_adult_effect <- Adult_Effect %>% rownames_to_column("ENSEMBL")
deg_adult_effect$ENSEMBL <- str_replace(deg_adult_effect$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
deg_adult_effect$SYMBOL <- mapIds(org.Hs.eg.db, keys = deg_adult_effect$ENSEMBL, keytype="ENSEMBL", column = "SYMBOL")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
deg_adult_effect <- deg_adult_effect %>% drop_na(SYMBOL)
#random sampling of duplicate entries
deg_adult_effect <- deg_adult_effect %>%
  group_by(SYMBOL) %>%
  sample_n(1) %>% as.data.frame()
rownames(deg_adult_effect) <- deg_adult_effect$SYMBOL
deg_adult_effect_all <- deg_adult_effect %>% dplyr::select(-SYMBOL, -ENSEMBL)
deg_adult_effect <- deg_adult_effect_all %>% dplyr::select(logFC, t, P.Value) %>% dplyr::filter(!is.na(t)) %>% as.matrix()


#preterm beta vs. adult_beta
deg_GLUC_preterm_effect <- GLUC_AdPre %>% rownames_to_column("ENSEMBL")
deg_GLUC_preterm_effect$ENSEMBL <- str_replace(deg_GLUC_preterm_effect$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
deg_GLUC_preterm_effect$SYMBOL <- mapIds(org.Hs.eg.db, keys = deg_GLUC_preterm_effect$ENSEMBL, keytype="ENSEMBL", column = "SYMBOL")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
deg_GLUC_preterm_effect <- deg_GLUC_preterm_effect %>% drop_na(SYMBOL)
#random sampling of duplicate entries
deg_GLUC_preterm_effect <- deg_GLUC_preterm_effect %>%
  group_by(SYMBOL) %>%
  sample_n(1) %>% as.data.frame()
rownames(deg_GLUC_preterm_effect) <- deg_GLUC_preterm_effect$SYMBOL
deg_GLUC_preterm_effect_all <- deg_GLUC_preterm_effect %>% dplyr::select(-SYMBOL, -ENSEMBL)
deg_GLUC_preterm_effect <- deg_GLUC_preterm_effect_all %>% dplyr::select(logFC, t, P.Value) %>% dplyr::filter(!is.na(t)) %>% as.matrix()


#term beta vs. adult beta
deg_GLUC_term_effect <- GLUC_AdTe %>% rownames_to_column("ENSEMBL")
deg_GLUC_term_effect$ENSEMBL <- str_replace(deg_GLUC_term_effect$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
deg_GLUC_term_effect$SYMBOL <- mapIds(org.Hs.eg.db, keys = deg_GLUC_term_effect$ENSEMBL, keytype="ENSEMBL", column = "SYMBOL")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
deg_GLUC_term_effect <- deg_GLUC_term_effect %>% drop_na(SYMBOL)
#random sampling of duplicate entries
deg_GLUC_term_effect <- deg_GLUC_term_effect %>%
  group_by(SYMBOL) %>%
  sample_n(1) %>% as.data.frame()
rownames(deg_GLUC_term_effect) <- deg_GLUC_term_effect$SYMBOL
deg_GLUC_term_effect_all <- deg_GLUC_term_effect %>% dplyr::select(-SYMBOL, -ENSEMBL)
deg_GLUC_term_effect <- deg_GLUC_term_effect_all %>% dplyr::select(logFC, t, P.Value) %>% dplyr::filter(!is.na(t)) %>% as.matrix()
```

### Heatmap of transcription factors

``` r
#get TF-interaction network
net <- get_collectri(organism='human', split_complexes=FALSE)

#calculate permutation test for all samples
#very intensive -> saved
#sample_acts_all <- run_wmean(mat=genes_new, net=net, .source='source', .target='target',
                         #.mor='mor', times = 10000, minsize = 5)
sample_acts <- rio::import(here::here("data", "sample_acts_all.rds"))

#select TFs to display in heatmap later on
n_tfs =150

#select TFs and get most variable ones
sample_acts_mat <- sample_acts %>%
  dplyr::filter(statistic == 'corr_wmean') %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

tfs <- sample_acts %>%
  group_by(source) %>%
  summarise(std = sd(score)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

#wrangle + get metadata
sample_acts_mat <- sample_acts_mat[,tfs] %>% as.data.frame()
sample_acts_mat$SampleName <- rownames(sample_acts_mat)
sample_acts_mat <- sample_acts_mat %>% left_join(metadata, by = "SampleName")

#transpose and scale selected number of samples for heatmap + define coldata etc.
sample_acts_mat_scaled <- scale(sample_acts_mat[1:n_tfs]) %>% (t)
sample_acts_mat_metadata <- sample_acts_mat %>% dplyr::select(Group, Treatment)

col_an = HeatmapAnnotation(Stimulation = sample_acts_mat_metadata$Treatment, Group = sample_acts_mat_metadata$Group, col = list(Group = c("adult" = "#b2c5b3", "term" = "#8dafd1", "preterm" = "#ce8793"), Stimulation = c("betaglucan" = "#4D3C7E", "mock" = "#DEC08B")))

my_palette <-  colorRampPalette(c("#0765A0", "#FFF1E4","#C61700"))(100)

#Visualise in heatmap
Heatmap(sample_acts_mat_scaled, show_row_names = FALSE, show_row_dend = TRUE, col = my_palette,
                                clustering_method_columns = "complete", clustering_method_rows = "complete",
                                column_dend_side = "top", column_dend_height = unit(4, "cm"), column_km = 3, column_gap =unit(3, "mm"),   
        column_title_gp = gpar(fontsize = 10), top_annotation = col_an)
```

![](Trained_immunity_GIT_files/figure-gfm/setup78-1.png)<!-- -->

``` r
#does not look entirely convincing



#do the same for MOCK-samples
genes_new <- t(genes_new_mock) %>% as.data.frame()
genes_new <- genes_new %>% rownames_to_column("SampleName")
genes_new <- genes_new %>% left_join(metadata, by ="SampleName")
genes_new <- genes_new %>% dplyr::filter(Treatment == "mock")
genes_new <- genes_new %>% column_to_rownames("SampleName")
genes_new <- genes_new %>% dplyr::select(-Treatment, -Group, -Index, -Index_2)
genes_new <- t(genes_new) %>% as.data.frame()

#calculate permutation test for MOCK samples
#sample_acts_mock <- run_wmean(mat=genes_new, net=net, .source='source', .target='target',
                         #.mor='mor', times = 10000, minsize = 5)

sample_acts_mock <- rio::import(here::here("data", "sample_acts_mock.rds"))

#select TFs to display in heatmap later on
n_tfs = 150

#select TFs and get most variable ones
sample_acts_mat_mock <- sample_acts_mock %>%
  dplyr::filter(statistic == 'corr_wmean') %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

tfs <- sample_acts_mock %>%
  group_by(source) %>%
  summarise(iqr = iqr(score)) %>%
  arrange(-abs(iqr)) %>%
  head(n_tfs) %>%
  pull(source)

#wrangle + get metadata
sample_acts_mat_mock <- sample_acts_mat_mock[,tfs] %>% as.data.frame()
sample_acts_mat_mock$SampleName <- rownames(sample_acts_mat_mock)
sample_acts_mat_mock <- sample_acts_mat_mock %>% left_join(metadata, by = "SampleName")


#transpose and scale selected number of samples for heatmap + define coldata etc.
sample_acts_mat_scaled_mock <- scale(sample_acts_mat_mock[1:n_tfs]) %>% (t)
sample_acts_mat_metadata_mock <- sample_acts_mat_mock %>% dplyr::select(Group, Treatment)

col_an_mock = HeatmapAnnotation( Group = sample_acts_mat_metadata_mock$Group, col = list(Group = c("adult" = "#b2c5b3", "term" = "#8dafd1", "preterm" = "#ce8793")))

my_palette <-  colorRampPalette(c("#0765A0", "#FFF1E4","#C61700"))(100)

#Visualise in heatmap
Heatmap(sample_acts_mat_scaled_mock, show_row_names = FALSE, show_row_dend = TRUE, col = my_palette,
                                clustering_method_columns = "complete", clustering_method_rows = "complete",
                                column_dend_side = "top", column_dend_height = unit(4, "cm"), column_km = 3, column_gap =unit(3, "mm"),   
        column_title_gp = gpar(fontsize = 10), top_annotation = col_an_mock)
```

![](Trained_immunity_GIT_files/figure-gfm/setup78-2.png)<!-- -->

``` r
#also does not look entirely convincing



#do the same for beta-Glucan samples
genes_new <- t(genes_new_beta) %>% as.data.frame()
genes_new <- genes_new %>% rownames_to_column("SampleName")
genes_new <- genes_new %>% left_join(metadata, by ="SampleName")
genes_new <- genes_new %>% dplyr::filter(Treatment == "betaglucan")
genes_new <- genes_new %>% column_to_rownames("SampleName")
genes_new <- genes_new %>% dplyr::select(-Treatment, -Group, -Index, -Index_2)
genes_new <- t(genes_new) %>% as.data.frame()

#calculate permutation test for Beta-Glucan samples
#sample_acts_beta <- run_wmean(mat=genes_new, net=net, .source='source', .target='target',
                        # .mor='mor', times = 10000, minsize = 5)

sample_acts_beta <- rio::import(here::here("data", "sample_acts_beta.rds"))

#select TFs to display in heatmap later on
n_tfs = 150

#select TFs and get most variable ones
sample_acts_mat_beta <- sample_acts_beta %>%
  dplyr::filter(statistic == 'corr_wmean') %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

tfs <- sample_acts_beta %>%
  group_by(source) %>%
  summarise(iqr = iqr(score)) %>%
  arrange(-abs(iqr)) %>%
  head(n_tfs) %>%
  pull(source)

#wrangle + get metadata
sample_acts_mat_beta <- sample_acts_mat_beta[,tfs] %>% as.data.frame()
sample_acts_mat_beta$SampleName <- rownames(sample_acts_mat_beta)
sample_acts_mat_beta <- sample_acts_mat_beta %>% left_join(metadata, by = "SampleName")


#transpose and scale selected number of samples for heatmap + define coldata etc.
sample_acts_mat_scaled_beta <- scale(sample_acts_mat_beta[1:n_tfs]) %>% (t)
sample_acts_mat_metadata_beta <- sample_acts_mat_beta %>% dplyr::select(Group, Treatment)

col_an_beta = HeatmapAnnotation( Group = sample_acts_mat_metadata_beta$Group, col = list(Group = c("adult" = "#b2c5b3", "term" = "#8dafd1", "preterm" = "#ce8793")))

my_palette <-  colorRampPalette(c("#0765A0", "#FFF1E4","#C61700"))(100)

#Visualise in heatmap
Heatmap(sample_acts_mat_scaled_beta, show_row_names = FALSE, show_row_dend = TRUE, col = my_palette,
        clustering_method_columns = "complete", clustering_method_rows = "complete",
        column_dend_side = "top", column_dend_height = unit(4, "cm"), column_km = 3, column_gap =unit(3, "mm"),   
        column_title_gp = gpar(fontsize = 10), top_annotation = col_an_beta)
```

![](Trained_immunity_GIT_files/figure-gfm/setup78-3.png)<!-- -->

``` r
#also here there is quite high donor-variability --> another approach needed
```

## TF-activity inferece

``` r
#for adult samples (change over treatment in adults)
#contrast_acts_ad <- run_wmean(mat=deg_adult_effect[, 't', drop=FALSE], net=net, .source='source', .target='target',
                           #.mor='mor', times = 10000, minsize = 5)

contrast_acts_ad <- rio::import(here::here("data", "contrast_acts_ad.rds"))

#get values
f_contrast_acts_ad <- contrast_acts_ad %>%
  dplyr::filter(statistic == 'corr_wmean') %>%
  mutate(rnk = NA)

msk <- f_contrast_acts_ad$score > 0
f_contrast_acts_ad[msk, 'rnk'] <- rank(-f_contrast_acts_ad[msk, 'score'])
f_contrast_acts_ad[!msk, 'rnk'] <- rank(-abs(f_contrast_acts_ad[!msk, 'score']))

#select number of Tfs
n_tfs <- 20

tfs_ad <- f_contrast_acts_ad %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)
f_contrast_acts_ad <- f_contrast_acts_ad %>%
  dplyr::filter(source %in% tfs_ad)

#visualize
tf_enrichment_adult_effect <- ggplot(f_contrast_acts_ad, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "#0765A0", high = "#C61700", 
                       mid = "#FFF1E4", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 90, hjust = 1, size =7),
        axis.text.y = element_text(size =7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

#export
graph2svg(tf_enrichment_adult_effect, file =  here::here("plots","tf_enrichment_adult_effect"), width = 3.5, height = 2.5)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/tf_enrichment_adult_effect.svg

``` r
#comparison between term and adult beta-glucan treated samples
#contrast_acts_term_ad <- run_wmean(mat=deg_GLUC_term_effect[, 't', drop=FALSE], net=net, .source='source', .target='target',
                           #.mor='mor', times = 10000, minsize = 5)

contrast_acts_term_ad <- rio::import(here::here("data", "contrast_acts_term_ad.rds"))

#get values
f_contrast_acts_term_ad <- contrast_acts_term_ad %>%
  dplyr::filter(statistic == 'corr_wmean') %>%
  mutate(rnk = NA)

msk <- f_contrast_acts_term_ad$score > 0
f_contrast_acts_term_ad[msk, 'rnk'] <- rank(-f_contrast_acts_term_ad[msk, 'score'])
f_contrast_acts_term_ad[!msk, 'rnk'] <- rank(-abs(f_contrast_acts_term_ad[!msk, 'score']))

#select number of Tfs
n_tfs <-20

tfs_term_ad <- f_contrast_acts_term_ad %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)
f_contrast_acts_term_ad <- f_contrast_acts_term_ad %>%
  dplyr::filter(source %in% tfs_term_ad)

#visualize
tf_enrichment_term_ad_BETA <- ggplot(f_contrast_acts_term_ad, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "#0765A0", high = "#C61700", 
                       mid = "#FFF1E4", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 90, hjust = 1, size =7),
        axis.text.y = element_text(size =7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

#export
graph2svg(tf_enrichment_term_ad_BETA, file =  here::here("plots","tf_enrichment_term_ad_BETA"), width = 3.5, height = 2.5)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/tf_enrichment_term_ad_BETA.svg

``` r
#comparison between preterm and adult beta-glucan treated samples
#contrast_acts_pre_ad <- run_wmean(mat=deg_GLUC_preterm_effect[, 't', drop=FALSE], net=net, .source='source', .target='target',
                                  #.mor='mor', times = 10000, minsize = 5)

contrast_acts_pre_ad <- rio::import(here::here("data", "contrast_acts_pre_ad.rds"))

#get values
f_contrast_acts_pre_ad <- contrast_acts_pre_ad %>%
  dplyr::filter(statistic == 'corr_wmean') %>%
  mutate(rnk = NA)

msk <- f_contrast_acts_pre_ad$score > 0
f_contrast_acts_pre_ad[msk, 'rnk'] <- rank(-f_contrast_acts_pre_ad[msk, 'score'])
f_contrast_acts_pre_ad[!msk, 'rnk'] <- rank(-abs(f_contrast_acts_pre_ad[!msk, 'score']))

#select number of Tfs
n_tfs <- 20

tfs_pre_ad <- f_contrast_acts_pre_ad %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)
f_contrast_acts_pre_ad <- f_contrast_acts_pre_ad %>%
  dplyr::filter(source %in% tfs_pre_ad)

#visualize
tf_enrichment_preterm_ad_BETA <- ggplot(f_contrast_acts_pre_ad, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "#0765A0", high = "#C61700", 
                       mid = "#FFF1E4", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 90, hjust = 1, size =7),
        axis.text.y = element_text(size =7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

#export
graph2svg(tf_enrichment_preterm_ad_BETA, file = here::here("plots","tf_enrichment_preterm_ad_BETA"), width = 3.5, height = 2.5)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/tf_enrichment_preterm_ad_BETA.svg

## Reverse mapping of Transcription factor activity

``` r
#prepare cpm-transformed gene counts with gene-annotations
x_lcpm <- cpm(x)
genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)
genes$ENSEMBL <- str_replace(genes$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
rownames(genes) <- genes$ENSEMBL
genes_new <- genes
genes_new$SYMBOL <- mapIds(org.Hs.eg.db, keys = genes_new$ENSEMBL, keytype="ENSEMBL", column = "SYMBOL")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
genes_new <- genes_new %>% drop_na(SYMBOL)

#random sampling of duplicate entries
genes_new <- genes_new %>%
  group_by(SYMBOL) %>%
  sample_n(1) %>% as.data.frame()

rownames(genes_new) <- genes_new$SYMBOL
genes_new <- genes_new %>% dplyr::select(-ENSEMBL)

#more lenient filter settings for subsequent exploratory analyses
Adult_Effect = topTable(fitDupCor, coef="Adult_Effect", number = Inf)
Preterm_Effect = topTable(fitDupCor, coef="Preterm_Effect", number = Inf)
Term_Effect = topTable(fitDupCor, coef="Term_Effect", number = Inf)
MOCK_AdTe = topTable(fitDupCor, coef="MOCK_AdTe", number = Inf)
MOCK_AdPre = topTable(fitDupCor, coef="MOCK_AdPre", number = Inf)
MOCK_PreTe = topTable(fitDupCor, coef="MOCK_PreTe", number = Inf)
GLUC_AdTe = topTable(fitDupCor, coef="GLUC_AdTe", number = Inf)
GLUC_AdPre = topTable(fitDupCor, coef="GLUC_AdPre", number = Inf)
GLUC_PreTe = topTable(fitDupCor, coef="GLUC_PreTe", number = Inf)


Adult_Effect_noFC = subset(Adult_Effect, adj.P.Val < 0.1)
Preterm_Effect_noFC = subset(Preterm_Effect, adj.P.Val < 0.1)
Term_Effect_noFC = subset(Term_Effect, adj.P.Val < 0.1)
MOCK_AdTe_noFC = subset(MOCK_AdTe, adj.P.Val < 0.1)
MOCK_AdPre_noFC = subset(MOCK_AdPre, adj.P.Val < 0.1)
MOCK_PreTe_noFC = subset(MOCK_PreTe, adj.P.Val < 0.1)
GLUC_AdTe_noFC = subset(GLUC_AdTe, adj.P.Val < 0.1)
GLUC_AdPre_noFC = subset(GLUC_AdPre, adj.P.Val < 0.1)
GLUC_PreTe_noFC = subset(GLUC_PreTe, adj.P.Val < 0.1)

#prepare vectors to subset the data

#adults
vec_adult_effect <- Adult_Effect_noFC %>% rownames_to_column("ENSEMBL")
vec_adult_effect$ENSEMBL <- str_replace(vec_adult_effect$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
vec_adult_effect$ENTREZID <- mapIds(org.Hs.eg.db, keys = vec_adult_effect$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
vec_adult_effect <- vec_adult_effect %>% drop_na(ENTREZID)
#random sampling of duplicate entries
vec_adult_effect <- vec_adult_effect %>%
  group_by(ENTREZID) %>%
  sample_n(1) %>% as.data.frame()
vec_adult_effect <- vec_adult_effect$ENTREZID %>% unique()

#term vs. adult beta
vec_adult_term_beta_effect <- GLUC_AdTe_noFC %>% rownames_to_column("ENSEMBL")
vec_adult_term_beta_effect$ENSEMBL <- str_replace(vec_adult_term_beta_effect$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
vec_adult_term_beta_effect$ENTREZID <- mapIds(org.Hs.eg.db, keys = vec_adult_term_beta_effect$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
vec_adult_term_beta_effect <- vec_adult_term_beta_effect %>% drop_na(ENTREZID)
#random sampling of duplicate entries
vec_adult_term_beta_effect <- vec_adult_term_beta_effect %>%
  group_by(ENTREZID) %>%
  sample_n(1) %>% as.data.frame()
vec_adult_term_beta_effect <- vec_adult_term_beta_effect$ENTREZID %>% unique()

#preterm vs. adult beta
vec_adult_preterm_beta_effect <- GLUC_AdPre_noFC %>% rownames_to_column("ENSEMBL")
vec_adult_preterm_beta_effect$ENSEMBL <- str_replace(vec_adult_preterm_beta_effect$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
vec_adult_preterm_beta_effect$ENTREZID <- mapIds(org.Hs.eg.db, keys = vec_adult_preterm_beta_effect$ENSEMBL, keytype="ENSEMBL", column = "ENTREZID")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
vec_adult_preterm_beta_effect <- vec_adult_preterm_beta_effect %>% drop_na(ENTREZID)
#random sampling of duplicate entries
vec_adult_preterm_beta_effect <- vec_adult_preterm_beta_effect %>%
  group_by(ENTREZID) %>%
  sample_n(1) %>% as.data.frame()
vec_adult_preterm_beta_effect <- vec_adult_preterm_beta_effect$ENTREZID %>% unique()

#import TFlink data (all, homo sapiens, downloaded on 11.01.24)
TF_link <- rio::import(here::here("data","TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv"))
```

## Analysis of enriched trancription factors in adults upon beta-glucan treatment

``` r
#get TF-link data for reverse mapping (is also a more extensive database = more suitable for exploratory analyses)
TF_link <- TF_link %>% dplyr::filter(Organism == "Homo sapiens")
TF_Link <- TF_link %>% dplyr::select(Name.TF, Name.Target)

#use top 3/bot 3 TFs from last analyses (most up/downregulated)
TF_link_adults_FOXC2 <-  subset(TF_link, (TF_link$Name.TF %in% c("FOXC2"))) #only 2 entries, excluded from subsequent analyses
TF_link_adults_ARNT <-  subset(TF_link, (TF_link$Name.TF %in% c("ARNT")))
TF_link_adults_BACH2 <-  subset(TF_link, (TF_link$Name.TF %in% c("BACH2")))
TF_link_adults_SMAD5 <- subset(TF_link, (TF_link$Name.TF %in% c("SMAD5"))) 
TF_link_adults_PANK7 <-  subset(TF_link, (TF_link$Name.TF %in% c("PANK7"))) #only 6 entries, excluded from subsequent analyses
TF_link_adults_ZNF300 <-  subset(TF_link, (TF_link$Name.TF %in% c("ZNF300"))) #only 6 entries, excluded from subsequent analyses


#ARNT
TF_link_adults_ARNT_GENEID <- TF_link_adults_ARNT %>% dplyr::select(NCBI.GeneID.TF, NCBI.GeneID.Target)
Adults_ARNT_GENEID_TF <- TF_link_adults_ARNT_GENEID$NCBI.GeneID.TF %>% unique() %>% as.data.frame()
Adults_ARNT_GENEID_Target <- TF_link_adults_ARNT_GENEID$NCBI.GeneID.Target %>% unique() %>% as.data.frame()
colnames(Adults_ARNT_GENEID_Target) <- c("NCBI.GeneID.Target")
Adults_ARNT_GENEID_Target <- subset(Adults_ARNT_GENEID_Target, (Adults_ARNT_GENEID_Target$NCBI.GeneID.Target %in% vec_adult_effect))
vec_Adults_ARNT_GENEID_Target <-  Adults_ARNT_GENEID_Target$NCBI.GeneID.Target  %>% unique()

#SMAD5
TF_link_adults_SMAD5_GENEID <- TF_link_adults_SMAD5 %>% dplyr::select(NCBI.GeneID.TF, NCBI.GeneID.Target)
Adults_SMAD5_GENEID_TF <- TF_link_adults_SMAD5_GENEID$NCBI.GeneID.TF %>% unique() %>% as.data.frame()
Adults_SMAD5_GENEID_Target <- TF_link_adults_SMAD5_GENEID$NCBI.GeneID.Target %>% unique() %>% as.data.frame()
colnames(Adults_SMAD5_GENEID_Target) <- c("NCBI.GeneID.Target")
Adults_SMAD5_GENEID_Target <- subset(Adults_SMAD5_GENEID_Target, (Adults_SMAD5_GENEID_Target$NCBI.GeneID.Target %in% vec_adult_effect))
vec_Adults_SMAD5_GENEID_Target <-  Adults_SMAD5_GENEID_Target$NCBI.GeneID.Target  %>% unique()

#BACH2
TF_link_adults_BACH2_GENEID <- TF_link_adults_BACH2 %>% dplyr::select(NCBI.GeneID.TF, NCBI.GeneID.Target)
Adults_BACH2_GENEID_TF <- TF_link_adults_BACH2_GENEID$NCBI.GeneID.TF %>% unique() %>% as.data.frame()
Adults_BACH2_GENEID_Target <- TF_link_adults_BACH2_GENEID$NCBI.GeneID.Target %>% unique() %>% as.data.frame()
colnames(Adults_BACH2_GENEID_Target) <- c("NCBI.GeneID.Target")
Adults_BACH2_GENEID_Target <- subset(Adults_BACH2_GENEID_Target, (Adults_BACH2_GENEID_Target$NCBI.GeneID.Target %in% vec_adult_effect))
vec_Adults_BACH2_GENEID_Target <-  Adults_BACH2_GENEID_Target$NCBI.GeneID.Target  %>% unique()


#perfrom ORA with target genes of respective transcription factors

#set background
background_genes <- readRDS(here::here("data","se.rds"))
background_genes <- background_genes@rowRanges@partitioning@NAMES %>% as.data.frame()

background_genes$ENSEMBL <- str_replace(background_genes$.,
                                        pattern = ".[0-9]+$",
                                        replacement = "")

background_genes <- background_genes %>% dplyr::select(-.)

"entrez_id" = mapIds(
  # Replace with annotation package for the organism relevant to your data
  org.Hs.eg.db,
  keys =  background_genes$ENSEMBL,
  # Replace with the type of gene identifiers in your data
  keytype = "ENSEMBL",
  # Replace with the type of gene identifiers you would like to map to
  column = "ENTREZID",
  # This will keep only the first mapped value for each Ensembl ID
  multiVals = "first"
)
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
background_genes$ENTREZ <- entrez_id
background_genes <- background_genes %>% dplyr::filter(!is.na(ENTREZ)) %>% dplyr::select(-ENSEMBL)
background_genes <- background_genes$ENTREZ %>% unique()


#loop GO:BP ORA
list_DF_adult_beta <- list(vec_Adults_ARNT_GENEID_Target = vec_Adults_ARNT_GENEID_Target, vec_Adults_SMAD5_GENEID_Target= vec_Adults_SMAD5_GENEID_Target,
                           vec_Adults_BACH2_GENEID_Target = vec_Adults_BACH2_GENEID_Target)

adult_beta_BP <- list()
for (i in 1:3) {
  enrich_adult_beta_BP <- enrichGO(gene = list_DF_adult_beta[[i]],
                                        universe = background_genes,
                                        OrgDb = org.Hs.eg.db,
                                        ont = "BP",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff = 0.1,
                                        qvalueCutoff = 0.1,
                                        readable = TRUE)
  adult_beta_BP[[i]] <- enrich_adult_beta_BP
}

#save files
#saveRDS(adult_beta_BP, "C:/Users/micha/OneDrive/Dokumente/R-Trained_Immunity/data/RDS_objects/TF_stuff/adult_beta_BP.rds")
#adult_beta_BP <- rio::import("C:/Users/micha/OneDrive/Dokumente/R-Trained_Immunity/data/RDS_objects/TF_stuff/adult_beta_BP.rds")


print(enrich_adult_beta_BP)
```

    ## #
    ## # over-representation test
    ## #
    ## #...@organism     Homo sapiens 
    ## #...@ontology     BP 
    ## #...@keytype      ENTREZID 
    ## #...@gene     chr [1:133] "5702" "54518" "9123" "54953" "2023" "1606" "23433" "7277" ...
    ## #...pvalues adjusted by 'BH' with cutoff <0.1 
    ## #...14 enriched terms found
    ## 'data.frame':    14 obs. of  9 variables:
    ##  $ ID         : chr  "GO:0043161" "GO:0042176" "GO:0061136" "GO:1901800" ...
    ##  $ Description: chr  "proteasome-mediated ubiquitin-dependent protein catabolic process" "regulation of protein catabolic process" "regulation of proteasomal protein catabolic process" "positive regulation of proteasomal protein catabolic process" ...
    ##  $ GeneRatio  : chr  "22/125" "12/125" "9/125" "7/125" ...
    ##  $ BgRatio    : chr  "456/18584" "364/18584" "196/18584" "117/18584" ...
    ##  $ pvalue     : num  3.75e-13 6.43e-06 7.16e-06 1.41e-05 3.04e-05 ...
    ##  $ p.adjust   : num  9.14e-10 5.82e-03 5.82e-03 8.61e-03 1.40e-02 ...
    ##  $ qvalue     : num  8.04e-10 5.12e-03 5.12e-03 7.58e-03 1.23e-02 ...
    ##  $ geneID     : chr  "PSMC3/CALR/SIAH2/PSMB5/PSMC5/PSMD4/VCP/ADRM1/PSMB7/PSMD1/CBFA2T3/PSMC2/PSMD3/PSMD11/PSMB4/ANAPC7/PSMC1/PSMA5/HE"| __truncated__ "PSMC3/PSMC5/NQO1/VCP/RHBDF1/PSMD1/CBFA2T3/PSMC2/PSMD3/PSMC1/HERPUD1/HSP90AB1" "PSMC3/PSMC5/VCP/RHBDF1/CBFA2T3/PSMC2/PSMC1/HERPUD1/HSP90AB1" "PSMC3/PSMC5/VCP/CBFA2T3/PSMC2/PSMC1/HERPUD1" ...
    ##  $ Count      : int  22 12 9 7 9 3 7 10 4 8 ...
    ## #...Citation
    ##  T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu.
    ##  clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.
    ##  The Innovation. 2021, 2(3):100141

``` r
#prepare data for heatmap
ARNT_BP_adult <- adult_beta_BP[[1]] %>% dplyr::filter(adult_beta_BP[[1]]@result[["p.adjust"]] <= 0.1)
enrichplot::dotplot(ARNT_BP_adult)
```

![](Trained_immunity_GIT_files/figure-gfm/setup81-1.png)<!-- -->

``` r
dotplot_ARNT <- enrichplot::dotplot(ARNT_BP_adult)

ARNT_ad <- dotplot_ARNT$data %>% as.data.frame()
ARNT_Proteasome <- ARNT_ad[c(1,2,4,5,7,10),]
ARNT_Vesicle <- ARNT_ad[c(3,9),]
ARNT_Nucleoside_metabolism <- ARNT_ad[c(8,6),]

vector_ARNT_Proteasome <- ARNT_Proteasome$geneID
vector_ARNT_Proteasome <- strsplit(vector_ARNT_Proteasome, "/") %>% unlist() %>% unique()

vector_ARNT_Vesicle <- ARNT_Vesicle$geneID
vector_ARNT_Vesicle <- strsplit(vector_ARNT_Vesicle, "/") %>% unlist() %>% unique()

vector_ARNT_Nucleoside_metabolism <- ARNT_Nucleoside_metabolism$geneID
vector_ARNT_Nucleoside_metabolism <- strsplit(vector_ARNT_Nucleoside_metabolism, "/") %>% unlist() %>% unique()

#here a filter can be set for a small heatmap!

SMAD5_BP_adult <- adult_beta_BP[[2]] %>% dplyr::filter(adult_beta_BP[[2]]@result[["p.adjust"]] <= 0.1)
enrichplot::dotplot(SMAD5_BP_adult)
```

![](Trained_immunity_GIT_files/figure-gfm/setup81-2.png)<!-- -->

``` r
dotplot_SMAD5 <- enrichplot::dotplot(SMAD5_BP_adult)

SMAD5_ad <- dotplot_SMAD5$data %>% as.data.frame()
SMAD5_ad_proteasome <- SMAD5_ad[c(1,3,4,5,6,9,10),]
SMAD5_ad_cellresp <- SMAD5_ad[c(7),]
SMAD5_ad_Vesicle <- SMAD5_ad[c(2),]
SMAD5_ad_Nucleoside_metabolism <- SMAD5_ad[c(8),]

vector_SMAD5_ad_proteasome <- SMAD5_ad_proteasome$geneID
vector_SMAD5_ad_proteasome <- strsplit(vector_SMAD5_ad_proteasome, "/") %>% unlist() %>% unique()

vector_SMAD5_ad_cell_resp <- SMAD5_ad_cellresp$geneID
vector_SMAD5_ad_cell_resp <- strsplit(vector_SMAD5_ad_cell_resp, "/") %>% unlist() %>% unique()

vector_SMAD5_ad_Vesicle <- SMAD5_ad_Vesicle$geneID
vector_SMAD5_ad_Vesicle <- strsplit(vector_SMAD5_ad_Vesicle, "/") %>% unlist() %>% unique()

vector_SMAD5_Nucleoside_metabolism <- SMAD5_ad_Nucleoside_metabolism$geneID
vector_SMAD5_Nucleoside_metabolism <- strsplit(vector_SMAD5_Nucleoside_metabolism, "/") %>% unlist() %>% unique()

#here a filter can be set for a small heatmap!
BACH2_BP_adult <- adult_beta_BP[[3]] %>% dplyr::filter(adult_beta_BP[[3]]@result[["p.adjust"]] <= 0.1)
enrichplot::dotplot(BACH2_BP_adult)
```

![](Trained_immunity_GIT_files/figure-gfm/setup81-3.png)<!-- -->

``` r
dotplot_BACH2 <- enrichplot::dotplot(BACH2_BP_adult)

BACH2_ad <- dotplot_BACH2$data %>% as.data.frame()
BACH2_ad_proteasome <- BACH2_ad[c(1,2,3,4,5,7,8),]

vector_BACH2_ad_proteasome <- BACH2_ad_proteasome$geneID
vector_BACH2_ad_proteasome <- strsplit(vector_BACH2_ad_proteasome, "/") %>% unlist() %>% unique()


#define overall vectors
Prot_ad_eff <- c(vector_ARNT_Proteasome,vector_SMAD5_ad_proteasome,vector_BACH2_ad_proteasome) %>% unique()
Ves_ad_eff <- c(vector_ARNT_Vesicle, vector_SMAD5_ad_Vesicle) %>% unique()
Resp_ad_eff <- c(vector_SMAD5_ad_cell_resp) %>% unique()
Nuc_ad_eff <- c(vector_ARNT_Nucleoside_metabolism, vector_SMAD5_Nucleoside_metabolism) %>% unique()
```

### Heatmap of treatment-related effects in adults

``` r
#get cpm-transformed count matrix
x_lcpm <- cpm(x)

genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)
genes$ENSEMBL <- str_replace(genes$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
rownames(genes) <- genes$ENSEMBL

genes_new <- genes
genes_new$SYMBOL <- mapIds(org.Hs.eg.db, keys = genes_new$ENSEMBL, keytype="ENSEMBL", column = "SYMBOL")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
genes_new <- genes_new %>% drop_na(SYMBOL)

#random sampling of duplicate entries
genes_new <- genes_new %>%
  group_by(SYMBOL) %>%
  sample_n(1) %>% as.data.frame()

#get metadata into dataframe
rownames(genes_new) <- genes_new$SYMBOL
genes_new <- t(genes_new) %>% as.data.frame()
genes_new <- genes_new %>% rownames_to_column("SampleName")
genes_new <- genes_new[-c(27:28),]
genes_new <- genes_new %>% left_join(metadata, by = "SampleName")
rownames(genes_new) <- genes_new$SampleName 
genes_new <- genes_new %>% mutate_at(c(2:14513), as.numeric)
genes_new <- genes_new %>% mutate_at(c(14514:14517), as.factor)

genes_new_adult_effect <- genes_new %>% dplyr::filter(genes_new$Group == "adult")
genes_new_adult_effect <- genes_new_adult_effect %>% pivot_longer(cols = 2:14513, values_to = "expression", names_to = "SYMBOL")

#Proteaseome related genes
genes_new_PROT <- subset(genes_new_adult_effect, (genes_new_adult_effect$SYMBOL %in% Prot_ad_eff))
genes_new_PROT$Pathway <- c("Proteasome related")
htmp_PROT_ad <- genes_new_PROT %>% pivot_wider(values_from = "expression", names_from = "SYMBOL") %>% as.data.frame()
rownames(htmp_PROT_ad) <- htmp_PROT_ad$SampleName
htmp_PROT_ad <- htmp_PROT_ad %>% dplyr::select(-SampleName)

htmp_PROT_ad_scaled <- scale(htmp_PROT_ad[,c(6:64)])
htmp_PROT_ad_metadata <- htmp_PROT_ad %>% dplyr::select(c(1:5))
htmp_PROT_ad_metadata <- htmp_PROT_ad_metadata %>% mutate_at(c(1:5), as.factor)
  
row_ha  = rowAnnotation(Treatment = htmp_PROT_ad_metadata$Treatment,col = list(Treatment = c("mock" = "#DEC08B", "betaglucan" = "#4D3C7E")))

col_fun_2 = colorRamp2(c(-4,-1.45, -0.12,0,0.12,1.45,4), c("#0765A0","#5d95b8","#FFF1E4", "#FFF1E4","#FFF1E4" ,"#de725f","#c61700"))
lgd = Legend(col_fun = col_fun_2, title = "z-score")


a <- Heatmap(htmp_PROT_ad_scaled, show_row_dend = TRUE, col = col_fun_2,
        clustering_method_columns = "complete", clustering_method_rows = "complete",
        column_dend_side = "top", column_dend_height = unit(0.5, "cm"),   
        column_title_gp = gpar(fontsize = 10), right_annotation = row_ha, column_names_gp = grid::gpar(fontsize = 5))

#Vesicle related genes
genes_new_VES <- subset(genes_new_adult_effect, (genes_new_adult_effect$SYMBOL %in% Ves_ad_eff))
genes_new_VES$Pathway <- c("Vesicle related")

htmp_VES_ad <- genes_new_VES %>% pivot_wider(values_from = "expression", names_from = "SYMBOL") %>% as.data.frame()
rownames(htmp_VES_ad) <- htmp_VES_ad$SampleName
htmp_VES_ad <- htmp_VES_ad %>% dplyr::select(-SampleName)

htmp_VES_ad_scaled <- scale(htmp_VES_ad[,c(6:43)])
htmp_VES_ad_metadata <- htmp_VES_ad %>% dplyr::select(c(1:5))
htmp_VES_ad_metadata <- htmp_VES_ad_metadata %>% mutate_at(c(1:5), as.factor)

row_ha  = rowAnnotation(Treatment = htmp_VES_ad_metadata$Treatment,col = list(Treatment = c("mock" = "#DEC08B", "betaglucan" = "#4D3C7E")))

col_fun_2 = colorRamp2(c(-4,-1.45, -0.12,0,0.12,1.45,4), c("#0765A0","#5d95b8","#FFF1E4", "#FFF1E4","#FFF1E4" ,"#de725f","#c61700"))
lgd = Legend(col_fun = col_fun_2, title = "z-score")

b <- Heatmap(htmp_VES_ad_scaled, show_row_dend = TRUE, col = my_palette,
             clustering_method_columns = "complete", clustering_method_rows = "complete",
             column_dend_side = "top", column_dend_height = unit(0.5, "cm"),   
             column_title_gp = gpar(fontsize = 10), right_annotation = row_ha, column_names_gp = grid::gpar(fontsize = 5))

#Ox-phos related genes
genes_new_RESP <- subset(genes_new_adult_effect, (genes_new_adult_effect$SYMBOL %in% Resp_ad_eff))
genes_new_RESP$Pathway <- c("Cell. respiration related")

htmp_RESP_ad <- genes_new_RESP %>% pivot_wider(values_from = "expression", names_from = "SYMBOL") %>% as.data.frame()
rownames(htmp_RESP_ad) <- htmp_RESP_ad$SampleName
htmp_RESP_ad <- htmp_RESP_ad %>% dplyr::select(-SampleName)

htmp_RESP_ad_scaled <- scale(htmp_RESP_ad[,c(6:23)])
htmp_RESP_ad_metadata <- htmp_RESP_ad %>% dplyr::select(c(1:5))
htmp_RESP_ad_metadata <- htmp_RESP_ad_metadata %>% mutate_at(c(1:5), as.factor)

row_ha  = rowAnnotation(Treatment = htmp_RESP_ad_metadata$Treatment,col = list(Treatment = c("mock" = "#DEC08B", "betaglucan" = "#4D3C7E")))

col_fun_2 = colorRamp2(c(-4,-1.45, -0.12,0,0.12,1.45,4), c("#0765A0","#5d95b8","#FFF1E4", "#FFF1E4","#FFF1E4" ,"#de725f","#c61700"))
lgd = Legend(col_fun = col_fun_2, title = "z-score")


c <- Heatmap(htmp_RESP_ad_scaled, show_row_dend = TRUE, col = col_fun_2,
             clustering_method_columns = "complete", clustering_method_rows = "complete",
             column_dend_side = "top", column_dend_height = unit(0.5, "cm"),   
             column_title_gp = gpar(fontsize = 10), right_annotation = row_ha, column_names_gp = grid::gpar(fontsize = 5))

#Nucleoside related genes
genes_new_nuc <- subset(genes_new_adult_effect, (genes_new_adult_effect$SYMBOL %in% Nuc_ad_eff))
genes_new_nuc$Pathway <- c("Nucleotide metabolism related")

htmp_Nuc_ad <- genes_new_nuc %>% pivot_wider(values_from = "expression", names_from = "SYMBOL") %>% as.data.frame()
rownames(htmp_Nuc_ad) <- htmp_Nuc_ad$SampleName
htmp_Nuc_ad <- htmp_Nuc_ad %>% dplyr::select(-SampleName)

htmp_Nuc_ad_scaled <- scale(htmp_Nuc_ad[,c(6:23)])
htmp_Nuc_ad_metadata <- htmp_Nuc_ad %>% dplyr::select(c(1:5))
htmp_Nuc_ad_metadata <- htmp_Nuc_ad_metadata %>% mutate_at(c(1:5), as.factor)

row_ha  = rowAnnotation(Treatment = htmp_Nuc_ad_metadata$Treatment,col = list(Treatment = c("mock" = "#DEC08B", "betaglucan" = "#4D3C7E")))

col_fun_2 = colorRamp2(c(-4,-1.45, -0.12,0,0.12,1.45,4), c("#0765A0","#5d95b8","#FFF1E4", "#FFF1E4","#FFF1E4" ,"#de725f","#c61700"))
lgd = Legend(col_fun = col_fun_2, title = "z-score")


d <- Heatmap(htmp_Nuc_ad_scaled, show_row_dend = TRUE, col = col_fun_2,
             clustering_method_columns = "complete", clustering_method_rows = "complete",
             column_dend_side = "top", column_dend_height = unit(0.5, "cm"),   
             column_title_gp = gpar(fontsize = 10), right_annotation = row_ha, column_names_gp = grid::gpar(fontsize = 5))

ht_list = a + b + c +d

ht_list_TF_adults <- draw(ht_list)
```

![](Trained_immunity_GIT_files/figure-gfm/setup82-1.png)<!-- -->

``` r
graph2svg(ht_list_TF_adults, file = here::here("plots","ht_list_TF_adults"), width = 11.01, height = 1.725)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/ht_list_TF_adults.svg

### Plot how much a TF-factor influences a pathway of interest (adult effect)

``` r
Prot_ad_eff <- c(vector_ARNT_Proteasome,vector_SMAD5_ad_proteasome,vector_BACH2_ad_proteasome) %>% unique()
Ves_ad_eff <- c(vector_ARNT_Vesicle, vector_SMAD5_ad_Vesicle) %>% unique()
Resp_ad_eff <- c(vector_SMAD5_ad_cell_resp) %>% unique()
Nuc_ad_eff <- c(vector_ARNT_Nucleoside_metabolism, vector_SMAD5_Nucleoside_metabolism) %>% unique()


ARNT_PROT <- 100*(length(vector_ARNT_Proteasome)/length(Prot_ad_eff))
SMAD5_PROT <- 100*(length(vector_SMAD5_ad_proteasome)/length(Prot_ad_eff))
BACH2_PROT <- 100*(length(vector_BACH2_ad_proteasome)/length(Prot_ad_eff))

ARNT_VES <- 100*(length(vector_ARNT_Vesicle)/length(Ves_ad_eff))
SMAD5_VES <- 100*(length(vector_SMAD5_ad_Vesicle)/length(Ves_ad_eff))
BACH2_VES <- 100*(0/length(Ves_ad_eff))

ARNT_resp <- 100*(0/length(Resp_ad_eff))
SMAD5_resp <- 100*(length(vector_SMAD5_ad_cell_resp)/length(Resp_ad_eff))
BACH2_resp <- 100*(0/length(Resp_ad_eff))

ARNT_nucleo <- 100*(length(vector_ARNT_Nucleoside_metabolism)/length(Nuc_ad_eff))
SMAD5_nucleo <- 100*(length(vector_SMAD5_Nucleoside_metabolism)/length(Nuc_ad_eff))
BACH2_nucleo <- 100*(0/length(Nuc_ad_eff))



PROT <- as.data.frame(1:3)
PROT$Pathway <- c("PROT")
PROT$value <- c(ARNT_PROT,SMAD5_PROT, BACH2_PROT )
PROT$TF <- c("ARNT", "SMAD5", "BACH2")
PROT <- PROT %>% dplyr::select(-"1:3")
PROT <- PROT %>% mutate_at(c(2), as.numeric)

dotplot_influence_PROT <- ggplot(PROT, aes(x = value, y = TF, size = 1)) + 
  geom_point(stat = 'identity') + scale_y_discrete(limits = c("SMAD5", "BACH2", "ARNT")) +  
  xlab("% regulated genes") + ylab("Transcription factor") + xlim(25,105)+
  theme_bw()
print(dotplot_influence_PROT)
```

![](Trained_immunity_GIT_files/figure-gfm/setup83-1.png)<!-- -->

``` r
graph2svg(dotplot_influence_PROT, file = here::here("plots","dotplot_influence_PROT"), width =3.1, height = 0.87)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/dotplot_influence_PROT.svg

``` r
VES <- as.data.frame(1:3)
VES$Pathway <- c("VES")
VES$value <- c(ARNT_VES,SMAD5_VES, BACH2_VES)
VES$TF <- c("ARNT", "SMAD5", "BACH2")
VES <- VES %>% dplyr::select(-"1:3")
VES <- VES %>% mutate_at(c(2), as.numeric)

dotplot_influence_VES <- ggplot(VES, aes(x = value, y = TF, size = 1)) + 
  geom_point(stat = 'identity') + scale_y_discrete(limits = c("SMAD5", "BACH2", "ARNT")) +  
  xlab("% regulated genes") + ylab("Transcription factor") + xlim(25,105)+
  theme_bw()
print(dotplot_influence_VES)
```

![](Trained_immunity_GIT_files/figure-gfm/setup83-2.png)<!-- -->

``` r
graph2svg(dotplot_influence_VES, file = here::here("plots","dotplot_influence_VES"), width =2.7, height = 0.87)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/dotplot_influence_VES.svg

``` r
RESP <- as.data.frame(1:3)
RESP$Pathway <- c("RESP")
RESP$value <- c(ARNT_resp,SMAD5_resp, BACH2_resp)
RESP$TF <- c("ARNT", "SMAD5", "BACH2")
RESP <- RESP %>% dplyr::select(-"1:3")
RESP <- RESP %>% mutate_at(c(2), as.numeric)

dotplot_influence_RESP <- ggplot(RESP, aes(x = value, y = TF, size = 1)) + 
  geom_point(stat = 'identity') + scale_y_discrete(limits = c("SMAD5", "BACH2", "ARNT")) +  
  xlab("% regulated genes") + ylab("Transcription factor") + xlim(25,105)+
  theme_bw()
print(dotplot_influence_RESP)
```

![](Trained_immunity_GIT_files/figure-gfm/setup83-3.png)<!-- -->

``` r
graph2svg(dotplot_influence_RESP, file = here::here("plots","dotplot_influence_RESP"), width =2.0, height = 0.87)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/dotplot_influence_RESP.svg

``` r
NUCLEO <- as.data.frame(1:3)
NUCLEO$Pathway <- c("NUCLEO")
NUCLEO$value <- c(ARNT_nucleo,SMAD5_nucleo, BACH2_nucleo)
NUCLEO$TF <- c("ARNT", "SMAD5", "BACH2")
NUCLEO <- NUCLEO %>% dplyr::select(-"1:3")
NUCLEO <- NUCLEO %>% mutate_at(c(2), as.numeric)

dotplot_influence_NUCLEO <- ggplot(NUCLEO, aes(x = value, y = TF, size = 1)) + 
  geom_point(stat = 'identity') + scale_y_discrete(limits = c("SMAD5", "BACH2", "ARNT")) +  
  xlab("% regulated genes") + ylab("Transcription factor") + xlim(25,105)+
  theme_bw()
print(dotplot_influence_NUCLEO)
```

![](Trained_immunity_GIT_files/figure-gfm/setup83-4.png)<!-- -->

``` r
graph2svg(dotplot_influence_NUCLEO, file = here::here("plots","dotplot_influence_NUCLEO"), width =2.0, height = 0.87)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/dotplot_influence_NUCLEO.svg

### Analysis of enriched trancription factors in preterm vs. adults upon beta-glucan treatment

``` r
TF_link <- rio::import(here::here("data","TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv"))
TF_link <- TF_link %>% dplyr::filter(Organism == "Homo sapiens")
TF_Link <- TF_link %>% dplyr::select(Name.TF, Name.Target)

TF_link_preterm_HCFC1 <-  subset(TF_link, (TF_link$Name.TF %in% c("HCFC1")))
TF_link_preterm_ZNF804A <- subset(TF_link, (TF_link$Name.TF %in% c("ZNF804A"))) # 0 entries, excluded from subsequent analyses
TF_link_preterm_TFDP1 <- subset(TF_link, (TF_link$Name.TF %in% c("TFDP1"))) 
TF_link_preterm_RFXAP <-  subset(TF_link, (TF_link$Name.TF %in% c("RFXAP"))) # only 15 entries, excluded from subsequent analyses
TF_link_preterm_RFXANK <- subset(TF_link, (TF_link$Name.TF %in% c("RFXANK")))
TF_link_preterm_ARNT <- subset(TF_link, (TF_link$Name.TF %in% c("ARNT"))) 

#HCFC1
TF_link_preterm_HCFC1_GENEID <- TF_link_preterm_HCFC1 %>% dplyr::select(NCBI.GeneID.TF, NCBI.GeneID.Target)
Preterm_HCFC1_GENEID_TF <- TF_link_preterm_HCFC1_GENEID$NCBI.GeneID.TF %>% unique() %>% as.data.frame()
Preterm_HCFC1_GENEID_Target <- TF_link_preterm_HCFC1_GENEID$NCBI.GeneID.Target %>% unique() %>% as.data.frame()
colnames(Preterm_HCFC1_GENEID_Target) <- c("NCBI.GeneID.Target")
Preterm_HCFC1_GENEID_Target <- subset(Preterm_HCFC1_GENEID_Target, (Preterm_HCFC1_GENEID_Target$NCBI.GeneID.Target %in% vec_adult_preterm_beta_effect))
vec_Preterm_HCFC1_GENEID_Target <-  Preterm_HCFC1_GENEID_Target$NCBI.GeneID.Target  %>% unique()

#TFDP1
TF_link_preterm_TFDP1_GENEID <- TF_link_preterm_TFDP1 %>% dplyr::select(NCBI.GeneID.TF, NCBI.GeneID.Target)
Preterm_TFDP1_GENEID_TF <- TF_link_preterm_TFDP1_GENEID$NCBI.GeneID.TF %>% unique() %>% as.data.frame()
Preterm_TFDP1_GENEID_Target <- TF_link_preterm_TFDP1_GENEID$NCBI.GeneID.Target %>% unique() %>% as.data.frame()
colnames(Preterm_TFDP1_GENEID_Target) <- c("NCBI.GeneID.Target")
Preterm_TFDP1_GENEID_Target <- subset(Preterm_TFDP1_GENEID_Target, (Preterm_TFDP1_GENEID_Target$NCBI.GeneID.Target %in% vec_adult_preterm_beta_effect))
vec_Preterm_TFDP1_GENEID_Target <-  Preterm_TFDP1_GENEID_Target$NCBI.GeneID.Target  %>% unique()

#RFXANK
TF_link_preterm_RFXANK_GENEID <- TF_link_preterm_RFXANK %>% dplyr::select(NCBI.GeneID.TF, NCBI.GeneID.Target)
Preterm_RFXANK_GENEID_TF <- TF_link_preterm_RFXANK_GENEID$NCBI.GeneID.TF %>% unique() %>% as.data.frame()
Preterm_RFXANK_GENEID_Target <- TF_link_preterm_RFXANK_GENEID$NCBI.GeneID.Target %>% unique() %>% as.data.frame()
colnames(Preterm_RFXANK_GENEID_Target) <- c("NCBI.GeneID.Target")
Preterm_RFXANK_GENEID_Target <- subset(Preterm_RFXANK_GENEID_Target, (Preterm_RFXANK_GENEID_Target$NCBI.GeneID.Target %in% vec_adult_preterm_beta_effect))
vec_Preterm_RFXANK_GENEID_Target <-  Preterm_RFXANK_GENEID_Target$NCBI.GeneID.Target %>% unique()

#ARNT
TF_link_preterm_ARNT_GENEID <- TF_link_preterm_ARNT %>% dplyr::select(NCBI.GeneID.TF, NCBI.GeneID.Target)
Preterm_ARNT_GENEID_TF <- TF_link_preterm_ARNT_GENEID$NCBI.GeneID.TF %>% unique() %>% as.data.frame()
Preterm_ARNT_GENEID_Target <- TF_link_preterm_ARNT_GENEID$NCBI.GeneID.Target %>% unique() %>% as.data.frame()
colnames(Preterm_ARNT_GENEID_Target) <- c("NCBI.GeneID.Target")
Preterm_ARNT_GENEID_Target <- subset(Preterm_ARNT_GENEID_Target, (Preterm_ARNT_GENEID_Target$NCBI.GeneID.Target %in% vec_adult_preterm_beta_effect))
vec_Preterm_ARNT_GENEID_Target <-  Preterm_ARNT_GENEID_Target$NCBI.GeneID.Target %>% unique()

#ORA

#set background
background_genes <- readRDS(here::here("data","se.rds"))
background_genes <- background_genes@rowRanges@partitioning@NAMES %>% as.data.frame()

background_genes$ENSEMBL <- str_replace(background_genes$.,
                                        pattern = ".[0-9]+$",
                                        replacement = "")

background_genes <- background_genes %>% dplyr::select(-.)

"entrez_id" = mapIds(
  # Replace with annotation package for the organism relevant to your data
  org.Hs.eg.db,
  keys =  background_genes$ENSEMBL,
  # Replace with the type of gene identifiers in your data
  keytype = "ENSEMBL",
  # Replace with the type of gene identifiers you would like to map to
  column = "ENTREZID",
  # This will keep only the first mapped value for each Ensembl ID
  multiVals = "first"
)
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
background_genes$ENTREZ <- entrez_id
background_genes <- background_genes %>% dplyr::filter(!is.na(ENTREZ)) %>% dplyr::select(-ENSEMBL)
background_genes <- background_genes$ENTREZ %>% unique()


list_DF_preterm_adult_beta <- list(vec_Preterm_HCFC1_GENEID_Target = vec_Preterm_HCFC1_GENEID_Target, vec_Preterm_TFDP1_GENEID_Target= vec_Preterm_TFDP1_GENEID_Target, 
                vec_Preterm_RFXANK_GENEID_Target = vec_Preterm_RFXANK_GENEID_Target, vec_Preterm_ARNT_GENEID_Target=vec_Preterm_ARNT_GENEID_Target )

preterm_adult_beta_BP <- list()
for (i in 1:4) {
  enrich_preterm_adult_beta_BP <- enrichGO(gene = list_DF_preterm_adult_beta[[i]],
                                           universe = background_genes,
                                           OrgDb = org.Hs.eg.db,
                                           ont = "BP",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff = 0.1,
                                           qvalueCutoff = 0.1,
                                           readable = TRUE)
  preterm_adult_beta_BP[[i]] <- enrich_preterm_adult_beta_BP
}
#saveRDS(preterm_adult_beta_BP, "C:/Users/micha/OneDrive/Dokumente/R-Trained_Immunity/data/RDS_objects/TF_stuff/preterm_adult_beta_BP.rds")
#preterm_adult_beta_BP <- rio::import("C:/Users/micha/OneDrive/Dokumente/R-Trained_Immunity/data/RDS_objects/TF_stuff/preterm_adult_beta_BP.rds")


#HCFC1
HCFC1_BP_preterm <- preterm_adult_beta_BP[[1]] %>% dplyr::filter(preterm_adult_beta_BP[[1]]@result[["p.adjust"]] <= 0.1)
dotplot_HCFC1 <- enrichplot::dotplot(HCFC1_BP_preterm)

HCFC1_adPre <- dotplot_HCFC1$data %>% as.data.frame()
HCFC1_adPre_cellcycle <- HCFC1_adPre[c(1,2,3,5:10),]
HCFC1_adPre_generalmetab <- HCFC1_adPre[c(4),]

vector_HCFC1_adPre_cellcycle <- HCFC1_adPre_cellcycle$geneID
vector_HCFC1_adPre_cellcycle <- strsplit(vector_HCFC1_adPre_cellcycle, "/") %>% unlist() %>% unique()

vector_HCFC1_adPre_generalmetab <- HCFC1_adPre_generalmetab$geneID
vector_HCFC1_adPre_generalmetab <- strsplit(vector_HCFC1_adPre_generalmetab, "/") %>% unlist() %>% unique()

#TFDP1
TFDP1_BP_preterm <- preterm_adult_beta_BP[[2]] %>% dplyr::filter(preterm_adult_beta_BP[[2]]@result[["p.adjust"]] <= 0.1)
dotplot_TFDP1 <- enrichplot::dotplot(TFDP1_BP_preterm)

TFDP1_adPre <- dotplot_TFDP1$data %>% as.data.frame()
TFDP1_adPre_cellcycle <- TFDP1_adPre[c(1:10),]

vector_TFDP1_adPre_cellcycle <- TFDP1_adPre_cellcycle$geneID
vector_TFDP1_adPre_cellcycle <- strsplit(vector_TFDP1_adPre_cellcycle, "/") %>% unlist() %>% unique()

#RFXANK
RFXANK_BP_preterm <- preterm_adult_beta_BP[[3]] %>% dplyr::filter(preterm_adult_beta_BP[[3]]@result[["p.adjust"]] <= 0.1)
dotplot_RFXANK <- enrichplot::dotplot(RFXANK_BP_preterm)

RFXANK_adPre <- dotplot_RFXANK$data %>% as.data.frame()
RFXANK_adPre_generalmetab <- RFXANK_adPre[c(2,3,4),]

vector_RFXANK_adPre_generalmetab <- RFXANK_adPre_generalmetab$geneID
vector_RFXANK_adPre_generalmetab <- strsplit(vector_RFXANK_adPre_generalmetab, "/") %>% unlist() %>% unique()

#ARNT
ARNT_BP_preterm <- preterm_adult_beta_BP[[4]] %>% dplyr::filter(preterm_adult_beta_BP[[4]]@result[["p.adjust"]] <= 0.1)
dotplot_ARNT <- enrichplot::dotplot(ARNT_BP_preterm)

ARNT_adPre <- dotplot_ARNT$data %>% as.data.frame()
ARNT_adPre_Cell_Cycle <- ARNT_adPre[c(1,4,5,9,10),]
ARNT_adPre_generalmetab <- ARNT_adPre[c(2,3,6,8),]
ARNT_adPre_Inflammation <- ARNT_adPre[c(7),]

vector_ARNT_adPre_Cell_Cycle <- ARNT_adPre_Cell_Cycle$geneID
vector_ARNT_adPre_Cell_Cycle <- strsplit(vector_ARNT_adPre_Cell_Cycle, "/") %>% unlist() %>% unique()

vector_ARNT_adPre_generalmetab <- ARNT_adPre_generalmetab$geneID
vector_ARNT_adPre_generalmetab <- strsplit(vector_ARNT_adPre_generalmetab, "/") %>% unlist() %>% unique()

vector_ARNT_adPre_Inflammation <- ARNT_adPre_Inflammation$geneID
vector_ARNT_adPre_Inflammation <- strsplit(vector_ARNT_adPre_Inflammation, "/") %>% unlist() %>% unique()


#define overall vectors
Cell_Cycle_Ad_Prem_eff <- c(vector_HCFC1_adPre_cellcycle, vector_TFDP1_adPre_cellcycle, vector_ARNT_adPre_Cell_Cycle) %>% unique()
Gen_metab_Ad_Prem_eff <- c(vector_HCFC1_adPre_generalmetab, vector_RFXANK_adPre_generalmetab, vector_ARNT_adPre_generalmetab) %>% unique()
Inflammation_Ad_Prem_eff <- c(vector_ARNT_adPre_Inflammation) %>% unique()
```

### Heatmap for adult vs. preterm (beta glucan-treated)

``` r
x_lcpm <- cpm(x)

genes <- x_lcpm %>% as.data.frame()
genes$ENSEMBL <- rownames(genes)
genes$ENSEMBL <- str_replace(genes$ENSEMBL, pattern = ".[0-9]+$", replacement = "")
rownames(genes) <- genes$ENSEMBL

genes_new <- genes
genes_new$SYMBOL <- mapIds(org.Hs.eg.db, keys = genes_new$ENSEMBL, keytype="ENSEMBL", column = "SYMBOL")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
genes_new <- genes_new %>% drop_na(SYMBOL)

#random sampling of duplicate entries
genes_new <- genes_new %>%
  group_by(SYMBOL) %>%
  sample_n(1) %>% as.data.frame()

#get metadata into dataframe
rownames(genes_new) <- genes_new$SYMBOL
genes_new <- t(genes_new) %>% as.data.frame()
genes_new <- genes_new %>% rownames_to_column("SampleName")
genes_new <- genes_new[-c(27:28),]
genes_new <- genes_new %>% left_join(metadata, by = "SampleName")
rownames(genes_new) <- genes_new$SampleName 
genes_new <- genes_new %>% mutate_at(c(2:14513), as.numeric)
genes_new <- genes_new %>% mutate_at(c(14514:14517), as.factor)

genes_new_adPre_effect <- genes_new %>% dplyr::filter(genes_new$Group == "adult" | genes_new$Group == "preterm")
genes_new_adPre_effect <- genes_new_adPre_effect %>% dplyr::filter(genes_new_adPre_effect$Treatment == "betaglucan")
genes_new_adPre_effect <- genes_new_adPre_effect %>% pivot_longer(cols = 2:14513, values_to = "expression", names_to = "SYMBOL")


#cellcycle
genes_new_cellcycle <- subset(genes_new_adPre_effect, (genes_new_adPre_effect$SYMBOL %in% Cell_Cycle_Ad_Prem_eff))
genes_new_cellcycle$Pathway <- c("Cell cycle")
htmp_cellcycle_AdPre <- genes_new_cellcycle %>% pivot_wider(values_from = "expression", names_from = "SYMBOL") %>% as.data.frame()
rownames(htmp_cellcycle_AdPre) <- htmp_cellcycle_AdPre$SampleName
htmp_cellcycle_AdPre <- htmp_cellcycle_AdPre %>% dplyr::select(-SampleName)

htmp_cellcycle_AdPre_scaled <- scale(htmp_cellcycle_AdPre[,c(6:77)])
htmp_cellcycle_AdPre_metadata <- htmp_cellcycle_AdPre %>% dplyr::select(c(1:5))
htmp_cellcycle_AdPre_metadata <- htmp_cellcycle_AdPre_metadata %>% mutate_at(c(1:5), as.factor)

row_ha  = rowAnnotation(Treatment = htmp_cellcycle_AdPre_metadata$Group,col = list(Treatment = c("adult" = "#b2c5b3", "preterm" = "#ce8793")))

col_fun_2 = colorRamp2(c(-4,-1.45, -0.12,0,0.12,1.45,4), c("#0765A0","#5d95b8","#FFF1E4", "#FFF1E4","#FFF1E4" ,"#de725f","#c61700"))
lgd = Legend(col_fun = col_fun_2, title = "z-score")


a <- Heatmap(htmp_cellcycle_AdPre_scaled, show_row_dend = TRUE, col = col_fun_2,
             clustering_method_columns = "complete", clustering_method_rows = "complete",
             column_dend_side = "top", column_dend_height = unit(0.5, "cm"),   
             column_title_gp = gpar(fontsize = 10), right_annotation = row_ha, column_names_gp = grid::gpar(fontsize = 5))


#general metabolism
genes_new_Gen_metab <- subset(genes_new_adPre_effect, (genes_new_adPre_effect$SYMBOL %in% Gen_metab_Ad_Prem_eff))
genes_new_Gen_metab$Pathway <- c("general metabolism")
htmp_Gen_metab_AdPre <- genes_new_Gen_metab %>% pivot_wider(values_from = "expression", names_from = "SYMBOL") %>% as.data.frame()
rownames(htmp_Gen_metab_AdPre) <- htmp_Gen_metab_AdPre$SampleName
htmp_Gen_metab_AdPre <- htmp_Gen_metab_AdPre %>% dplyr::select(-SampleName)

htmp_Gen_metab_AdPre_scaled <- scale(htmp_Gen_metab_AdPre[,c(6:46)])
htmp_Gen_metab_AdPre_metadata <- htmp_Gen_metab_AdPre %>% dplyr::select(c(1:5))
htmp_Gen_metab_AdPre_metadata <- htmp_Gen_metab_AdPre_metadata %>% mutate_at(c(1:5), as.factor)

row_ha  = rowAnnotation(Treatment = htmp_Gen_metab_AdPre_metadata$Group,col = list(Treatment = c("adult" = "#b2c5b3", "preterm" = "#ce8793")))

col_fun_2 = colorRamp2(c(-4,-1.45, -0.12,0,0.12,1.45,4), c("#0765A0","#5d95b8","#FFF1E4", "#FFF1E4","#FFF1E4" ,"#de725f","#c61700"))
lgd = Legend(col_fun = col_fun_2, title = "z-score")


b <- Heatmap(htmp_Gen_metab_AdPre_scaled, show_row_dend = TRUE, col = col_fun_2,
             clustering_method_columns = "complete", clustering_method_rows = "complete",
             column_dend_side = "top", column_dend_height = unit(0.5, "cm"),   
             column_title_gp = gpar(fontsize = 10), right_annotation = row_ha, column_names_gp = grid::gpar(fontsize = 5))


#Inflammation
genes_new_Inflammation <- subset(genes_new_adPre_effect, (genes_new_adPre_effect$SYMBOL %in% Inflammation_Ad_Prem_eff))
genes_new_Inflammation$Pathway <- c("Inflammation")
htmp_Gen_Inflammation_AdPre <- genes_new_Inflammation %>% pivot_wider(values_from = "expression", names_from = "SYMBOL") %>% as.data.frame()
rownames(htmp_Gen_Inflammation_AdPre) <- htmp_Gen_Inflammation_AdPre$SampleName
htmp_Gen_Inflammation_AdPre <- htmp_Gen_Inflammation_AdPre %>% dplyr::select(-SampleName)

htmp_Gen_Inflammation_AdPre_scaled <- scale(htmp_Gen_Inflammation_AdPre[,c(6:40)])
htmp_Gen_Inflammation_AdPre_metadata <- htmp_Gen_Inflammation_AdPre %>% dplyr::select(c(1:5))
htmp_Gen_Inflammation_AdPre_metadata <- htmp_Gen_Inflammation_AdPre_metadata %>% mutate_at(c(1:5), as.factor)

row_ha  = rowAnnotation(Treatment = htmp_Gen_Inflammation_AdPre_metadata$Group,col = list(Treatment = c("adult" = "#b2c5b3", "preterm" = "#ce8793")))

col_fun_2 = colorRamp2(c(-4,-1.45, -0.12,0,0.12,1.45,4), c("#0765A0","#5d95b8","#FFF1E4", "#FFF1E4","#FFF1E4" ,"#de725f","#c61700"))
lgd = Legend(col_fun = col_fun_2, title = "z-score")


c <- Heatmap(htmp_Gen_Inflammation_AdPre_scaled, show_row_dend = TRUE, col = col_fun_2,
             clustering_method_columns = "complete", clustering_method_rows = "complete",
             column_dend_side = "top", column_dend_height = unit(0.5, "cm"),   
             column_title_gp = gpar(fontsize = 10), right_annotation = row_ha, column_names_gp = grid::gpar(fontsize = 5))



ht_list_ad_pre = a + b + c

ht_list_TF_ad_pre <- draw(ht_list_ad_pre)
```

![](Trained_immunity_GIT_files/figure-gfm/setup85-1.png)<!-- -->

``` r
graph2svg(ht_list_TF_ad_pre, file =  here::here("plots","ht_list_TF_ad_pre"), width = 10.170, height = 1.585)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/ht_list_TF_ad_pre.svg

### TF-factor influence on a pathway of interest (adult vs. preterm, beta glucan treated)

``` r
Cell_Cycle_Ad_Prem_eff <- c(vector_HCFC1_adPre_cellcycle, vector_TFDP1_adPre_cellcycle, vector_ARNT_adPre_Cell_Cycle) %>% unique()
Gen_metab_Ad_Prem_eff <- c(vector_HCFC1_adPre_generalmetab, vector_RFXANK_adPre_generalmetab, vector_ARNT_adPre_generalmetab) %>% unique()
Inflammation_Ad_Prem_eff <- c(vector_ARNT_adPre_Inflammation) %>% unique()


HCFC1_Cycle <- 100*(length(vector_HCFC1_adPre_cellcycle)/length(Cell_Cycle_Ad_Prem_eff))
TFDP1_Cycle <- 100*(length(vector_TFDP1_adPre_cellcycle)/length(Cell_Cycle_Ad_Prem_eff))
RFXANK_Cycle <- 100*(0/length(Cell_Cycle_Ad_Prem_eff))
ARNT_Cycle <- 100*(length(vector_ARNT_adPre_Cell_Cycle)/length(Cell_Cycle_Ad_Prem_eff))

HCFC1_met <- 100*(length(vector_HCFC1_adPre_generalmetab)/length(Gen_metab_Ad_Prem_eff))
TFDP1_met <- 100*(0/length(Gen_metab_Ad_Prem_eff))
RFXANK_met <- 100*(length(vector_RFXANK_adPre_generalmetab)/length(Gen_metab_Ad_Prem_eff))
ARNT_met <- 100*(length(vector_ARNT_adPre_generalmetab)/length(Gen_metab_Ad_Prem_eff))

HCFC1_Inflamm <- 100*(0/length(Inflammation_Ad_Prem_eff))
TFDP1_Inflamm <- 100*(0/length(Inflammation_Ad_Prem_eff))
RFXANK_Inflamm <- 100*(0/length(Inflammation_Ad_Prem_eff))
ARNT_Inflamm <- 100*(length(vector_ARNT_adPre_Inflammation)/length(Inflammation_Ad_Prem_eff))

CYCLE <- as.data.frame(1:4)
CYCLE$Pathway <- c("CYCLE")
CYCLE$value <- c(HCFC1_Cycle,TFDP1_Cycle, RFXANK_Cycle, ARNT_Cycle)
CYCLE$TF <- c("HCFC1", "TFDP1", "RFXANK", "ARNT")
CYCLE <- CYCLE %>% dplyr::select(-"1:4")
CYCLE <- CYCLE %>% mutate_at(c(2), as.numeric)

dotplot_influence_CYCLE <- ggplot(CYCLE, aes(x = value, y = TF, size = 1)) + 
  geom_point(stat = 'identity') + scale_y_discrete(limits = c("RFXANK", "ARNT", "TFDP1", "HCFC1")) + 
  xlab("% regulated genes") + ylab("Transcription factor") + xlim(25,105)+
  theme_bw()
print(dotplot_influence_CYCLE)
```

![](Trained_immunity_GIT_files/figure-gfm/setup86-1.png)<!-- -->

``` r
graph2svg(dotplot_influence_CYCLE, file = here::here("plots","dotplot_influence_CYCLE"), width = 3.6, height = 0.85)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/dotplot_influence_CYCLE.svg

``` r
MET <- as.data.frame(1:4)
MET$Pathway <- c("MET")
MET$value <- c(HCFC1_met,TFDP1_met, RFXANK_met, ARNT_met)
MET$TF <- c("HCFC1", "TFDP1", "RFXANK", "ARNT")
MET <- MET %>% dplyr::select(-"1:4")
MET <- MET %>% mutate_at(c(2), as.numeric)

dotplot_influence_MET <- ggplot(MET, aes(x = value, y = TF, size = 1)) + 
  geom_point(stat = 'identity') + scale_y_discrete(limits = c("RFXANK", "ARNT", "TFDP1", "HCFC1")) + 
  xlab("% regulated genes") + ylab("Transcription factor") + xlim(25,105)+
  theme_bw()
print(dotplot_influence_MET)
```

![](Trained_immunity_GIT_files/figure-gfm/setup86-2.png)<!-- -->

``` r
graph2svg(dotplot_influence_MET, file =  here::here("plots","dotplot_influence_MET"), width = 2.90, height = 0.85)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/dotplot_influence_MET.svg

``` r
Inflamm <- as.data.frame(1:4)
Inflamm$Pathway <- c("Inflammation")
Inflamm$value <- c(HCFC1_Inflamm,TFDP1_Inflamm, RFXANK_Inflamm, ARNT_Inflamm)
Inflamm$TF <- c("HCFC1", "TFDP1", "RFXANK", "ARNT")
Inflamm <- Inflamm %>% dplyr::select(-"1:4")
Inflamm <- Inflamm %>% mutate_at(c(2), as.numeric)

dotplot_influence_Inflamm <- ggplot(Inflamm, aes(x = value, y = TF, size = 1)) + 
  geom_point(stat = 'identity') + scale_y_discrete(limits = c("RFXANK", "ARNT", "TFDP1", "HCFC1")) +
  xlab("% regulated genes") + ylab("Transcription factor") + xlim(25,105)+
  theme_bw()
print(dotplot_influence_Inflamm)
```

![](Trained_immunity_GIT_files/figure-gfm/setup86-3.png)<!-- -->

``` r
graph2svg(dotplot_influence_Inflamm, file = here::here("plots","dotplot_influence_Inflamm"), width = 2.70, height = 0.85)
```

    ## Exported graph as C:/Users/Michi/Documents/Trained_Immunity_GIT/Trained_immunity_GIT/plots/dotplot_influence_Inflamm.svg
