library(data.table)
library(jsonlite)
library(dplyr)
library(ggplot2)

biospecimen.file <- "biospecimen.cart.2022-05-18.tar.gz"
clinical.file <- "clinical.cart.2022-05-18.tar.gz"
manifest.file <- "gdc_manifest_20220518_055903.txt"
metadata.file <- "metadata.cart.2022-05-18.json"
sample_sheet.file <- "gdc_sample_sheet.2022-05-18.tsv"

####################################################
# Merge manifest and sample sheet
####################################################
# read manifest
manifest <- fread(manifest.file)

# read sample_sheet
sample_sheet <- fread(sample_sheet.file)
names(sample_sheet) <- gsub(" ", "_", tolower(names(sample_sheet)))
names(sample_sheet)[which(names(sample_sheet) == "case_id")] <- "case_submitter_id"
names(sample_sheet)[which(names(sample_sheet) == "sample_id")] <- "sample_submitter_id"

# merge manifest and sample_sheet
manifest <- merge(manifest, sample_sheet[, -"file_name"], by.x = "id", by.y = "file_id")

# it is possible that multiple case/sample entries link to the same file. This is normally 
# due to tumor-normal paired analysis, or project with RPPA data when multiple samples were
# merged into one aliquot. We should clean up cases

# validation: check if all case entries for the same file is unique. "0" is expected.
print(sum(sapply(gsub(" ", "", manifest$case_submitter_id), function(x) length(unique(strsplit(x, ",")[[1]])) > 1)))
manifest$case_submitter_id <- sapply(manifest$case, function(x) strsplit(x, ",")[[1]][1])

# categorize data
manifest <- manifest %>%
              mutate(type = case_when(
                              grepl("wxs.aliquot_ensemble_raw.maf.gz", filename) ~ "wxs_protected_maf",
                              grepl("wxs.aliquot_ensemble_masked.maf.gz", filename) ~ "wxs_public_maf",
                              grepl("wgs.ASCAT.gene_level.copy_number_variation.tsv", filename) ~ "ascat_gene_level_cnv",
                              grepl("arriba.rna_fusion.bedpe", filename) ~ "arriba_bedpe",
                              grepl("star_fusion.rna_fusion.bedpe", filename) ~ "star_fusion_bedpe",
                              grepl("wgs.sanger_raw_pindel.raw_somatic_mutation.vcf.gz", filename) ~ "wgs_pindel_vcf",
                              grepl("wgs.CaVEMan.raw_somatic_mutation.vcf.gz", filename) ~ "wgs_caveman_vcf",
                              grepl("rna_seq.star_splice_junctions.tsv.gz", filename) ~ "splice_junction",
                              grepl("rna_seq.augmented_star_gene_counts.tsv", filename) ~ "star_count",
                              grepl("qc_filtered_feature_bc_matrix.tar.gz", filename) ~ "scrna_filtered_count",
                              grepl("seurat.deg.tsv", filename) ~ "scrna_deg",
                              grepl("seurat.analysis.tsv", filename) ~ "scrna_cluster",
                              grepl("seurat.loom", filename) ~ "seurat_loom",
                              TRUE ~ "Unknown"
                              )) %>%
              setDT()


# inspection of sample_type
manifest$sample_type_unique <- sapply(manifest$sample_type, function(x) paste(sort(unique(strsplit(gsub(", ", ",", x), ",")[[1]])), collapse=","))
print(table(manifest$sample_type_unique))
# remove "Solid Tissue Normal"
manifest <- manifest[sample_type_unique != "Solid Tissue Normal"]

# add aliquot to the manifest
metadata <- fromJSON(metadata.file)
d <- data.table(id = metadata$file_id, aliquot = sapply(metadata$associated_entities, function(x) paste(x$entity_submitter_id, collapse=","))
)
manifest <- merge(manifest, d, by="id")


####################################################
# Generate manifest and mapping files for Alec project
####################################################

# subset of wgs/rna-seq manifest for Alec
alec.rna <- manifest[type %in% c("star_count")]
alec.wgs <- manifest[type %in% c("ascat_gene_level_cnv", "wgs_caveman_vcf", "wgs_pindel_vcf")]

# remove those aliquot from multiple samples
alec.rna <- alec.rna[nchar(sample_submitter_id) == 12]
alec.wgs <- alec.wgs[nchar(sample_submitter_id) == 26]

# filter for common aliquots
w <- which(alec.wgs$sample_type == "Blood Derived Normal, Primary Tumor")
alec.wgs$tumor_sample <- sapply(alec.wgs$sample_submitter_id, function(x) strsplit(gsub(" ", "", x), ",")[[1]][1])
alec.wgs$tumor_sample[w] <- sapply(alec.wgs$sample_submitter_id[w], function(x) strsplit(gsub(" ", "", x), ",")[[1]][2])
common.sample <- unique(intersect(alec.wgs$tumor_sample, alec.rna$sample_submitter_id))
alec.wgs <- alec.wgs[tumor_sample %in% common.sample]
temp <- table(alec.wgs$tumor_sample)
alec.wgs <- alec.wgs[(!tumor_sample %in% names(temp[temp<3])) & 
                     (! aliquot %in% c("CPT0008800006_1,CPT0008850002", "CPT0065640009_1,CPT0070720002", "CPT0078000006,CPT0078080002", "CPT0088320009_1,CPT0088350002", "CPT0124360006_1,CPT0124390002", "CPT0125090006_1,CPT0125120002", "CPT0125130006_1,CPT0125150002"))]
common.sample <- unique(alec.wgs$tumor_sample)
alec.rna <- alec.rna[sample_submitter_id %in% common.sample]

alec.mapping <- alec.rna[match(common.sample, sample_submitter_id), c("case_submitter_id", "sample_submitter_id", "id")]
names(alec.mapping)[3] <- "rna_expression_file_id"
alec.mapping$ascat_cnv_file_id <- alec.wgs[type=="ascat_gene_level_cnv"][match(common.sample, tumor_sample)]$id
alec.mapping$caveman_vcf_file_id <- alec.wgs[type=="wgs_caveman_vcf"][match(common.sample, tumor_sample)]$id
alec.mapping$pindel_vcf_file_id <- alec.wgs[type=="wgs_pindel_vcf"][match(common.sample, tumor_sample)]$id

alec.manifest <- rbind(alec.rna[, 1:5], alec.wgs[, 1:5])

write.table(alec.manifest, "alec.wgs_rna.manifest.dr33.tsv", col.names=T, row.names=F, sep="\t", quote=F)
write.table(alec.mapping, "alec.wgs_rna.mapping.dr33.tsv", col.names=T, row.names=F, sep="\t", quote=F)


####################################################
# Generate manifest and mapping files v22 vs v36 project
####################################################

# read cptac-3 DR32 vs DR31 mapping file, and filter for those DR32 files with a v22 version
mapping <- fread("https://github.com/NCI-GDC/gdc-docs/blob/develop/docs/Data/Release_Notes/GCv36_Manifests/CPTAC-3.tsv?raw=true")
# these removed WGS VCFs and scRNA-Seq. There is DR33 release of scRNA-Seq.
manifest <- manifest[id %in% mapping$new_id]

# remove duplicate data on the same case
manifest$key <- paste(manifest$case_submitter_id, manifest$type)
manifest <- manifest[order(case_submitter_id, type, aliquot)]
manifest$duplicate <- duplicated(manifest$key)
manifest <- manifest[duplicate == F]
# filter for complete case (case with all 7 data types)
temp <- table(manifest$case_submitter_id)
manifest <- manifest[case_submitter_id %in% names(temp[temp==7])]
manifest <- manifest[, -c("duplicate", "key", "sample_type_unique")]

mapping <- mapping[, c("new_id", "id", "filename", "md5", "size")]
names(mapping) <- c("id", "id_v22", "filename_v22", "md5_v22", "size_v22")
manifest <- merge(manifest, mapping)
write.table(manifest, "cptac-3_augmented_manifest.dr32.tsv", col.names=T, row.names=F, sep="\t", quote=F)




