#!/usr/bin/env RScript --vanilla

suppressPackageStartupMessages({
  source('global.R')
})

# Inputs ------------------------------------------------------------------

# example:
# ./run-batch.R ~/Downloads/multieditr_test/hGli1_3h_control_F.ab1 ~/Downloads/multieditr_test/hGli1_3h_0.25nM_F.ab1,~/Downloads/multieditr_test/hGli1_3h_0.5nM_F.ab1,~/Downloads/multieditr_test/hGli1_3hr_0.1nM_F.ab1 "TGCTAGAGG" "A|G"

args = commandArgs(trailingOnly = TRUE)
ctrl = args[1] # control file
samples = args[2] # comma-separated input files
motif = args[3] # motif of interest
wt = strsplit(args[4], '\\|')[[1]][1] # wt|edit
edit = strsplit(args[4], '\\|')[[1]][2] # wt|edit
phred = 1e-4
p_value = 1e-4
adjust_p = TRUE

# Process -----------------------------------------------------------------

# there is only one control // reading from .ab1 file
ctrl_sanger = readsangerseq(ctrl)
ctrl_df = make_ctrl_sanger_df(ctrl_sanger)
init_ctrl_seq = ctrl_df$base_call %>% paste0(., collapse = "")
ctrl_fastq = abif_to_fastq(path = ctrl, cutoff = phred)
ctrl_alignment = pairwiseAlignment(ctrl_fastq$seq, subject = init_ctrl_seq)

# multiple inputs, just loop through
res = map(strsplit(samples, ',')[[1]], function(sample) {

  ctrl_df = make_ctrl_sanger_df(ctrl_sanger) # need to reinitiate this, gets overwritten

  ### read
  sample_sanger = readsangerseq(sample)
  sample_df = suppressWarnings(make_samp_sanger_df(sample_sanger, init_ctrl_seq))
  init_sample_seq = paste0(sample_df$primary_base_call, collapse = '')
  sample_fastq = abif_to_fastq(path = sample, cutoff = phred)
  sample_alignment = pairwiseAlignment(sample_fastq$seq, subject = init_sample_seq)

  ### filter
  sample_df_orig = sample_df
  filteredData = list(
    'sample_df' = (sample_df %>%
                     filter(index >= sample_alignment@subject@range@start) %>%
                     filter(index <= sample_alignment@subject@range@start +
                              sample_alignment@subject@range@width - 1) %>%
                     mutate(post_filter_index = 1:NROW(index))),
    'ctrl_df' = (ctrl_df %>%
                   filter(index >= ctrl_alignment@subject@range@start) %>%
                   filter(index <= ctrl_alignment@subject@range@start +
                            ctrl_alignment@subject@range@width - 1) %>%
                   mutate(post_filter_index = 1:NROW(index)))
  )

  filteredData$ctrl_seq = filteredData$ctrl_df$base_call %>% paste0(., collapse = "")
  filteredData$sample_seq = filteredData$sample_df$max_base %>% paste0(., collapse = "")
  filteredData$pre_cross_align_sample_df = filteredData$sample_df

  ### trim
  trimmed_alignment = align_and_trim(filteredData$sample_seq,
                                     filteredData$ctrl_seq,
                                     min_continuity = 15)

  samp_alignment_seq = trimmed_alignment$alignment@pattern %>% as.character()
  ctrl_alignment_seq = trimmed_alignment$alignment@subject %>% as.character()

  sample_trimmed_alignment = pairwiseAlignment(pattern = trimmed_alignment$pattern,
                                               subject = filteredData$sample_seq,
                                               gapOpening = 1000,
                                               gapExtension = 1000,
                                               type = "local")

  ctrl_trimmed_alignment = pairwiseAlignment(pattern = trimmed_alignment$subject,
                                             subject = filteredData$ctrl_seq,
                                             gapOpening = 1000,
                                             gapExtension = 1000,
                                             type = "local")

  sample_df = filteredData$sample_df %>%
    dplyr::select(-(A_area:max_base_height), -post_filter_index, index) %>%
    bind_rows(., ., ., .) %>%
    mutate(base = rep(ACGT, each = NROW(index)/4), max_base = base) %>%
    mutate_if(is.character, nucleotide_factor) %>%
    mutate(eta = 1) %>%
    dplyr::select(-base) %>%
    spread(max_base, eta) %>%
    inner_join(filteredData$sample_df, .) %>%
    dplyr::rename(A_eta = A, C_eta = C, G_eta = G, T_eta = `T`) %>%
    mutate(pred_height = {ifelse(max_base == "A", A_eta,
                                 ifelse(max_base == "C", C_eta,
                                        ifelse(max_base == "G", G_eta,
                                               ifelse(max_base == "T", `T_eta`,NA))))}) %>%
    dplyr::select(A_area:T_perc,
                  max_base, Tot.Area, index, max_base_height, post_filter_index, A_eta:T_eta, pred_height)

  sample_df %<>%
    filter(post_filter_index >= sample_trimmed_alignment@subject@range@start) %<>%
    filter(post_filter_index <= sample_trimmed_alignment@subject@range@start + sample_trimmed_alignment@subject@range@width - 1) %>%
    mutate(post_aligned_index = 1:NROW(index))

  ctrl_df = filteredData$ctrl_df %>%
    filter(post_filter_index >= ctrl_trimmed_alignment@subject@range@start) %>%
    filter(post_filter_index <= ctrl_trimmed_alignment@subject@range@start + ctrl_trimmed_alignment@subject@range@width - 1) %>%
    mutate(post_aligned_index = 1:NROW(index))

  ctrl_seq = ctrl_df$base_call %>% paste0(., collapse = "")
  samp_seq = sample_df$max_base %>% paste0(., collapse = "")

  pre_aligned_sample_df = sample_df

  samp_indel = samp_alignment_seq %>% gregexpr("-", .) %>% unlist
  ctrl_indel = ctrl_alignment_seq %>% gregexpr("-", .) %>% unlist

  if (samp_indel == -1) { samp_indel = 0 } else {}
  if (ctrl_indel == -1) { ctrl_indel = 0 } else {}

  ctrl_df = ctrl_df %>%
    mutate(.,
           indel_filter = ctrl_alignment_seq %>%
             subchar(., samp_indel, "_") %>%
             gsub("-", "", .) %>%
             strsplit(x = ., split = "") %>%
             .[[1]]
    ) %>%
    filter(indel_filter != "_") %>%
    dplyr::select(-indel_filter) %>%
    mutate(filtered_index = 1:NROW(max_base))

  sample_df = sample_df %>%
    mutate(.,
           indel_filter = samp_alignment_seq %>%
             subchar(., ctrl_indel, "_") %>%
             gsub("-", "", .) %>%
             strsplit(x = ., split = "") %>%
             .[[1]]
    ) %>%
    filter(indel_filter != "_") %>%
    dplyr::select(-indel_filter) %>%
    mutate(filtered_index = 1:NROW(max_base)) %>%
    mutate(ctrl_post_aligned_index = ctrl_df$post_aligned_index)

  sample_df = ctrl_df %>%
    dplyr::select(post_aligned_index, index) %>%
    dplyr::rename(ctrl_post_aligned_index = post_aligned_index, ctrl_index = index) %>%
    inner_join(., sample_df)

  sample_df %<>% mutate(sample_file = sample, ctrl_file = ctrl)
  ctrl_df %<>% mutate(sample_file = sample, ctrl_file = ctrl)

  trimmedData = list("sample_df" = sample_df,
                     "ctrl_df" = ctrl_df,
                     "ctrl_seq" = ctrl_seq,
                     "pre_aligned_sample_df" = pre_aligned_sample_df)

  ### alignment stats
  motif_alignment = matchPattern(pattern = DNAString(motif),
                                 subject = DNAString(trimmedData$ctrl_seq),
                                 fixed = FALSE)

  n_alignments = motif_alignment@ranges %>% length()

  motif_positions = mapply(FUN = seq,
                           from = motif_alignment@ranges@start,
                           to = (motif_alignment@ranges@start + nchar(motif) - 1)) %>%
    as.vector()

  names(motif_positions) = rep(x = c(1:n_alignments), each = nchar(motif))

  motif_positions = data.frame(ctrl_post_aligned_index = motif_positions,
                               motif_id = names(motif_positions))

  sample_df = trimmedData$sample_df %>%
    mutate(ctrl_max_base = trimmedData$ctrl_df$max_base,
           ctrl_base_call = trimmedData$ctrl_df$base_call
    ) %>%
    # Add a column for the percent base for the ctrl basecall
    mutate(ctrl_max_base_perc = {ifelse(ctrl_max_base == "A", A_perc,
                                        ifelse(ctrl_max_base == "C", C_perc,
                                               ifelse(ctrl_max_base == "G", G_perc,
                                                      ifelse(ctrl_max_base == "T", T_perc,0))))})

  sample_null = sample_df %>%
    filter(!(ctrl_post_aligned_index %in% motif_positions$ctrl_post_aligned_index))

  # This dataframe consists of all bases in the sample where the motif is found in the ctrl sequence
  sample_alt = motif_positions %>%
    inner_join(., sample_df)

  # Find all potential events of significant noise
  filtered_sample_alt = sample_alt %>%
    filter(grepl(wt, ctrl_max_base)) %>% # Use the ctrl wt base for determining data
    dplyr::rename(A = A_area, C = C_area, G = G_area, `T` = T_area) %>%
    gather(base, height, A:`T`) %>%
    filter(grepl(edit, base)) # Filter out hypothetical mutations that are not of interest

  n_comparisons = NROW(filtered_sample_alt)

  if (adjust_p) {
    p_adjust = 1-((1-p_value)^(1/n_comparisons))
  } else {
    p_adjust = p_value
  }

  zaga_parameters = make_ZAGA_df(sample_null, p_adjust = p_adjust)
  critical_values = zaga_parameters$crit

  output_sample_alt = gbm_adjust(sample_alt, wt,
                                 edit, motif,
                                 sample, critical_values)

  # Calculate the EI_index for each base
  sanger_EI = output_sample_alt %>%
    dplyr::group_by(ctrl_max_base) %>% # technically shouldn't change anything as it's the same across samples
    dplyr::mutate(AEI_sanger = (sum(G_perc) / (sum(!! sym(paste0(wt, "_perc"))) + sum(G_perc))),
                  CEI_sanger = (sum(T_perc) / (sum(!! sym(paste0(wt, "_perc"))) + sum(T_perc))),
                  GEI_sanger = (sum(A_perc) / (sum(!! sym(paste0(wt, "_perc"))) + sum(A_perc))),
                  TEI_sanger = (sum(C_perc) / (sum(!! sym(paste0(wt, "_perc"))) + sum(C_perc)))
    ) %>%
    dplyr::select(AEI_sanger:TEI_sanger) %>%
    dplyr::distinct() %>%
    ungroup() %>%
    # Keep only the EI_index for each reference base
    gather(EI_base, EI_sanger, AEI_sanger:TEI_sanger) %>%
    mutate(EI_base = gsub("EI_sanger", "", EI_base)) %>%
    dplyr::select(EI_base, EI_sanger)  %>%
    dplyr::rename(Edit =  EI_base, MEI = EI_sanger) %>%
    filter(grepl(wt, Edit)) %>%
    distinct()

  statisticalAnalysis = list("sample_df" = sample_df,
                             "sample_alt" = sample_alt,
                             "output_sample_alt" = output_sample_alt,
                             "zaga_parameters" = zaga_parameters,
                             "motif_positions" = motif_positions,
                             "sanger_EI" =  sanger_EI)

  editingData = calculateEditingData(statisticalAnalysis$output_sample_alt,
                                     edit)

  ### outputs
  list(
    sample_file = sample,
    raw_data = (statisticalAnalysis$output_sample_alt %>%
                  mutate(A = round(A_perc*100),
                         C = round(C_perc*100),
                         G = round(G_perc*100),
                         `T` = round(T_perc*100)) %>%
                  dplyr::select(index, ctrl_index, max_base, A:`T`, A_sig:T_sig,
                                motif, sample_file)),
    raw_sample_plot = (plotRawSample(sample_df_orig,
                                     statisticalAnalysis$sample_alt,
                                     filteredData$pre_cross_align_sample_df)),
    trimmed_sample_plot = (plotTrimmedSample(statisticalAnalysis$sample_df,
                                             filteredData$pre_cross_align_sample_df,
                                             statisticalAnalysis$output_sample_alt,
                                             sample_df_orig,
                                             statisticalAnalysis$sample_alt)),
    editing_plot = plotEditingData(editingData),
    editing_table = tableEditingData(editingData),
    mei_table = statisticalAnalysis$sanger_EI,
    stats_table = (statisticalAnalysis$zaga_parameters %>%
      filter(grepl(edit, Base)))
  )
})

map_dfr(res, \(x) mutate(x[['raw_data']],
                         sample_file = x[['sample_file']],
                         .before = 1)) %>%
  write_tsv('raw_edits.tsv')

map_dfr(res, \(x) mutate(x[['editing_table']],
                         sample_file = x[['sample_file']],
                         .before = 1)) %>%
  write_tsv('editing_table.tsv')

map_dfr(res, \(x) mutate(x[['mei_table']],
                         sample_file = x[['sample_file']],
                         .before = 1)) %>%
  write_tsv('mei_table.tsv')

map_dfr(res, \(x) mutate(x[['stats_table']],
                         sample_file = x[['sample_file']],
                         .before = 1)) %>%
  write_tsv('stats_table.tsv')

pdf('raw_sample_plot.pdf', h = 5, w = 12)
map(res, \(x) x[['raw_sample_plot']] + ggtitle(x[['sample_file']]))
dev.off()

pdf('trimmed_sample_plot.pdf', h = 5, w = 12)
map(res, \(x) x[['trimmed_sample_plot']] + ggtitle(x[['sample_file']]))
dev.off()

pdf('editing_plot.pdf', h = 5, w = 5)
map(res, \(x) x[['editing_plot']] + ggtitle(x[['sample_file']]))
dev.off()
