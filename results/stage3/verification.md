# Stage-3 loader verification

data_version: `sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74`

**42/42 checks passed.**

| Check | Result | Detail |
|---|---|---|
| protein_raw_row_count | PASS | expected 4344, got 4344 |
| peptide_raw_row_count | PASS | expected 33358, got 33358 |
| protein_header_n_intensity_columns | PASS | expected 8, got 8 |
| peptide_header_n_intensity_and_detection_columns | PASS | expected 8/8, got 8/8 |
| pairing_n_samples | PASS | expected 8, got 8 |
| pairing_condition_4_4 | PASS | condition_counts={'raloxifene-d0': 4, 'control': 4} |
| pairing_batch_6_2 | PASS | batch_counts={'B2021-05-06': 6, 'B2022-03-18': 2} |
| pairing_experimental_8_control_0 | PASS | no QC/pool/reference/blank marker in any of the 8 Replicate filenames |
| pairing_matches_samples_tsv_row_by_row | PASS | independent re-derivation == samples.tsv for every sample_id (search_scan_file_id, condition, batch) |
| pairing_search_id_has_data_columns | PASS | every independently-derived search_scan_file_id has an Intensity column in both quant files |
| protein_orientation_rows_are_samples | PASS | abundances.shape=(8, 4344), expected (8, 4344) |
| protein_abundances_dtype_float | PASS | dtype=float64 |
| protein_feature_names_dtype_str | PASS | dtype=<U39 |
| protein_finite_values_all_positive | PASS | n_finite=21654, n_nonpositive=0 |
| protein_finite_values_all_finite | PASS | no inf among finite-flagged values |
| peptide_orientation_rows_are_samples | PASS | abundances.shape=(8, 33358), expected (8, 33358) |
| peptide_abundances_dtype_float | PASS | dtype=float64 |
| peptide_feature_names_dtype_str | PASS | dtype=<U64 |
| peptide_finite_values_all_positive | PASS | n_finite=177683, n_nonpositive=0 |
| peptide_finite_values_all_finite | PASS | no inf among finite-flagged values |
| protein_identifier_roundtrip | PASS | feature_names == raw column order/text exactly |
| peptide_identifier_roundtrip | PASS | feature_names == raw column order/text exactly |
| protein_nan_token_count | PASS | expected 14, got 14 |
| protein_zero_token_count | PASS | expected 13084, got 13084 |
| peptide_zero_count | PASS | expected 89181, got 89181 |
| peptide_detection_type_tally | PASS | counts={'MSMS': 152578, 'MBR': 25105, 'NotDetected': 75363, 'MSMSIdentifiedButNotQuantified': 10135, 'MSMSAmbiguousPeakfinding': 3683}, unknown=0, expected={'MSMS': 152578, 'MBR': 25105, 'NotDetected': 75363, 'MSMSIdentifiedButNotQuantified': 10135, 'MSMSAmbiguousPeakfinding': 3683}, taxonomy_agrees_with_loader=True |
| protein_contaminant_count | PASS | expected 33, got 33 |
| protein_complete_noncontaminant_count | PASS | expected 1801, got 1801 |
| peptide_contaminant_count | PASS | expected 449, got 449 |
| peptide_complete_noncontaminant_count | PASS | expected 10229, got 10229 |
| protein_spot_reconciliation | PASS | checked=250, mismatches=0 |
| peptide_spot_reconciliation | PASS | checked=250, mismatches=0 |
| protein_two_way_total_median | PASS | raw-parse and Dataset totals/medians agree for every sample |
| peptide_two_way_total_median | PASS | raw-parse and Dataset totals/medians agree for every sample |
| limelight_label_map_reproduced | PASS | resolved 8 label(s), all matched results/stage2/limelight_label_map.tsv |
| psms_nsaf_n_group_rows | PASS | expected 4344, got 4344 |
| psms_nsaf_n_comma_joined_groups | PASS | expected 27, got 27 |
| psms_nsaf_sample_ids_match_samples_tsv_order | PASS | limelight.sample_ids=[np.str_('AZ905'), np.str_('AZ907'), np.str_('AZ909'), np.str_('AZ906'), np.str_('AZ908'), np.str_('AZ910'), np.str_('AZ941'), np.str_('AZ942')] vs samples.tsv order ['AZ905', 'AZ907', 'AZ909', 'AZ906', 'AZ908', 'AZ910', 'AZ941', 'AZ942'] |
| psms_full_reconciliation | PASS | checked=34752 (expected 34752), mismatches=0 |
| nsaf_full_reconciliation | PASS | checked=34752 (expected 34752), mismatches=0 |
| nsaf_column_sums_in_range | PASS | sums={'1_0': 0.9862, '1_05': 0.9838, '1_050': 0.9835, '1_0506': 0.9885, '1_0506_': 0.9932, '1_0506_A': 0.9916, '2_0318_': 0.9907, '2_0318_A': 0.9821} |
| data_version_matches_metadata_stamp | PASS | computed=sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74, expected (results/metadata/data_version.json)=sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74 |
