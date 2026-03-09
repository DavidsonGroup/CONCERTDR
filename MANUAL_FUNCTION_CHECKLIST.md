# CONCERTDR Manual Function Check Checklist

Use this checklist to manually validate each exported function for Bioconductor readiness.

## Per-function checks

For each function below, verify all items:

- [ ] `?function_name` opens and documentation is clear
- [ ] `\value{}` matches actual returned object
- [ ] At least one runnable example exists outside `\donttest{}`
- [ ] Heavy examples (large files/interactive/network) are inside `\donttest{}`
- [ ] Example code runs without editing when package data is available
- [ ] Error messages are informative for common bad inputs

## Exported function list

### Configuration / setup
- [ ] `create_cmap_config_template`
- [ ] `create_console_config`
- [ ] `create_interactive_config`
- [ ] `interactive_cmap_setup`
- [ ] `read_cmap_config`

### Data extraction / processing
- [ ] `extract_cmap_data_from_config`
- [ ] `extract_cmap_data_from_siginfo`
- [ ] `extract_cmap_parameters`
- [ ] `subset_siginfo_beta`
- [ ] `process_combinations`
- [ ] `process_combinations_file`

### Signature matching / workflow
- [ ] `process_signature_with_df`
- [ ] `run_cmap_workflow`
- [ ] `demonstrate_workflow`

### Result annotation / utilities
- [ ] `annotate_drug_results`
- [ ] `extract_compound_id`
- [ ] `fuzzy_drug_match`

### Visualization
- [ ] `extract_signature_zscores`
- [ ] `plot_signature_direction_tile_barcode`

### Combination generation
- [ ] `generate_combinations_from_config`
- [ ] `generate_combinations_from_selections`

## Sign-off

- [ ] All exported functions manually checked
- [ ] All examples reviewed and categorized correctly (runnable vs `\donttest{}`)
- [ ] Re-ran package checks after final edits
