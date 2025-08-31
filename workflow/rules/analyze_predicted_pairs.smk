
rule visualize_synergistic_pairs:
  input:
    synergy_predictions = "results/synergy_predictions_{cell_type}.tsv",
    hic = config["scratch_storage"] + "{cell_type}.hic"
  output:
    completed_file = config["scratch_storage"] + "{cell_type}_synergistic_pairs_hic_plots/" + "completed.txt"
  params:
    hic_res = config['thresholds']['hic_res'],
    hic_res_small = config['thresholds']['hic_res_small'],
    plot_window = config['thresholds']['plot_window'],
    plots_dir = config["scratch_storage"] + "{cell_type}_synergistic_pairs_hic_plots/"
  conda:
    "../envs/hic_contact.yml"
  resources:
    mem = "64G",
    time = "4:00:00"
  script:
    "../scripts/visualize_synergistic_pairs.py"
