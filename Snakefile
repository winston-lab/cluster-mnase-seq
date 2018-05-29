#!/usr/bin/env python

configfile: "config.yaml"

localrules:
    all,

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

rule all:
    input:
        "figures/mnase_data-heatmap.svg"

rule cluster_mnase_data:
    input:
        mnase_data = config["mnase_data"],
        annotation = config["annotation"],
    output:
        annotations = expand("spt6-1004-37C-induced-intragenic-MNase-cluster-{x}.bed", x=[1,2]),
        loss_plot = "figures/som_loss-plot.svg",
        dendrogram = "figures/som_dendrogram.svg",
        code_vectors = "figures/som_code-vectors.svg",
        data_by_unit = "figures/mnase_data-by-som-unit.svg",
        data_by_unit_and_cluster = "figures/mnase_data-by-som-unit-and-cluster.svg",
        data_by_cluster = "figures/mnase_data-by-som-cluster.svg",
        data_by_cluster_group = "figures/mnase_data-by-som-cluster-group.svg",
        heatmap = "figures/mnase_data-heatmap.svg"
    params:
        sample_list = config["samples"],
        left_limit = config["left_limit"],
        right_limit = config["right_limit"],
        epochs = config["epochs"]
    script:
        "scripts/cluster_mnase.R"

