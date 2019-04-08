library(tidyverse)
library(kohonen)
library(ggthemes)
library(ggdendro)
library(viridis)

theme_default = theme_light() +
    theme(text = element_text(size=12, color="black"),
          axis.text = element_text(size=10, color="black"),
          plot.title = element_text(size=12),
          strip.background = element_blank(),
          strip.text = element_text(color="black"))

build_som_input = function(df, condition, left, right){
    df %>%
        filter(between(position, left, right)) %>%
        group_by(index) %>%
        mutate(signal = (signal-mean(signal))/sd(signal)) %>%
        filter(group == condition) %>%
        ungroup() %>%
        spread(position, signal) %>%
        select(-c(group, index)) %>%
        as.matrix() %>%
        return()
}

get_codes_df = function(results, condition="WT"){
    results %>%
        getCodes(idx=ifelse(condition=="WT", 1, 2)) %>%
        as_tibble() %>%
        rowid_to_column(var="unit") %>%
        gather(key = position, value=signal, -unit, convert=TRUE) %>%
        mutate(group=condition) %>%
        return()
}


main = function(mnase_data, annotation, sample_list, som_width, som_height, left_limit, right_limit,
                epochs, n_clusters, outputs, meta_data_out, seed,
                loss_out, dendrogram_out, code_vector_plot_out, data_by_unit_plot_out,
                data_by_unit_and_cluster_plot_out, data_by_cluster_plot_out, data_by_cluster_group_plot_out,
                heatmap_out){
    set.seed(seed)

    df = read_tsv(mnase_data,
                  col_names = c('group', 'sample', 'annotation', 'index', 'position', 'signal')) %>%
        filter(sample %in% sample_list) %>%
        group_by(group, index, position) %>%
        summarise(signal = mean(signal)) %>%
        ungroup()

    #convert dataframes into one matrix each for wild-type and mutant for SOM input
    wt_input = build_som_input(df, condition="WT-37C", left=left_limit, right=right_limit)
    mutant_input = build_som_input(df, condition="spt6-1004-37C", left=left_limit, right=right_limit)

    #train SOM
    som_results = xyf(wt_input, mutant_input,
                      rlen=epochs,
                      grid = somgrid(xdim = som_width, ydim=som_height,
                                     topo = "rectangular"),
                      user.weights=1,
                      normalizeDataLayers=FALSE,
                      mode="batch")

    #assess SOM training
    loss_df = som_results[["changes"]] %>%
        as_tibble() %>%
        magrittr::set_colnames(c("WT", "mutant")) %>%
        rowid_to_column(var="epoch") %>%
        gather(key=condition, value=loss, -epoch)

    loss_plot = ggplot(data = loss_df,
                       aes(x=epoch, y=loss, color=condition)) +
        geom_line() +
        scale_color_ptol() +
        scale_x_continuous(expand = c(0,0),
                           name = "training epoch") +
        ylab("mean average deviations from code vectors") +
        ggtitle("SOM loss during training") +
        theme_default +
        theme(legend.title = element_blank())
    ggsave(loss_out, plot=loss_plot, width=12, height=8, units="cm")


    # hierarchically cluster code vectors
    codes_hclust = som_results %>%
        object.distances(type="codes") %>%
        hclust(method="ward.D2")

    codes_dendrogram = ggdendrogram(codes_hclust) +
        ggtitle("hierarchical clustering of SOM code vectors (Ward criterion)") +
        theme(plot.title = element_text(size=12))
    ggsave(dendrogram_out, plot=codes_dendrogram, width=12, height=6, units="cm")

    codes_cluster_df = codes_hclust %>%
        cutree(k=n_clusters) %>%
        as_tibble() %>%
        magrittr::set_colnames("cluster") %>%
        rowid_to_column(var="unit")

    #df of code vector information, with cluster assignments
    codes_df = get_codes_df(som_results, condition="WT") %>%
        bind_rows(get_codes_df(som_results, condition="mutant")) %>%
        left_join(som_results[["unit.classif"]] %>%
                      table(dnn="unit") %>%
                      as_tibble() %>%
                      mutate(unit = as.integer(unit)),
                  by="unit") %>%
        left_join(codes_cluster_df, by="unit") %>%
        mutate(group = ordered(group, levels = c("WT", "mutant")))

    #df of data assignments to SOM units and clusters
    som_classifications = som_results[["unit.classif"]] %>%
        as_tibble() %>%
        magrittr::set_colnames("unit") %>%
        rowid_to_column(var="index") %>%
        mutate(unit = as.integer(unit)) %>%
        left_join(codes_cluster_df, by="unit")

    #df with number of datapoints in each unit/cluster (for plot labeling)
    unit_numbers_df = codes_df %>%
        group_by(unit) %>%
        summarise(n = first(n),
                  cluster = first(cluster)) %>%
        mutate(cluster = factor(cluster))

    #plot the code vectors
    code_vector_plot = ggplot() +
        geom_rect(data = unit_numbers_df,
                  aes(alpha=n),
                  xmin=left_limit, xmax=right_limit,
                  ymin=min(codes_df[["signal"]])*1.05,
                  ymax=max(codes_df[["signal"]])*1.05,
                  fill="grey") +
        geom_text(data = unit_numbers_df,
                  aes(label = n),
                  x = right_limit, y=max(codes_df[["signal"]]),
                  hjust=1, vjust=1) +
        geom_vline(xintercept = 0,
                   color="grey65",
                   size=0.5) +
        geom_line(data = codes_df,
                  aes(x=position, y=signal, color=group)) +
        scale_x_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=3),
                           labels = function(x)if_else(x==0, "TSS", as.character(x))) +
        scale_y_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=2),
                           name = "standard scores") +
        scale_alpha(guide=FALSE) +
        facet_wrap(~unit, ncol=som_width) +
        scale_color_ptol() +
        ggtitle("SOM code vectors") +
        theme_default +
        theme(strip.background = element_blank(),
              strip.text = element_blank(),
              panel.spacing = unit(5, "pt"),
              legend.title = element_blank(),
              axis.title.x = element_blank())
    ggsave(code_vector_plot_out, plot=code_vector_plot, width=24, height=12, units="cm")

    #plot the actual data facetted by which SOM unit they were
    #assigned to
    meta_df = df %>%
        filter(position %>% between(left_limit, right_limit)) %>%
        left_join(som_classifications, by="index") %>%
        group_by(group, position, unit, cluster) %>%
        summarise(mid = median(signal),
                  low = quantile(signal, 0.1),
                  high = quantile(signal, 0.9),
                  n = n()) %>%
        ungroup() %>%
        mutate(group = ordered(group, levels = c("WT-37C", "spt6-1004-37C"),
                               labels = c("WT", "spt6-1004"))) %>%
        write_tsv(meta_data_out)

    data_by_unit_plot = ggplot() +
        geom_rect(data = unit_numbers_df,
                  aes(alpha=n),
                  xmin=left_limit, xmax=right_limit,
                  ymin=0, ymax=max(meta_df[["high"]])*1.05,
                  fill="grey") +
        geom_vline(xintercept = 0,
                   color="grey65",
                   size=0.3) +
        geom_text(data = unit_numbers_df,
                  aes(label = n),
                  x = right_limit, y=max(meta_df[["high"]]),
                  hjust=1, vjust=1) +
        geom_ribbon(data = meta_df,
                    aes(x=position, ymin=low, ymax=high,
                        fill=group),
                    alpha=0.4, size=0) +
        geom_line(data = meta_df,
                  aes(x=position, y=mid,
                      color=group)) +
        facet_wrap(~unit, ncol=som_width) +
        scale_x_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=3),
                           labels = function(x)if_else(x==0, "TSS", as.character(x))) +
        scale_y_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=2),
                           name = "normalized counts") +
        scale_alpha(guide=FALSE) +
        scale_fill_ptol() +
        scale_color_ptol() +
        ggtitle("MNase dyad signal grouped by SOM unit") +
        theme_default +
        theme(strip.background = element_blank(),
              strip.text = element_blank(),
              panel.spacing = unit(5, "pt"),
              legend.title = element_blank(),
              axis.title.x = element_blank())
    ggsave(data_by_unit_plot_out, plot=data_by_unit_plot, width=24, height=12, units="cm")

    data_by_unit_and_cluster_plot = ggplot(data = meta_df,
                                           aes(x=position, ymin=low, ymax=high, y=mid,
                                               color=group, fill=group, group=interaction(group, unit))) +
        geom_vline(xintercept = 0) +
        geom_ribbon(alpha=0.2, size=0) +
        # geom_line(alpha=0.2) +
        scale_fill_ptol() +
        scale_color_ptol() +
        facet_grid(group~cluster) +
        scale_x_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=3),
                           labels = function(x)if_else(x==0, "TSS", as.character(x))) +
        scale_y_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=2),
                           name = "normalized counts") +
        ggtitle("MNase-seq dyad signal, grouped by SOM unit and facetted by cluster") +
        theme_default +
        theme(legend.title = element_blank())
    ggsave(data_by_unit_and_cluster_plot_out, data_by_unit_and_cluster_plot,
           width=16, height=10, units="cm")

    # plot the actual data facetted by which cluster
    # they were assigned to
    meta_cluster_df = df %>%
        left_join(som_classifications, by="index") %>%
        group_by(group, position, cluster) %>%
        summarise(mid = median(signal),
                  low = quantile(signal, 0.1),
                  high = quantile(signal, 0.9),
                  n = n()) %>%
        ungroup() %>%
        mutate(group = ordered(group, levels = c("WT-37C", "spt6-1004-37C"),
                               labels = c("WT", "spt6-1004")))

    data_by_cluster_plot = ggplot(data = meta_cluster_df,
                                  aes(x=position, ymin=low, ymax=high, y=mid,
                                      color=group, fill=group)) +
        geom_vline(xintercept = 0) +
        geom_ribbon(alpha=0.4, size=0) +
        geom_line(alpha=0.8) +
        scale_fill_ptol() +
        scale_color_ptol() +
        facet_wrap(~cluster, ncol=2) +
        scale_x_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=3),
                           labels = function(x)if_else(x==0, "TSS", as.character(x))) +
        scale_y_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=2),
                           name = "normalized counts") +
        ggtitle("MNase-seq dyad signal, facetted by cluster") +
        theme_default +
        theme(legend.title = element_blank())
    ggsave(data_by_cluster_plot_out, plot=data_by_cluster_plot, width=16, height=8, units="cm")

    data_by_cluster_group_plot = ggplot(data = meta_cluster_df,
           aes(x=position, ymin=low, ymax=high, y=mid,
               group = interaction(group, cluster),
               color=factor(cluster), fill=factor(cluster))) +
        geom_vline(xintercept = 0) +
        geom_ribbon(alpha=0.4, size=0) +
        geom_line(alpha=0.8) +
        scale_fill_tableau(name="cluster") +
        scale_color_tableau(name="cluster") +
        scale_x_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=3),
                           labels = function(x)if_else(x==0, "TSS", as.character(x))) +
        scale_y_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=2),
                           name = "normalized counts") +
        ggtitle("MNase-seq dyad signal, grouped by cluster") +
        facet_wrap(~group, ncol=2) +
        theme_default
    ggsave(data_by_cluster_group_plot_out, data_by_cluster_group_plot, width=16, height=8, units="cm")

    #plot heatmap
    df_heatmap = som_classifications %>%
        group_by(cluster) %>%
        arrange(unit, index) %>%
        mutate(new_index = row_number()) %>%
        left_join(df, by="index")

    heatmap = ggplot(data = df_heatmap, aes(x=position, y=new_index, fill=signal)) +
        geom_raster() +
        geom_vline(xintercept = c(left_limit, 0, right_limit)) +
        scale_fill_viridis(limits = c(NA, quantile(df$signal, 0.9)),
                           oob=scales::squish) +
        facet_grid(cluster~group, space="free_y", scales="free_y") +
        scale_x_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=3),
                           labels = function(x)if_else(x==0, "TSS", as.character(x))) +
        scale_y_continuous(expand = c(0,0)) +
        theme_default +
        theme(axis.title = element_blank())
    ggsave(heatmap_out, plot=heatmap, width=16, height=28, units="cm")

    #write out annotations
    bed = read_tsv(annotation, col_names = c('chrom', 'start', 'end', 'name', 'score', 'strand')) %>%
        rowid_to_column(var="index") %>%
        left_join(som_classifications, by=c("index")) %>%
        arrange(cluster, unit, score)

    for (i in 1:n_clusters){
        bed %>%
            filter(cluster==i) %>%
            select(-c(index, unit, cluster)) %>%
            write_tsv(outputs[i], col_names=FALSE)
    }
}

main(mnase_data = snakemake@input[["mnase_data"]],
     annotation = snakemake@input[["annotation"]],
     sample_list = snakemake@params[["sample_list"]],
     som_width = 8,
     som_height = 6,
     left_limit = snakemake@params[["left_limit"]],
     right_limit = snakemake@params[["right_limit"]],
     epochs = snakemake@params[["epochs"]],
     n_clusters = 3,
     outputs = snakemake@output[["annotations"]],
     meta_data_out = snakemake@output[["metagene_data"]],
     seed = 1,
     loss_out = snakemake@output[["loss_plot"]],
     dendrogram_out = snakemake@output[["dendrogram"]],
     code_vector_plot_out = snakemake@output[["code_vectors"]],
     data_by_unit_plot_out = snakemake@output[["data_by_unit"]],
     data_by_unit_and_cluster_plot_out = snakemake@output[["data_by_unit_and_cluster"]],
     data_by_cluster_plot_out = snakemake@output[["data_by_cluster"]],
     data_by_cluster_group_plot_out = snakemake@output[["data_by_cluster_group"]],
     heatmap_out = snakemake@output[["heatmap"]])

