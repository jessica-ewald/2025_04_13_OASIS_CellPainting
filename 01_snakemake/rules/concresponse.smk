features = config["features"]
scenario = config["workflow"]
name = config["name"]

rule compute_distances_R:
    input:
        f"outputs/{features}/{name}/profiles/{scenario}.parquet",
    output:
        expand("outputs/{features}/{name}/distances/{method}.parquet", method=config["distances_R"], features=config["features"], scenario=config["workflow"], name=config["name"]),
    params:
        cover_var=config["cover_var"],
        treatment=config["treatment"],
        categories=",".join(config["categories"]),
        distances=config["distances_R"],
    shell:
        """
        for method in {params.distances}; do
            method_name=$(echo $method | tr -d '[],"') 
            Rscript concresponse/compute_distances.R {input} outputs/{features}/{name}/distances/${{method_name}}.parquet {params.cover_var} {params.treatment} {params.categories} ${{method_name}}
        done
        """

rule compute_distances_python:
    input:
        f"outputs/{features}/{name}/profiles/{scenario}.parquet",
    output:
        expand("outputs/{features}/{name}/distances/{method}.parquet", method=config["distances_python"], features=config["features"], name=config["name"]),
    params:
        distances=config["distances_python"],
    run:
        for method in config["distances_python"]:
            output_file = f"outputs/{features}/{name}/distances/{method}.parquet"
            cr.ap.calculate_distances(input[0], output_file, method)


distances = config["distances_R"] + config["distances_python"]
rule compile_distances:
    input:
        [f"outputs/{features}/{name}/distances/{method}.parquet" for method in distances],
    output:
        f"outputs/{features}/{name}/distances/distances.parquet",
    params:
        transform=config["dist_transform"],
    run:
        input_files = list(input)
        cr.compile_dist.compile_dist(input_files, params.transform, *output)


rule fit_curves:
    input:
        f"outputs/{features}/{name}/distances/distances.parquet",
    output:
        f"outputs/{features}/{name}/curves/bmds.parquet",
    params:
        num_sds = config['num_sds']
    shell:
        "Rscript concresponse/fit_curves.R {input} {output} {params.num_sds}"

rule fit_curves_cc:
    input:
        f"outputs/{features}/{name}/profiles/{scenario}.parquet",
    output:
        f"outputs/{features}/{name}/curves/ccpods.parquet",
    params:
        num_sds = config['num_sds'],
        meta_nm = "Metadata_Count_Cells"
    shell:
        "Rscript concresponse/fit_curves_meta.R {input} {output} {params.num_sds} {params.meta_nm}"

rule fit_curves_mtt:
    input:
        f"outputs/{features}/{name}/profiles/{scenario}.parquet",
    output:
        f"outputs/{features}/{name}/curves/mttpods.parquet",
    params:
        num_sds = config['num_sds'],
        meta_nm = "Metadata_mtt_normalized"
    shell:
        "Rscript concresponse/fit_curves_meta.R {input} {output} {params.num_sds} {params.meta_nm}"

rule fit_curves_ldh:
    input:
        f"outputs/{features}/{name}/profiles/{scenario}.parquet",
    output:
        f"outputs/{features}/{name}/curves/ldhpods.parquet",
    params:
        num_sds = config['num_sds'],
        meta_nm = "Metadata_ldh_abs_signal"
    shell:
        "Rscript concresponse/fit_curves_meta.R {input} {output} {params.num_sds} {params.meta_nm}"

rule select_pod:
    input:
        f"outputs/{features}/{name}/curves/bmds.parquet",
        f"outputs/{features}/{name}/curves/ccpods.parquet",
    output:
        f"outputs/{features}/{name}/curves/pods.parquet",
    shell:
        "Rscript concresponse/select_pod.R {input} {output}"


rule plot_cc_curve_fits:
    input:
        f"outputs/{features}/{name}/curves/ccpods.parquet",
        f"outputs/{features}/{name}/profiles/{scenario}.parquet",
    output:
        f"outputs/{features}/{name}/curves/plots/cc_plots.pdf",
    params:
        meta_nm = "Metadata_Count_Cells"
    shell:
        "Rscript concresponse/plot_meta_curve.R {input} {output} {params.meta_nm}"

rule plot_mtt_curve_fits:
    input:
        f"outputs/{features}/{name}/curves/mttpods.parquet",
        f"outputs/{features}/{name}/profiles/{scenario}.parquet",
    output:
        f"outputs/{features}/{name}/curves/plots/mtt_plots.pdf",
    params:
        meta_nm = "Metadata_mtt_normalized"
    shell:
        "Rscript concresponse/plot_meta_curve.R {input} {output} {params.meta_nm}"

rule plot_ldh_curve_fits:
    input:
        f"outputs/{features}/{name}/curves/ldhpods.parquet",
        f"outputs/{features}/{name}/profiles/{scenario}.parquet",
    output:
        f"outputs/{features}/{name}/curves/plots/ldh_plots.pdf",
    params:
        meta_nm = "Metadata_ldh_abs_signal"
    shell:
        "Rscript concresponse/plot_meta_curve.R {input} {output} {params.meta_nm}"


rule plot_cp_curve_fits:
    input:
        f"outputs/{features}/{name}/curves/pods.parquet",
        f"outputs/{features}/{name}/curves/ccpods.parquet",
        f"outputs/{features}/{name}/distances/distances.parquet",
    output:
        f"outputs/{features}/{name}/curves/plots/cp_plots.pdf",
    shell:
        "Rscript concresponse/plot_cp_curve.R {input} {output}"