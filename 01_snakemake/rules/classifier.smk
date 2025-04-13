features = config["features"]
scenario = config["workflow"]
name = config["name"]

rule create_classifier_profiles:
    input:
        f"outputs/{features}/{name}/profiles/{scenario}.parquet",
        f"outputs/{features}/{name}/curves/pods.parquet",
    output:
        f"outputs/{features}/{name}/aggregated_profiles/agg.parquet",
    run:
        cl.aggregate_profiles.aggregate_profiles(*input, *output)


rule axiom_assay_hitcall:
    input:
        f"outputs/{features}/{name}/curves/mttpods.parquet",
        f"outputs/{features}/{name}/curves/ldhpods.parquet",
        f"outputs/{features}/{name}/curves/ccpods.parquet",
    output:
        f"outputs/{features}/{name}/curves/axiom_hits.parquet",
    run:
        cl.hitcalls.call_hits(*input, *output)


rule toxcast_cellbased_binary:
    input:
        f"outputs/{features}/{name}/aggregated_profiles/agg.parquet",
        f"inputs/annotations/toxcast_cellbased_binary.parquet",
    output:
        f"outputs/{features}/{name}/classifier_results/toxcast_cellbased_binary_predictions.parquet",
    run:
        cl.classify.predict_binary(*input, *output)


rule toxcast_cellfree_binary:
    input:
        f"outputs/{features}/{name}/aggregated_profiles/agg.parquet",
        f"inputs/annotations/toxcast_cellfree_binary.parquet",
    output:
        f"outputs/{features}/{name}/classifier_results/toxcast_cellfree_binary_predictions.parquet",
    run:
        cl.classify.predict_binary(*input, *output)


rule toxcast_cytotox_binary:
    input:
        f"outputs/{features}/{name}/aggregated_profiles/agg.parquet",
        f"inputs/annotations/toxcast_cytotox_binary.parquet",
    output:
        f"outputs/{features}/{name}/classifier_results/toxcast_cytotox_binary_predictions.parquet",
    run:
        cl.classify.predict_binary(*input, *output)


rule predict_axiom_binary:
    input:
        f"outputs/{features}/{name}/aggregated_profiles/agg.parquet",
        f"outputs/{features}/{name}/curves/axiom_hits.parquet",
    output:
        f"outputs/{features}/{name}/classifier_results/axiom_binary_predictions.parquet",
    run:
        cl.classify.predict_binary(*input, *output)


rule predict_axiom_continuous:
    input:
        f"outputs/{features}/{name}/profiles/{scenario}.parquet"
    output:
        f"outputs/{features}/{name}/classifier_results/axiom_continuous_predictions.parquet",
        f"outputs/{features}/{name}/classifier_results/axiom_continuous_metrics.parquet",
    run:
        cl.regression.predict_axiom_assays(*input, *output)