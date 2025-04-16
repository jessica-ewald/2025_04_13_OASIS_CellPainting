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