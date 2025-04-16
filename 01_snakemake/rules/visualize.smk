features = config["features"]
scenario = config["workflow"]
name = config["name"]

rule make_umaps:
    input:
        f"outputs/{features}/{name}/profiles/{scenario}.parquet",
        f"outputs/{features}/{name}/curves/pods.parquet",
        f"outputs/{features}/{name}/curves/ccpods.parquet",
    output:
        f"outputs/{features}/{name}/figures/umaps.pdf",
    run:
        vs.umaps.make_umaps(*input, *output)