"""Download Axiom metadata.

Download all Axiom platemap and biochem metadata.

"""  # noqa: CPY001, INP001

import re
import os

from sh import aws
from pathlib import Path


def main() -> None:
    """Download metadata.

    Read in index file, download data.

    """
    aws_path = "s3://cellpainting-gallery/cpg0037-oasis/axiom/workspace/metadata"
    local_meta = "../1_snakemake/inputs/metadata"
    batches = ["prod_25", "prod_26", "prod_27", "prod_30"]

    os.makedirs(f"{local_meta}/metadata", exist_ok=True)
    os.makedirs(f"{local_meta}/biochem", exist_ok=True)

    # get metadata (both biochem.parquet and metadata.parquet)
    for batch in batches:
        batch_path = f"{aws_path}/{batch}/"
        aws_output = aws("s3", "ls", batch_path)
        plates = re.findall(r"plate_\d{8}", aws_output)

        for plate in plates:
            metadata_file = Path(f"{local_meta}/metadata/metadata_{plate}.parquet")
            biochem_file = Path(f"{local_meta}/biochem/biochem_{plate}.parquet")

            # only download if the files doesn't already exist
            if not metadata_file.exists():
                aws("s3", "cp", f"{batch_path}{plate}/metadata.parquet", str(metadata_file))

            if not biochem_file.exists():
                aws("s3", "cp", f"{batch_path}{plate}/biochem.parquet", str(biochem_file))


if __name__ == "__main__":
    main()
