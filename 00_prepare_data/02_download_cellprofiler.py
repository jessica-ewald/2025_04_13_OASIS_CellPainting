"""Download CellProfiler features.

Use AWS cli to download all CellProfiler profiles.

"""  # noqa: CPY001, INP001

import re

from sh import aws
from pathlib import Path


def main() -> None:
    """Download data.

    Read in index file, download data.

    """
    aws_path = "s3://cellpainting-gallery/cpg0037-oasis/broad/workspace/profiles"
    batches = [
        "2025_01_16_U2OS_Batch1",
        "2025_01_12_U2OS_Batch2",
        "2025_01_27_U2OS_Batch3",
        "2025_01_29_U2OS_Batch4",
    ]

    prof_dir = Path("../01_snakemake/inputs/profiles/cellprofiler/plates")
    prof_dir.mkdir(parents=True, exist_ok=True)

    # get Dino embedding paths
    for batch in batches:
        batch_path = f"{aws_path}/{batch}/"
        aws_output = aws("s3", "ls", batch_path)
        plates = re.findall(r"BR\d{8}", aws_output)

        for plate in plates:
            aws(
                "s3",
                "cp",
                f"{batch_path}{plate}/{plate}.csv.gz",
                f"{prof_dir}/{plate}.csv.gz",
            )


if __name__ == "__main__":
    main()
