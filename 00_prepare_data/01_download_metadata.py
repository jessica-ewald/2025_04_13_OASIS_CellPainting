"""Download Axiom metadata.

Download all Axiom platemap and biochem metadata.

"""  # noqa: CPY001, INP001

import re
from pathlib import Path

from sh import aws


def main() -> None:
    """Download metadata.

    Read in index file, download data.
    """
    aws_path = "s3://cellpainting-gallery/cpg0037-oasis/broad/workspace/metadata"

    metadata_dir = Path("../01_snakemake/inputs/metadata/plates/")
    metadata_dir.mkdir(parents=True, exist_ok=True)

    batches = [
        "2025_01_16_U2OS_Batch1",
        "2025_01_12_U2OS_Batch2",
        "2025_01_27_U2OS_Batch3",
        "2025_01_29_U2OS_Batch4",
    ]

    # get metadata (both biochem.parquet and metadata.parquet)
    for batch in batches:
        batch_path = f"{aws_path}/{batch}/platemap/"
        aws_output = aws("s3", "ls", batch_path)
        plates = re.findall(r"BR\d{8}.txt", aws_output)

        for plate in plates:
            metadata_file = metadata_dir / plate

            # only download if the files don't already exist
            if not metadata_file.exists():
                aws("s3", "cp", f"{batch_path}{plate}", str(metadata_file))


if __name__ == "__main__":
    main()
