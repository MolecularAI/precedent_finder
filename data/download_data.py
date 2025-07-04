"""Download pre-processed USPTO data from Zenodo."""
import requests
from pathlib import Path

from tqdm import tqdm


def download_file(url: str, file_path: Path) -> None:
    """Download a file in chunks.

    :param url: File download URL
    :param file_path: File save path
    """
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        total_size = int(response.headers.get("content-length", 0))
        pbar = tqdm(
            total=total_size,
            desc=file_path.name,
            unit="B",
            unit_scale=True,
        )
        with file_path.open(mode="wb") as fileobj:
            for chunk in response.iter_content(chunk_size=1024):
                fileobj.write(chunk)
                pbar.update(len(chunk))
        pbar.close()

    return


if __name__ == "__main__":
    datasets = [
        {
            "filename": "metadata_cleaned.csv.gz",
            "url": "https://zenodo.org/records/15798247/files/metadata_cleaned.csv.gz",
        },
        {
            "filename": "rfp_precomp.npz",
            "url": "https://zenodo.org/records/15798247/files/rfp_precomp.npz",
        },
        {
            "filename": "rcfp_precomp.npz",
            "url": "https://zenodo.org/records/15798247/files/rcfp_precomp.npz",
        },
    ]

    download_folder = Path(__file__).parent
    for dataset in datasets:
        download_file(
            url=dataset["url"],
            file_path=download_folder.joinpath(dataset["filename"]),
        )
