#!/usr/bin/env python3
"""Build the GitHub/Zenodo release archive without redistributing raw GEO data."""
from pathlib import Path
import hashlib
import zipfile

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT.parent / "release"
ARCHIVE = RELEASE / "ich-hfcg-raw-reanalysis-v1.0.0.zip"
EXCLUDED_DIRS = {
    ".git", ".renv-cache", ".venv", "data/raw", "data/extracted",
    "r-library", "results/objects", "vendor", "release",
}
EXCLUDED_SUFFIXES = {".part", ".tbi"}


def include(path: Path) -> bool:
    rel = path.relative_to(ROOT).as_posix()
    if any(rel == item or rel.startswith(item + "/") for item in EXCLUDED_DIRS):
        return False
    if path.suffix in EXCLUDED_SUFFIXES or ".parts" in path.parts:
        return False
    return path.is_file()


def add(zf: zipfile.ZipFile, source: Path, arcname: str) -> None:
    zf.write(source, arcname=arcname, compress_type=zipfile.ZIP_DEFLATED,
             compresslevel=9)


def main() -> None:
    RELEASE.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(ARCHIVE, "w", allowZip64=True) as zf:
        for path in sorted(ROOT.rglob("*")):
            if include(path):
                add(zf, path, f"ich-hfcg-raw-reanalysis/{path.relative_to(ROOT).as_posix()}")
        submission_root = ROOT.parent / "outputs" / "Frontiers_submission_package"
        submission_names = [
            "Manuscript_Frontiers_Inflammation_revised_raw_reanalysis.docx",
            "Manuscript_Frontiers_Inflammation_revised_raw_reanalysis.pdf",
            "Cover_Letter_Frontiers_Inflammation_revised.docx",
            "Cover_Letter_Frontiers_Inflammation_revised.pdf",
            "Scope_Statement_revised.docx",
            "Scope_Statement_revised.pdf",
            "Figure_Alt_Text_revised.docx",
            "Figure_Alt_Text_revised.pdf",
            "Supplementary_Material_Figures_and_Table_Legends_revised.docx",
            "Supplementary_Material_Figures_and_Table_Legends_revised.pdf",
            "Main_Figures_revised.zip",
            "Supplementary_Figures_revised.zip",
            "Supplementary_Tables_revised.zip",
        ]
        for name in submission_names:
            path = submission_root / name
            if path.exists():
                add(zf, path, f"ich-hfcg-raw-reanalysis/submission/{name}")
    digest = hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
    checksum = ARCHIVE.with_suffix(ARCHIVE.suffix + ".sha256")
    checksum.write_text(f"{digest}  {ARCHIVE.name}\n", encoding="ascii")
    print(ARCHIVE)
    print(checksum)


if __name__ == "__main__":
    main()
