#!/usr/bin/env python3
"""Assemble revised Frontiers figure and supplementary-table ZIP files."""
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT.parent / "outputs" / "Frontiers_submission_package"


def zip_directory(source: Path, target: Path) -> None:
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for path in sorted(source.glob("*")):
            if path.is_file():
                zf.write(path, path.name)


def main() -> None:
    zip_directory(PACKAGE / "Figures_main", PACKAGE / "Main_Figures_revised.zip")
    zip_directory(PACKAGE / "Supplementary_figures", PACKAGE / "Supplementary_Figures_revised.zip")

    original = PACKAGE / "Supplementary_Tables_S1-S9.zip"
    revised = PACKAGE / "Supplementary_Tables_revised.zip"
    retained = {
        "Table_S1_HFCG_gene_list_sources.xlsx",
        "Table_S3_Enrichment_worksheets.xlsx",
        "Table_S7_Annotated_cell_type_markers.xlsx",
        "Table_S8_Myeloid_subcluster_markers.xlsx",
        "Table_S9_CellChat_interactions_and_pathways.xlsx",
    }
    with zipfile.ZipFile(original) as source, zipfile.ZipFile(
        revised, "w", zipfile.ZIP_DEFLATED, compresslevel=9
    ) as target:
        for name in sorted(retained):
            target.writestr(name, source.read(name))
        replacements = {
            ROOT / "results" / "tables" / "TableS2_S4_raw_GEO_reanalysis.xlsx":
                "Table_S2_S4_Raw_GEO_reanalysis.xlsx",
            ROOT / "results" / "tables" / "TableS5_S6_MR_colocalisation_complete.xlsx":
                "Table_S5_S6_MR_colocalisation_complete.xlsx",
            ROOT / "results" / "tables" / "R_package_versions.tsv":
                "Table_S10_R_package_versions.tsv",
        }
        for source_path, archive_name in replacements.items():
            if not source_path.exists():
                raise FileNotFoundError(source_path)
            target.write(source_path, archive_name)
    print(PACKAGE / "Main_Figures_revised.zip")
    print(PACKAGE / "Supplementary_Figures_revised.zip")
    print(revised)


if __name__ == "__main__":
    main()
