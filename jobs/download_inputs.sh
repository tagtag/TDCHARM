#!/bin/bash
# ============================================================
# jobs/download_inputs.sh
#
# Downloads the public input files into the working directory.
# Run on a node with internet access (login node), not as a job.
#
# Usage:
#   GSE303006_WORKDIR=/path/to/workdir bash jobs/download_inputs.sh
#
# Files:
#   GEO GSE303006 supplementary files used by scripts 01-04
#   GENCODE vM23 annotation (mm10)
#   GSE303006_brain.20k.tdg.tar.gz -> tdg20k/geo_submit/brain_20k_tdg/
#     (reconstructed 20-kb 3D structures, *.20k.3dg.gz, one per cell)
# Already complete files are not downloaded again.
# ============================================================

set -euo pipefail
: "${GSE303006_WORKDIR:?Set GSE303006_WORKDIR}"
mkdir -p "${GSE303006_WORKDIR}"
cd "${GSE303006_WORKDIR}"

GEO="${TDCHARM_GEO_URL:-https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303006/suppl}"
GENCODE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23"

fetch() {  # fetch <url> <file>; skips files that are already complete
    local url="$1" f="$2"
    if [ -s "$f" ]; then
        echo "exists: $f"
    else
        echo "downloading: $f"
        curl -fL --retry 5 --retry-delay 10 -C - -o "$f.part" "$url"
        mv "$f.part" "$f"
    fi
}

for f in \
    GSE303006_brain.atac.fragment.tsv.gz \
    GSE303006_brain.h3k27ac.fragment.tsv.gz \
    GSE303006_brain.RNAcounts.gene.total.format.tsv.gz \
    GSE303006_charm_metadata_qcpass.tsv.gz
do
    fetch "${GEO}/${f}" "$f"
done
fetch "${GENCODE}/gencode.vM23.annotation.gtf.gz" gencode.vM23.annotation.gtf.gz

echo
echo "Checking gzip integrity..."
for f in GSE303006_brain.atac.fragment.tsv.gz GSE303006_brain.h3k27ac.fragment.tsv.gz \
         GSE303006_brain.RNAcounts.gene.total.format.tsv.gz \
         GSE303006_charm_metadata_qcpass.tsv.gz gencode.vM23.annotation.gtf.gz; do
    gzip -t "$f" && echo "ok: $f"
done

# ------------------------------------------------------------
# Reconstructed 20-kb 3D structures
# ------------------------------------------------------------
TDG_DIR="tdg20k/geo_submit/brain_20k_tdg"
TDG_TAR="GSE303006_brain.20k.tdg.tar.gz"

n_tdg() { find "$1/" -maxdepth 1 -name '*.20k.3dg.gz' 2>/dev/null | wc -l; }

if [ -d "${TDG_DIR}" ] && [ "$(n_tdg "${TDG_DIR}")" -gt 0 ]; then
    echo
    echo "3D directory present: ${TDG_DIR} ($(n_tdg "${TDG_DIR}") files)"
else
    fetch "${GEO}/${TDG_TAR}" "${TDG_TAR}"
    gzip -t "${TDG_TAR}" && echo "ok: ${TDG_TAR}"
    echo
    echo "First entries of ${TDG_TAR}:"
    tar -tzf "${TDG_TAR}" | head -5 || true
    mkdir -p tdg20k
    echo "Extracting ${TDG_TAR} into tdg20k/ ..."
    tar -xzf "${TDG_TAR}" -C tdg20k
    if [ "$(n_tdg "${TDG_DIR}")" -eq 0 ]; then
        # The archive layout differs from the expected path: link the
        # directory that holds the *.20k.3dg.gz files to the expected path.
        src="$(find "$(pwd)/tdg20k" -name '*.20k.3dg.gz' -print -quit | xargs -r dirname)"
        if [ -z "$src" ]; then
            echo "No *.20k.3dg.gz files found after extraction; check ${TDG_TAR}." >&2
            exit 1
        fi
        mkdir -p "$(dirname "${TDG_DIR}")"
        ln -sfn "$src" "${TDG_DIR}"
        echo "Linked ${TDG_DIR} -> ${src}"
    fi
    echo "3D directory ready: ${TDG_DIR} ($(n_tdg "${TDG_DIR}") files)"
fi

echo
echo "Inputs in $(pwd):"
ls -lh GSE303006_brain.atac.fragment.tsv.gz GSE303006_brain.h3k27ac.fragment.tsv.gz \
       GSE303006_brain.RNAcounts.gene.total.format.tsv.gz \
       GSE303006_charm_metadata_qcpass.tsv.gz gencode.vM23.annotation.gtf.gz
