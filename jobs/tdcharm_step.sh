#!/bin/bash
# ============================================================
# jobs/tdcharm_step.sh
#
# Runs ONE step of the TDCHARM workflow. Submitted by
# jobs/submit_all.sh through qsub (SGE/UGE, PBS Pro or Torque), but
# can also be run directly with bash for testing.
#
# Required environment variables (passed by submit_all.sh with -v):
#   TDCHARM_STEP        synthetic | preprocessing | analysis | revision | script
#   TDCHARM_SCRIPT      file name under R/ (only for TDCHARM_STEP=script)
#   TDCHARM_REPO        repository directory (patched clone)
#   GSE303006_WORKDIR   working directory with the GEO input files
#                       (for "synthetic": the directory to create)
# Optional:
#   TDCHARM_R_SETUP     shell command run before R, e.g.
#                       "module load R/4.3.3" or
#                       "source ~/miniconda3/bin/activate r43"
#   TDCHARM_RSCRIPT     Rscript executable (default: Rscript)
#   any TDCHARM_* setting documented in run_revision.R
# ============================================================

set -euo pipefail

# SGE/UGE: -hold_jid releases dependent jobs even when this job fails,
# unless the job exits with status 100 (job goes to error state and
# its dependents keep waiting). Translate any failure into 100 there.
# PBS/Torque use "afterok", so the normal non-zero status is enough.

: "${TDCHARM_STEP:?TDCHARM_STEP is not set}"
: "${TDCHARM_REPO:?TDCHARM_REPO is not set}"
: "${GSE303006_WORKDIR:?GSE303006_WORKDIR is not set}"
RSCRIPT="${TDCHARM_RSCRIPT:-Rscript}"

if [ -n "${TDCHARM_R_SETUP:-}" ]; then
    # conda/module activation scripts reference unset variables and may
    # return non-zero internally: run them without -e/-u.
    set +eu
    eval "${TDCHARM_R_SETUP}"
    rc=$?
    set -eu
    if [ "$rc" -ne 0 ]; then
        echo "TDCHARM_R_SETUP failed (status $rc): ${TDCHARM_R_SETUP}" >&2
        exit 1
    fi
fi
if ! command -v "${TDCHARM_RSCRIPT:-Rscript}" >/dev/null 2>&1; then
    echo "Rscript not found after TDCHARM_R_SETUP: ${TDCHARM_R_SETUP:-<none>}" >&2
    exit 1
fi

# Temporary directory for R (and other tools). The default TMPDIR of a
# batch job may not exist or not be writable on the compute node, which
# makes R stop with "cannot create 'R_TempDir'". Use a directory inside
# the working directory instead, one per job.
TDCHARM_TMP="${GSE303006_WORKDIR}/tmp/${TDCHARM_STEP}_${PBS_JOBID:-${JOB_ID:-$$}}"
mkdir -p "${TDCHARM_TMP}"
export TMPDIR="${TDCHARM_TMP}" TMP="${TDCHARM_TMP}" TEMP="${TDCHARM_TMP}"
trap 'rc=$?; rm -rf "${TDCHARM_TMP}" 2>/dev/null; if [ "$rc" -ne 0 ] && [ -n "${SGE_ROOT:-}" ] && [ -n "${JOB_ID:-}" ]; then exit 100; fi' EXIT

# Keep BLAS threads equal to the allocated cores
NCORES="${NSLOTS:-${NCPUS:-${PBS_NUM_PPN:-1}}}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$NCORES}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-$NCORES}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-$NCORES}"

echo "============================================================"
echo "TDCHARM step : ${TDCHARM_STEP}"
echo "Repository   : ${TDCHARM_REPO}"
echo "Workdir      : ${GSE303006_WORKDIR}"
echo "Host         : $(hostname)"
echo "Job ID       : ${JOB_ID:-${PBS_JOBID:-none}}"
echo "Cores        : ${NCORES}"
echo "TMPDIR       : ${TMPDIR}"
echo "Start        : $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "Rscript      : $(command -v "${RSCRIPT}")"
"${RSCRIPT}" --version 2>&1 | head -1
echo "Git state    :"
git -C "${TDCHARM_REPO}" rev-parse --short HEAD 2>/dev/null || true
git -C "${TDCHARM_REPO}" status --short 2>/dev/null | head -50 || true
echo "============================================================"

export TDCHARM_REPO
export GSE303006_WORKDIR

case "${TDCHARM_STEP}" in
    synthetic)
        # Smoke test on synthetic data: original scripts 05-11, 17-18
        # (12-16 need GO annotation and are skipped) and the revision
        # scripts 19-27.
        "${RSCRIPT}" "${TDCHARM_REPO}/tests/make_synthetic_workdir.R" "${GSE303006_WORKDIR}"
        cd "${GSE303006_WORKDIR}"
        for s in 05 06 07 08 09 10 11 17 18; do
            f=$(ls "${TDCHARM_REPO}"/R/${s}_*.R)
            echo "--- ${f}"
            "${RSCRIPT}" -e "source('${f}')" > "log_synthetic_${s}.txt" 2>&1 || {
                echo "FAILED: ${f}"; tail -30 "log_synthetic_${s}.txt"; exit 1; }
        done
        # Synthetic gene names (Gene001, ...) have no GO annotation, so GSEA
        # is skipped here; it runs in the real "analysis"/"revision" steps.
        TDCHARM_RUN_GSEA=0 TDCHARM_N_ROT="${TDCHARM_N_ROT:-2}" TDCHARM_N_PERM="${TDCHARM_N_PERM:-10}" \
            "${RSCRIPT}" "${TDCHARM_REPO}/run_revision.R"
        ;;
    preprocessing)
        cd "${GSE303006_WORKDIR}"
        "${RSCRIPT}" "${TDCHARM_REPO}/run_preprocessing.R"
        ;;
    analysis)
        cd "${GSE303006_WORKDIR}"
        "${RSCRIPT}" "${TDCHARM_REPO}/run_analysis.R"
        ;;
    revision)
        cd "${GSE303006_WORKDIR}"
        "${RSCRIPT}" "${TDCHARM_REPO}/run_revision.R"
        ;;
    script)
        # One script only, e.g. TDCHARM_SCRIPT=23b_CHARM_link_concordance.R
        : "${TDCHARM_SCRIPT:?TDCHARM_SCRIPT is not set}"
        cd "${GSE303006_WORKDIR}"
        "${RSCRIPT}" "${TDCHARM_REPO}/R/${TDCHARM_SCRIPT}"
        ;;
    *)
        echo "Unknown TDCHARM_STEP: ${TDCHARM_STEP}" >&2
        exit 2
        ;;
esac

echo "============================================================"
echo "Finished ${TDCHARM_STEP}: $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "============================================================"
