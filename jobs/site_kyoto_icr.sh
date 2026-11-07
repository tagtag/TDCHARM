# ============================================================
# jobs/site_kyoto_icr.sh
#
# Defaults for the Kyoto University ICR supercomputer system
# (化研スパコン, PBS Professional). Loaded by submit_all.sh when
# SITE=kyoto_icr. Every value can still be overridden by setting the
# variable before calling submit_all.sh.
#
# Based on the ICR user documentation (qsub options, job script
# examples, queue information):
#   https://www.scl.kyoto-u.ac.jp/Attention/batch_op.html
#   https://www.scl.kyoto-u.ac.jp/Attention/batch_script_sample.html
#   https://www.scl.kyoto-u.ac.jp/Attention/batch_stat.html
# Queue names and limits can change; check them with
#   qstat -Q ; qstat -Qf APC ; qstatmyjobs
# ============================================================

: "${SCHEDULER:=pbspro}"

# Short smoke test on synthetic data: SMALL queue
# (documented limits: 12 cores, 48 GB, 12 h).
: "${QUEUE_synthetic:=SMALL}"
: "${NCPU_synthetic:=1}"
: "${MEM_synthetic:=8G}"
: "${WALLTIME_synthetic:=02:00:00}"

# Real data: APC queue (nodes with 56 cores and 240/490/1000 GB;
# documented per-job limit 56 cores, 980 GB; no walltime maximum).
: "${QUEUE:=APC}"
: "${NCPU:=8}"
: "${MEM:=128G}"
: "${WALLTIME:=336:00:00}"

# R inside batch jobs: R is provided by the conda environment r_env.
# Batch shells do not read ~/.bashrc, so conda is initialized from its
# installation directory, which is looked up here on the login node.
: "${TDCHARM_CONDA_ENV:=r_env}"
if [ -z "${TDCHARM_R_SETUP:-}" ]; then
    if [ -z "${CONDA_BASE:-}" ]; then
        if command -v conda >/dev/null 2>&1; then
            CONDA_BASE="$(conda info --base 2>/dev/null)"
        elif [ -n "${CONDA_EXE:-}" ]; then
            CONDA_BASE="$(dirname "$(dirname "${CONDA_EXE}")")"
        fi
    fi
    if [ -z "${CONDA_BASE:-}" ] || [ ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
        echo "Cannot locate conda; run 'conda activate ${TDCHARM_CONDA_ENV}' once on the login node," \
             "or set CONDA_BASE=/path/to/miniconda3 (or TDCHARM_R_SETUP)." >&2
        exit 1
    fi
    TDCHARM_R_SETUP="source ${CONDA_BASE}/etc/profile.d/conda.sh; conda activate ${TDCHARM_CONDA_ENV}"
fi
export TDCHARM_R_SETUP

# MACS3 (needed only when the peak files are regenerated, script 01).
# It lives in its own conda environment because it needs a different
# Python than the R environment:
#   conda create -n macs3 -c conda-forge -c bioconda python=3.12 macs3=3.1.0
: "${TDCHARM_MACS3_ENV:=macs3}"
if [ -z "${TDCHARM_MACS3:-}" ] && [ -n "${CONDA_BASE:-}" ] && \
   [ -x "${CONDA_BASE}/envs/${TDCHARM_MACS3_ENV}/bin/macs3" ]; then
    TDCHARM_MACS3="${CONDA_BASE}/envs/${TDCHARM_MACS3_ENV}/bin/macs3"
fi
if [ -n "${TDCHARM_MACS3:-}" ]; then export TDCHARM_MACS3; fi
