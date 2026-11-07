#!/bin/bash
# ============================================================
# jobs/submit_all.sh
#
# Submits the TDCHARM workflow as a chain of qsub jobs. Each job
# starts only if the previous one finished successfully.
#
#   synthetic  ->  preprocessing  ->  analysis  ->  revision
#   (smoke test     (01-04)            (05-18)      (19-27)
#    on synthetic
#    data)
#
# Supported schedulers: SGE / Univa-Altair Grid Engine ("sge"),
# PBS Pro ("pbspro") and Torque ("torque"). The scheduler is
# detected automatically or set with SCHEDULER=.
#
# Usage (from anywhere):
#   GSE303006_WORKDIR=/path/to/workdir bash jobs/submit_all.sh
#
# Settings (environment variables; defaults in brackets):
#   GSE303006_WORKDIR  working directory with the GEO input files (required)
#   SYNTH_DIR          directory for the synthetic smoke test
#                      [$GSE303006_WORKDIR/../tdcharm_synthetic_test]
#   STEPS              steps to submit, in order
#                      ["synthetic preprocessing analysis revision"]
#   SCHEDULER          sge | pbspro | torque [auto]
#   SITE               site profile, e.g. kyoto_icr (jobs/site_kyoto_icr.sh)
#   QUEUE              queue name [none]; QUEUE_<STEP> per step
#   NCPU               cores per job [4]
#   MEM                memory per job, e.g. 64G [64G]
#                      (SGE: divided by NCPU and given as h_vmem per slot)
#   WALLTIME           wall time per job, hh:mm:ss [72:00:00]
#   MEM_<STEP>, WALLTIME_<STEP>, NCPU_<STEP>
#                      per-step overrides, e.g. MEM_revision=96G
#   SGE_PE             SGE parallel environment [smp]
#   SGE_MEM_RESOURCE   SGE memory resource name [h_vmem]
#   EXTRA_QSUB_OPTS    extra options added to every qsub call
#   TDCHARM_R_SETUP    command to make R available inside the job,
#                      e.g. "module load R/4.3.3"
#   TDCHARM_RSCRIPT    Rscript executable inside the job [Rscript]
#   DRY_RUN=1          print the qsub commands without submitting
#
# Any other TDCHARM_* variable (e.g. TDCHARM_RUN_GSEA,
# TDCHARM_EXTERNAL_LINKS) set when calling this script is passed on
# to the jobs.
# ============================================================

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STEP_SCRIPT="${REPO}/jobs/tdcharm_step.sh"

# Site profile (e.g. SITE=kyoto_icr loads jobs/site_kyoto_icr.sh).
# Values already set in the environment take precedence.
if [ -n "${SITE:-}" ]; then
    SITE_FILE="${REPO}/jobs/site_${SITE}.sh"
    [ -f "$SITE_FILE" ] || { echo "Unknown SITE: ${SITE} (${SITE_FILE} not found)" >&2; exit 1; }
    # shellcheck disable=SC1090
    source "$SITE_FILE"
fi

: "${GSE303006_WORKDIR:?Set GSE303006_WORKDIR to the analysis working directory}"
WORKDIR="$(cd "${GSE303006_WORKDIR}" && pwd)"
SYNTH_DIR="${SYNTH_DIR:-$(dirname "${WORKDIR}")/tdcharm_synthetic_test}"
STEPS="${STEPS:-synthetic preprocessing analysis revision}"
NCPU="${NCPU:-4}"
MEM="${MEM:-64G}"
WALLTIME="${WALLTIME:-72:00:00}"
SGE_PE="${SGE_PE:-smp}"
SGE_MEM_RESOURCE="${SGE_MEM_RESOURCE:-h_vmem}"
EXTRA_QSUB_OPTS="${EXTRA_QSUB_OPTS:-}"
DRY_RUN="${DRY_RUN:-0}"
LOGDIR="${WORKDIR}/qsub_logs"
mkdir -p "${LOGDIR}"

# ---------- scheduler detection ----------
if [ -z "${SCHEDULER:-}" ]; then
    if [ -n "${SGE_ROOT:-}" ] || qstat -help 2>&1 | grep -qi "grid engine\|GE [0-9]"; then
        SCHEDULER=sge
    elif qstat --version 2>&1 | grep -qi "pbs"; then
        SCHEDULER=pbspro
    elif qstat --version 2>&1 | grep -qi "torque\|version: [0-9]"; then
        SCHEDULER=torque
    else
        echo "Cannot detect the scheduler; set SCHEDULER=sge|pbspro|torque" >&2
        exit 1
    fi
fi
case "${SCHEDULER}" in sge|pbspro|torque) ;; *)
    echo "Unknown SCHEDULER: ${SCHEDULER}" >&2; exit 1 ;; esac

# ---------- helpers ----------
mem_gb() {   # "64G" / "64gb" / "64000M" -> integer GB
    local m="${1^^}"; m="${m%B}"
    case "$m" in
        *T) echo $(( ${m%T} * 1024 )) ;;
        *G) echo "${m%G}" ;;
        *M) echo $(( (${m%M} + 1023) / 1024 )) ;;
        *)  echo "$m" ;;
    esac
}

pass_vars() {  # comma-separated VAR=value list for qsub -v
    local step="$1" wd="$2" out="TDCHARM_STEP=${step},TDCHARM_REPO=${REPO},GSE303006_WORKDIR=${wd}"
    local v
    for v in $(compgen -e | grep '^TDCHARM_' || true); do
        case "$v" in TDCHARM_STEP|TDCHARM_REPO) continue ;; esac
        local val="${!v}"
        if [[ "$val" == *","* ]]; then
            echo "WARNING: ${v} contains a comma and cannot be passed with qsub -v;" \
                 "set it inside TDCHARM_R_SETUP instead." >&2
            continue
        fi
        out="${out},${v}=${val}"
    done
    echo "$out"
}

submit() {  # submit <step> <dependency-id or empty>; prints job id
    local step="$1" dep="$2"
    local wd="${WORKDIR}"
    [ "$step" = "synthetic" ] && wd="${SYNTH_DIR}"
    local ncpu_var="NCPU_${step}" mem_var="MEM_${step}" wt_var="WALLTIME_${step}"
    local q_var="QUEUE_${step}"
    local ncpu="${!ncpu_var:-$NCPU}" mem="${!mem_var:-$MEM}" wt="${!wt_var:-$WALLTIME}"
    local queue="${!q_var:-${QUEUE:-}}"
    [ "$step" = "synthetic" ] && { ncpu="${NCPU_synthetic:-1}"; mem="${MEM_synthetic:-8G}"; wt="${WALLTIME_synthetic:-02:00:00}"; }
    local name="tdcharm_${step}"
    local vars; vars="$(pass_vars "$step" "$wd")"
    local log="${LOGDIR}/${name}.log"
    local -a cmd=(qsub)

    case "${SCHEDULER}" in
        sge)
            local per_slot=$(( ( $(mem_gb "$mem") + ncpu - 1 ) / ncpu ))
            cmd+=(-terse -N "$name" -S /bin/bash -cwd -j y -o "$log"
                  -pe "$SGE_PE" "$ncpu"
                  -l "${SGE_MEM_RESOURCE}=${per_slot}G,h_rt=${wt}"
                  -v "$vars")
            [ -n "$queue" ] && cmd+=(-q "$queue")
            [ -n "$dep" ] && cmd+=(-hold_jid "$dep")
            ;;
        pbspro)
            cmd+=(-N "$name" -S /bin/bash -j oe -o "$log"
                  -l "select=1:ncpus=${ncpu}:mem=$(mem_gb "$mem")gb"
                  -l "walltime=${wt}" -v "$vars")
            [ -n "$queue" ] && cmd+=(-q "$queue")
            [ -n "$dep" ] && cmd+=(-W "depend=afterok:${dep}")
            ;;
        torque)
            cmd+=(-N "$name" -S /bin/bash -j oe -o "$log"
                  -l "nodes=1:ppn=${ncpu},mem=$(mem_gb "$mem")gb,walltime=${wt}"
                  -v "$vars")
            [ -n "$queue" ] && cmd+=(-q "$queue")
            [ -n "$dep" ] && cmd+=(-W "depend=afterok:${dep}")
            ;;
    esac
    # shellcheck disable=SC2206
    [ -n "$EXTRA_QSUB_OPTS" ] && cmd+=($EXTRA_QSUB_OPTS)
    cmd+=("$STEP_SCRIPT")

    if [ "$DRY_RUN" = "1" ]; then
        printf '%q ' "${cmd[@]}" >&2; echo >&2
        echo "DRYRUN_${step}"
    else
        local id
        id="$(cd "$WORKDIR" && "${cmd[@]}")"
        # SGE -terse prints "123" (or "123.1-10:1" for arrays): keep the number.
        # PBS/Torque print "123.server": keep the full id for -W depend.
        [ "${SCHEDULER}" = "sge" ] && id="${id%%.*}"
        echo "$id"
    fi
}

echo "Site      : ${SITE:-none}"
echo "Scheduler : ${SCHEDULER}"
echo "Repository: ${REPO}"
echo "Workdir   : ${WORKDIR}"
echo "Synthetic : ${SYNTH_DIR}"
echo "Logs      : ${LOGDIR}"
echo "Steps     : ${STEPS}"
echo

prev=""
for step in ${STEPS}; do
    id="$(submit "$step" "$prev")"
    echo "submitted ${step}: job ${id}${prev:+ (after ${prev})}"
    prev="$id"
done
