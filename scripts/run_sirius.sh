#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Runs the SIRIUS CLI inside the rformassspectrometry/rusirius container on a
# folder of .ms files produced by scripts/convert_to_ms.py
#
#   ./scripts/run_sirius.sh work/ms_files work/project
#
# Subcommand names differ between SIRIUS 5 and 6 (e.g. 'structure' does not
# exist in 6.3), so the requested tool chain is validated against the image's
# own --help output before running, and unknown names are remapped or dropped.
#
# Environment variables
#   IMAGE            docker image                (default rformassspectrometry/rusirius)
#   SIRIUS_TOOLS     tool chain                  (default "formula fingerprint structure canopus")
#   SIRIUS_ARGS      extra global sirius flags   (default "--maxmz 800")
#   PPM_MAX          mass accuracy in ppm        (default 10)
#   INSTRUMENT       orbitrap | qtof             (default orbitrap)
#   ADDUCTS          e.g. "[M+H]+,[M+Na]+"       (default [M+H]+)
#   SIRIUS_USER /
#   SIRIUS_PASSWORD  CSI:FingerID login (needed for fingerprint/structure/canopus)
#   SIRIUS_BIN       override path to the sirius executable inside the image
# ---------------------------------------------------------------------------
set -euo pipefail

MS_DIR="${1:-work/ms_files}"
PROJECT="${2:-work/project}"

IMAGE="${IMAGE:-rformassspectrometry/rusirius}"
SIRIUS_TOOLS="${SIRIUS_TOOLS:-formula fingerprint structure canopus}"
SIRIUS_ARGS="${SIRIUS_ARGS:---maxmz 800}"
PPM_MAX="${PPM_MAX:-10}"
INSTRUMENT="${INSTRUMENT:-orbitrap}"
ADDUCTS="${ADDUCTS:-[M+H]+}"

if [ ! -d "$MS_DIR" ] || [ -z "$(ls -A "$MS_DIR"/*.ms 2>/dev/null)" ]; then
  echo "::error::no .ms files in $MS_DIR - run scripts/convert_to_ms.py first"
  exit 1
fi
mkdir -p "$PROJECT"

echo "==> pulling $IMAGE"
docker pull "$IMAGE"

echo "==> locating sirius binary"
SIRIUS_BIN="${SIRIUS_BIN:-$(docker run --rm --entrypoint /bin/sh "$IMAGE" -c '
  command -v sirius 2>/dev/null && exit 0
  for d in /opt /usr/local /usr/lib /root /home /srv; do
    p=$(find "$d" -maxdepth 6 -type f -name sirius -perm -u+x 2>/dev/null | head -n1)
    [ -n "$p" ] && echo "$p" && exit 0
  done
  echo ""' | tr -d "\r" | tail -n1)}"

if [ -z "$SIRIUS_BIN" ]; then
  echo "::error::could not find a sirius executable inside $IMAGE. Set SIRIUS_BIN."
  exit 2
fi
echo "    sirius binary: $SIRIUS_BIN"

DOCKER_ENV=()
[ -n "${SIRIUS_USER:-}" ]     && DOCKER_ENV+=(-e "SIRIUS_USER=${SIRIUS_USER}")
[ -n "${SIRIUS_PASSWORD:-}" ] && DOCKER_ENV+=(-e "SIRIUS_PASSWORD=${SIRIUS_PASSWORD}")

echo "==> running SIRIUS on $(ls "$MS_DIR"/*.ms | wc -l) .ms file(s)"

# the inner script is written to a file to keep the quoting sane
cat > /tmp/_sirius_inner.sh <<'INNER'
set -uo pipefail
SIR="$SIRIUS_BIN_IN"

"$SIR" --version 2>&1 | grep -iE "^SIRIUS|lib:" || true

echo "==> available subcommands in this SIRIUS build"
HELP="$("$SIR" --help 2>&1 || true)"
echo "$HELP" | sed -n '/[Cc]ommands:/,$p' | head -n 45

have() { echo "$HELP" | grep -qE "^[[:space:]]*$1([[:space:],]|$)"; }

# map each requested tool onto a name this build actually knows
resolve() {
  case "$1" in
    formula)      cands="formula formulas sirius" ;;
    zodiac)       cands="zodiac" ;;
    fingerprint)  cands="fingerprint fingerprints fingerprint-search" ;;
    structure)    cands="structure structures structure-db-search structuredb" ;;
    canopus)      cands="canopus compound-classes compound-class" ;;
    *)            cands="$1" ;;
  esac
  for c in $cands; do
    if have "$c"; then echo "$c"; return 0; fi
  done
  return 1
}

TOOLS=""
for t in $SIRIUS_TOOLS_IN; do
  if pick="$(resolve "$t")"; then
    TOOLS="$TOOLS $pick"
    [ "$pick" != "$t" ] && echo "::warning::'$t' is not a subcommand here - using '$pick' instead"
  else
    echo "::warning::subcommand '$t' does not exist in this SIRIUS build - skipping it"
  fi
done
TOOLS="$(echo "$TOOLS" | xargs || true)"
[ -z "$TOOLS" ] && TOOLS="formula"
echo "==> resolved tool chain: $TOOLS"

SUMMARY=""
if have "write-summaries"; then
  SUMMARY="write-summaries"
elif have "write-summary"; then
  SUMMARY="write-summary"
else
  echo "::warning::no write-summaries subcommand - the project will still contain results"
fi

if [ -n "${SIRIUS_USER:-}" ] && [ -n "${SIRIUS_PASSWORD:-}" ]; then
  echo "==> sirius login"
  "$SIR" login --user-env=SIRIUS_USER --password-env=SIRIUS_PASSWORD \
    || "$SIR" login -u "$SIRIUS_USER" -p "$SIRIUS_PASSWORD" \
    || echo "::warning::sirius login failed - only offline tools will work"
else
  echo "::warning::SIRIUS_USER and/or SIRIUS_PASSWORD not set - fingerprint/structure/canopus need BOTH"
fi

set -x
"$SIR" $SIRIUS_ARGS_IN \
  --input "$MS_DIR_IN" \
  --project "$PROJECT_IN" \
  config --IsotopeSettings.filter=true \
         --FormulaSearchDB=none \
         --StructureSearchDB=BIO \
         --MS1MassDeviation.allowedMassDeviation=${PPM_MAX_IN}ppm \
         --MS2MassDeviation.allowedMassDeviation=${PPM_MAX_IN}ppm \
         --AdductSettings.enforced="$ADDUCTS_IN" \
         --AlgorithmProfile=$INSTRUMENT_IN \
  $TOOLS $SUMMARY
RC=$?
set +x
echo "==> sirius exit code: $RC"
exit $RC
INNER

docker run --rm \
  -v "$PWD:/work" -w /work \
  -v /tmp/_sirius_inner.sh:/tmp/inner.sh:ro \
  "${DOCKER_ENV[@]}" \
  -e "SIRIUS_BIN_IN=$SIRIUS_BIN" \
  -e "SIRIUS_TOOLS_IN=$SIRIUS_TOOLS" \
  -e "SIRIUS_ARGS_IN=$SIRIUS_ARGS" \
  -e "MS_DIR_IN=$MS_DIR" \
  -e "PROJECT_IN=$PROJECT" \
  -e "PPM_MAX_IN=$PPM_MAX" \
  -e "INSTRUMENT_IN=$INSTRUMENT" \
  -e "ADDUCTS_IN=$ADDUCTS" \
  --entrypoint /bin/bash "$IMAGE" /tmp/inner.sh

echo "==> verifying SIRIUS actually wrote a project"
if [ -z "$(ls -A "$PROJECT" 2>/dev/null)" ]; then
  echo "::error::SIRIUS finished but '$PROJECT' is empty - no results were produced."
  echo "Check the 'available subcommands' list above and adjust SIRIUS_TOOLS."
  exit 3
fi

N_SUM=$(find "$PROJECT" -name "*.tsv" 2>/dev/null | wc -l)
echo "==> SIRIUS finished. Project: $PROJECT ($N_SUM summary tables)"
ls -R "$PROJECT" | head -n 40
