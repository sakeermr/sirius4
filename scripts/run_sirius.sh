#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Runs the SIRIUS CLI inside the rformassspectrometry/rusirius container on a
# folder of .ms files produced by scripts/convert_to_ms.py
#
#   ./scripts/run_sirius.sh work/ms_files work/project
#
# Environment variables
#   IMAGE            docker image                (default rformassspectrometry/rusirius)
#   SIRIUS_TOOLS     tool chain                  (default "formula fingerprint structure canopus")
#   SIRIUS_ARGS      extra global sirius flags   (default "--maxmz 800")
#   PPM_MAX          MS1 mass accuracy in ppm    (default 10)
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

# --- locate the sirius executable inside the image ------------------------
echo "==> locating sirius binary"
SIRIUS_BIN="${SIRIUS_BIN:-$(docker run --rm --entrypoint /bin/sh "$IMAGE" -c '
  command -v sirius 2>/dev/null && exit 0
  for d in /opt /usr/local /usr/lib /root /home /srv; do
    p=$(find "$d" -maxdepth 6 -type f -name sirius -perm -u+x 2>/dev/null | head -n1)
    [ -n "$p" ] && echo "$p" && exit 0
  done
  echo ""' | tr -d "\r" | tail -n1)}"

if [ -z "$SIRIUS_BIN" ]; then
  echo "::error::could not find a sirius executable inside $IMAGE."
  echo "Inspect the image and set SIRIUS_BIN, e.g.:"
  echo "  docker run --rm -it --entrypoint /bin/bash $IMAGE"
  exit 2
fi
echo "    sirius binary: $SIRIUS_BIN"

DOCKER_ENV=()
[ -n "${SIRIUS_USER:-}" ]     && DOCKER_ENV+=(-e "SIRIUS_USER=${SIRIUS_USER}")
[ -n "${SIRIUS_PASSWORD:-}" ] && DOCKER_ENV+=(-e "SIRIUS_PASSWORD=${SIRIUS_PASSWORD}")

echo "==> running SIRIUS on $(ls "$MS_DIR"/*.ms | wc -l) .ms file(s)"
docker run --rm \
  -v "$PWD:/work" -w /work \
  "${DOCKER_ENV[@]}" \
  --entrypoint /bin/bash "$IMAGE" -lc "
    set -euo pipefail
    SIR='$SIRIUS_BIN'
    \$SIR --version || true

    # CSI:FingerID / CANOPUS need an account (SIRIUS >= 5.8)
    if [ -n \"\${SIRIUS_USER:-}\" ] && [ -n \"\${SIRIUS_PASSWORD:-}\" ]; then
      echo '==> sirius login'
      \$SIR login --user-env=SIRIUS_USER --password-env=SIRIUS_PASSWORD || \
      \$SIR login -u \"\$SIRIUS_USER\" -p \"\$SIRIUS_PASSWORD\" || \
        echo '::warning::sirius login failed - continuing (formula step still works)'
    else
      echo '::warning::no SIRIUS_USER/SIRIUS_PASSWORD secret - webservice tools may fail'
    fi

    \$SIR $SIRIUS_ARGS \
      --input '$MS_DIR' \
      --project '$PROJECT' \
      config --IsotopeSettings.filter=true \
             --FormulaSearchDB=none \
             --StructureSearchDB=BIO \
             --MS1MassDeviation.allowedMassDeviation=${PPM_MAX}ppm \
             --MS2MassDeviation.allowedMassDeviation=${PPM_MAX}ppm \
             --AdductSettings.enforced='$ADDUCTS' \
             --AlgorithmProfile=$INSTRUMENT \
      $SIRIUS_TOOLS \
      write-summaries --output '$PROJECT/summaries'
  "

echo "==> SIRIUS finished. Project: $PROJECT"
ls -R "$PROJECT" | head -n 50
