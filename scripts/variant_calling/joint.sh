# #!/bin/bash
# set -e

# module load singularity
# mkdir -p "${snakemake_params[gl_nexus_prefix]}"
# export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/CEPH-K1463-TandemRepeats/deep_variant_tmp/

# singularity exec --cleanenv -H $SINGULARITYENV_TMPDIR -B /usr/lib/locale/:/usr/lib/locale/ \
#     ${snakemake_input[sif]} \
#     /usr/local/bin/glnexus_cli \
#     --dir ${snakemake_params[gl_nexus_prefix]} \
#     --config DeepVariant_unfiltered \
#     --mem-gbytes 64 \
#     --threads ${snakemake[threads]} \
#     ${snakemake_input[gvcfs]} > ${snakemake_output}


#!/bin/bash
#!/bin/bash
set -euo pipefail

module load singularity

# keep your tmpdir env for the container user
export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/CEPH-K1463-TandemRepeats/deep_variant_tmp/

# Snakemake-provided variables
RAW_DB="${snakemake_params[gl_nexus_prefix]}"   # final desired DB dir (under out_dir)
OUT_VCF="${snakemake_output}"
SIF="${snakemake_input[sif]}"

# parent of final DB (ensure it exists)
DB_PARENT="$(dirname "$RAW_DB")"
mkdir -p "$DB_PARENT"

# create a unique name for the temporary DB in the same parent, but DO NOT create it.
# using timestamp+random reduces race chance and avoids mktemp -u issues
TMP_DB_NAME=".gl_nexus_db_tmp.$(date +%s%N)_$RANDOM"
TMP_DB="${DB_PARENT}/${TMP_DB_NAME}"

# ensure output dir exists
mkdir -p "$(dirname "$OUT_VCF")"

# build bind list:
BIND_OPTS=( )
# bind parent of DB (so container can create TMP_DB inside it)
BIND_OPTS+=( -B "${DB_PARENT}:${DB_PARENT}" )
# bind output dir
BIND_OPTS+=( -B "$(realpath "$(dirname "$OUT_VCF")")":"$(realpath "$(dirname "$OUT_VCF")")" )

# bind each gVCF parent dir (avoid duplicates)
declare -A seen
for g in ${snakemake_input[gvcfs]}; do
    d=$(realpath "$(dirname "$g")")
    if [[ -z "${seen[$d]:-}" ]]; then
        BIND_OPTS+=( -B "$d:$d" )
        seen[$d]=1
    fi
done

# also bind locale
BIND_OPTS+=( -B /usr/lib/locale/:/usr/lib/locale/ )

# run GLnexus writing into TMP_DB (which does NOT exist yet; GLnexus will create it)
set -x
singularity exec --cleanenv -H "$SINGULARITYENV_TMPDIR" "${BIND_OPTS[@]}" \
    "$SIF" \
    /usr/local/bin/glnexus_cli \
    --dir "$TMP_DB" \
    --config DeepVariant_unfiltered \
    --mem-gbytes 64 \
    --threads ${snakemake[threads]} \
    ${snakemake_input[gvcfs]} > "${OUT_VCF}" || {
        echo "GLnexus failed. If it created something, check: ${TMP_DB}"
        exit 1
    }

# At this point, GLnexus should have created "$TMP_DB" inside DB_PARENT.
if [[ ! -d "$TMP_DB" ]]; then
    echo "Expected GLnexus to create $TMP_DB but it is missing. Listing $DB_PARENT:"
    ls -la "$DB_PARENT"
    exit 1
fi

# If a DB already exists at the final path, move it away to a timestamped backup
if [[ -d "$RAW_DB" ]]; then
    BACKUP="${RAW_DB}.bak.$(date +%s)"
    echo "Existing DB at $RAW_DB -> moving to $BACKUP"
    mv "$RAW_DB" "$BACKUP"
fi

# move the freshly-created tmp DB into place (atomic when on same fs)
mv "$TMP_DB" "$RAW_DB"

# set permissions sensibly (do not fail if chown not permitted)
chown "$(whoami)" -R "$RAW_DB" 2>/dev/null || true
chmod -R u+rwX "$RAW_DB" 2>/dev/null || true

echo "GLnexus DB placed at $RAW_DB and VCF at $OUT_VCF"
