#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 SYSTEM_DIR [REPLICATE ...]" >&2
    echo "Example: $0 /path/to/PGA34A 1 2 3 4 5" >&2
    exit 2
}

[[ $# -ge 1 ]] || usage

system_dir=$(cd "$1" && pwd)
shift
replicates=("$@")
if [[ ${#replicates[@]} -eq 0 ]]; then
    replicates=(1 2 3 4 5)
fi

input_dir="${system_dir}/input"
gmx_bin=${GMX_BIN:-gmx}
maxwarn=${GMX_MAXWARN:-5}
mdrun_args=()
if [[ -n ${GMX_MDRUN_ARGS:-} ]]; then
    read -r -a mdrun_args <<< "${GMX_MDRUN_ARGS}"
fi

command -v "${gmx_bin}" >/dev/null 2>&1 || {
    echo "ERROR: ${gmx_bin} is not available in PATH" >&2
    exit 1
}

required=(
    step3_input.gro
    step4.0_minimization.mdp
    step4.1_equilibration.mdp
    step5_production.mdp
    topol.top
    index.ndx
)
for file in "${required[@]}"; do
    [[ -f "${input_dir}/${file}" ]] || {
        echo "ERROR: missing ${input_dir}/${file}" >&2
        exit 1
    }
done
[[ -d "${input_dir}/toppar" ]] || {
    echo "ERROR: missing ${input_dir}/toppar" >&2
    exit 1
}

run_mdrun() {
    local name=$1
    shift
    "${gmx_bin}" mdrun -v -deffnm "${name}" -ntmpi 1 \
        "$@" "${mdrun_args[@]}"
}

run_replica() {
    local replica=$1
    local output_dir="${system_dir}/rep${replica}"
    local mini=step4.0_minimization
    local equi=step4.1_equilibration
    local prod=step5_1

    mkdir -p "${output_dir}"
    cd "${output_dir}"

    if [[ ! -f "${mini}.gro" ]]; then
        echo "[rep${replica}] energy minimization"
        "${gmx_bin}" grompp \
            -f "${input_dir}/${mini}.mdp" \
            -o "${mini}.tpr" \
            -c "${input_dir}/step3_input.gro" \
            -r "${input_dir}/step3_input.gro" \
            -p "${input_dir}/topol.top" \
            -n "${input_dir}/index.ndx" \
            -maxwarn "${maxwarn}"
        run_mdrun "${mini}"
    fi

    if [[ ! -f "${equi}.gro" ]]; then
        echo "[rep${replica}] 125 ps restrained equilibration"
        "${gmx_bin}" grompp \
            -f "${input_dir}/${equi}.mdp" \
            -o "${equi}.tpr" \
            -c "${mini}.gro" \
            -r "${mini}.gro" \
            -p "${input_dir}/topol.top" \
            -n "${input_dir}/index.ndx" \
            -maxwarn "${maxwarn}"
        run_mdrun "${equi}" -cpi "${equi}.cpt"
    fi

    if [[ ! -f "${prod}.gro" ]]; then
        echo "[rep${replica}] 100 ns production"
        "${gmx_bin}" grompp \
            -f "${input_dir}/step5_production.mdp" \
            -o "${prod}.tpr" \
            -c "${equi}.gro" \
            -r "${equi}.gro" \
            -t "${equi}.cpt" \
            -p "${input_dir}/topol.top" \
            -n "${input_dir}/index.ndx" \
            -maxwarn "${maxwarn}"
        run_mdrun "${prod}" -cpi "${prod}.cpt"
    fi

    if [[ ! -f md_fit_20_100ns.xtc ]]; then
        echo "[rep${replica}] prepare aligned 20–100 ns trajectory"
        printf '%s\n%s\n' "${CENTER_GROUP:-1}" "${OUTPUT_GROUP:-0}" | \
            "${gmx_bin}" trjconv -s "${prod}.tpr" -f "${prod}.xtc" \
                -o prod_100ns_noPBC.xtc -pbc mol -center
        printf '%s\n%s\n' "${FIT_GROUP:-4}" "${OUTPUT_GROUP:-0}" | \
            "${gmx_bin}" trjconv -s "${prod}.tpr" -f prod_100ns_noPBC.xtc \
                -o prod_100ns_fit.xtc -fit rot+trans
        printf '%s\n' "${OUTPUT_GROUP:-0}" | \
            "${gmx_bin}" trjconv -s "${prod}.tpr" -f prod_100ns_fit.xtc \
                -o md_fit_20_100ns.xtc -b 20000 -e 100000
    fi

    [[ -f step5_ref.tpr ]] || cp -p "${prod}.tpr" step5_ref.tpr

    "${gmx_bin}" check -f md_fit_20_100ns.xtc > check_20_100ns.txt 2>&1
    echo "[rep${replica}] complete"
}

for replica in "${replicates[@]}"; do
    [[ ${replica} =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: invalid replicate number: ${replica}" >&2
        exit 1
    }
    run_replica "${replica}"
done
