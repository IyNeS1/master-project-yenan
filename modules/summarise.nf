process SUMMARISE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path "benchmark_summary.csv"

    script:
    """
    first=\$(ls *.csv | head -1)
    head -1 \$first | sed 's/^/model,kvalue,/' > benchmark_summary.csv
    for f in *.csv; do
        base=\$(basename \$f .csv)
        model=\$(echo \$base | sed 's/_k[0-9]*//')
        kval=\$(echo \$base | grep -o 'k[0-9]*' | sed 's/k//')
        tail -n +2 \$f | while IFS= read -r line; do
            echo "\${model},\${kval},\${line}"
        done
    done >> benchmark_summary.csv
    """
}
