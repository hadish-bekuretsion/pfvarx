# Stage 4 — unified QC, dashboard data and report export.

rule unified_qc:
    input:
        evidence=expand(f"{RESULTS}/{{sample}}/{{sample}}.evidence.tsv", sample=SAMPLES),
    output:
        run_qc=f"{RESULTS}/run_qc.json",
    params:
        cli=pfvarx_cli(),
        cfg=CONFIG_FILE,
    log:
        f"{LOGS}/qc.log",
    shell:
        """
        {params.cli} qc \
            --config {params.cfg} \
            --results {RESULTS} \
            --out {output.run_qc} 2> {log}
        """


rule export_report:
    input:
        manifest=f"{RESULTS}/run_manifest.json",
        run_qc=f"{RESULTS}/run_qc.json",
        evidence=expand(f"{RESULTS}/{{sample}}/{{sample}}.evidence.tsv", sample=SAMPLES),
    output:
        pdf=f"{RESULTS}/reports/pfvarx_results.pdf",
        xlsx=f"{RESULTS}/reports/pfvarx_results.xlsx",
    params:
        cli=pfvarx_cli(),
        cfg=CONFIG_FILE,
    log:
        f"{LOGS}/report.log",
    shell:
        """
        {params.cli} report all \
            --config {params.cfg} \
            --results {RESULTS} 2> {log}
        """
