# Stage 0 — input validation and routing.
#
# Emits a machine-readable routing plan and per-sample validation records.
# This rule is cheap and always runs first; it fails closed on any validation
# error so that no downstream stage consumes an invalid input.

rule validate_inputs:
    output:
        plan=f"{RESULTS}/routing_plan.json",
        sheet=f"{RESULTS}/validation/samplesheet.json",
    params:
        cli=pfvarx_cli(),
        cfg=CONFIG_FILE,
        samplesheet=config.get("inputs", {}).get("samplesheet"),
    log:
        f"{LOGS}/validate.log",
    resources:
        mem_mb=1024,
        runtime=10,
    shell:
        """
        mkdir -p {RESULTS}/validation {LOGS}
        {params.cli} validate \
            --config {params.cfg} \
            --samplesheet {params.samplesheet} \
            --out {output.plan} \
            --json > {output.sheet} 2> {log}
        """
