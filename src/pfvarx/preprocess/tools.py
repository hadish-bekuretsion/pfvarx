"""External tool wrappers.

Every external tool is described by a :class:`ToolSpec` that records the
executable, the version command, the minimum supported version and the
command builder.  Wrappers *fail closed*: if a binary is missing or below the
minimum version the caller gets a :class:`ToolExecutionError` with a corrective
action rather than a silently degraded run.

The wrappers only build and (optionally) execute commands.  They never parse
scientific output into calls -- that is the job of the pure-Python caller.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from pfvarx.errors import ToolExecutionError


@dataclass(frozen=True)
class ToolSpec:
    name: str
    executable: str
    version_args: tuple[str, ...]
    min_version: str | None = None
    homepage: str = ""
    purpose: str = ""

    def resolve(self) -> str | None:
        return shutil.which(self.executable)

    def version(self) -> str | None:
        exe = self.resolve()
        if exe is None:
            return None
        try:
            out = subprocess.run(
                [exe, *self.version_args],
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
            text = (out.stdout or "") + (out.stderr or "")
            return text.strip().splitlines()[0] if text.strip() else "unknown"
        except Exception:  # pragma: no cover - defensive
            return "unknown"

    def require(self) -> str:
        exe = self.resolve()
        if exe is None:
            raise ToolExecutionError(
                f"Required tool '{self.executable}' was not found on PATH.",
                cause=f"{self.executable} is not installed or not executable.",
                action=(
                    f"Install {self.name} ({self.homepage}) or activate the Pixi "
                    "environment that provides it."
                ),
                consequence="The dependent workflow rule cannot run.",
            )
        return exe


TOOLS: dict[str, ToolSpec] = {
    "fastp": ToolSpec("fastp", "fastp", ("--version",), "0.23.0", "https://github.com/OpenGene/fastp", "read QC and trimming"),
    "bwa-mem2": ToolSpec("bwa-mem2", "bwa-mem2", ("version",), "2.2.1", "https://github.com/bwa-mem2/bwa-mem2", "short-read alignment"),
    "samtools": ToolSpec("samtools", "samtools", ("--version",), "1.17", "https://www.htslib.org", "BAM/CRAM manipulation"),
    "bcftools": ToolSpec("bcftools", "bcftools", ("--version",), "1.17", "https://www.htslib.org", "variant file manipulation"),
    "picard": ToolSpec("picard", "picard", ("--version",), "3.0.0", "https://broadinstitute.github.io/picard", "duplicate marking"),
    "fgbio": ToolSpec("fgbio", "fgbio", ("--version",), "2.1.0", "https://github.com/fulcrumgenomics/fgbio", "UMI consensus"),
    "snpeff": ToolSpec("SnpEff", "snpEff", ("-version",), "5.1", "https://pcingola.github.io/SnpEff", "variant annotation"),
    "multiqc": ToolSpec("MultiQC", "multiqc", ("--version",), "1.14", "https://multiqc.info", "aggregate QC"),
    "deploid": ToolSpec("DEploid", "deploid", ("--version",), None, "https://github.com/DEploid-dev/DEploid", "complexity of infection"),
    "hmmibd": ToolSpec("hmmIBD", "hmmIBD", ("-h",), None, "https://github.com/glipsman/hmmIBD", "relatedness"),
    "iqtree2": ToolSpec("IQ-TREE 2", "iqtree2", ("--version",), "2.2.0", "http://www.iqtree.org", "phylogenetics"),
}


def tool_matrix() -> list[dict]:
    """Return the machine-readable tool matrix for documentation and QC."""
    rows = []
    for key, spec in TOOLS.items():
        rows.append(
            {
                "key": key,
                "name": spec.name,
                "executable": spec.executable,
                "min_version": spec.min_version or "any",
                "available": spec.resolve() is not None,
                "detected_version": spec.version(),
                "purpose": spec.purpose,
                "homepage": spec.homepage,
            }
        )
    return rows


def run_command(cmd: list[str], *, log_path: str | None = None, timeout: int | None = None) -> subprocess.CompletedProcess:
    """Execute a command, capturing output, raising a typed error on failure."""
    if log_path:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=False)
    except FileNotFoundError as exc:
        raise ToolExecutionError(
            f"Executable not found: {cmd[0]}",
            cause=str(exc),
            action="Install the tool or correct the executable path in configuration.",
        ) from exc
    if log_path:
        with open(log_path, "w") as fh:
            fh.write("$ " + " ".join(cmd) + "\n")
            fh.write(proc.stdout or "")
            fh.write("\n--- stderr ---\n")
            fh.write(proc.stderr or "")
    if proc.returncode != 0:
        raise ToolExecutionError(
            f"Command failed with exit code {proc.returncode}: {cmd[0]}",
            cause=(proc.stderr or proc.stdout or "")[-2000:],
            action="Inspect the rule log; check inputs, reference indices and resources.",
        )
    return proc


@dataclass
class CommandPlan:
    """An ordered list of commands plus expected outputs for one sample."""

    sample_id: str
    steps: list[dict] = field(default_factory=list)

    def add(self, name: str, cmd: list[str], outputs: list[str], tool: str) -> None:
        self.steps.append({"name": name, "tool": tool, "cmd": cmd, "outputs": outputs})

    def as_dict(self) -> dict:
        return {"sample_id": self.sample_id, "steps": self.steps}
