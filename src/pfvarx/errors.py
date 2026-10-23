"""Typed exception hierarchy.

Every error carries the stage, optional sample, cause, corrective action and
downstream consequence so that logs and the dashboard can render actionable
diagnostics (specification section 16).
"""

from __future__ import annotations


class PfVarXError(Exception):
    """Base class for all PfVarX errors."""

    stage: str = "unknown"

    def __init__(
        self,
        message: str,
        *,
        sample: str | None = None,
        cause: str | None = None,
        action: str | None = None,
        consequence: str | None = None,
    ) -> None:
        super().__init__(message)
        self.message = message
        self.sample = sample
        self.cause = cause or message
        self.action = action or "Inspect the stage log and configuration."
        self.consequence = consequence or "Stage output is not produced; downstream stages cannot run."

    def as_dict(self) -> dict:
        return {
            "error_type": type(self).__name__,
            "stage": self.stage,
            "sample": self.sample,
            "message": self.message,
            "cause": self.cause,
            "corrective_action": self.action,
            "downstream_consequence": self.consequence,
        }

    def __str__(self) -> str:  # pragma: no cover - formatting only
        bits = [f"[{self.stage}] {self.message}"]
        if self.sample:
            bits.append(f"sample={self.sample}")
        bits.append(f"cause={self.cause}")
        bits.append(f"action={self.action}")
        return " | ".join(bits)


class ConfigError(PfVarXError):
    stage = "configuration"


class InputValidationError(PfVarXError):
    stage = "input_validation"


class RoutingError(PfVarXError):
    stage = "routing"


class ReferenceMismatchError(InputValidationError):
    stage = "reference_validation"


class ToolExecutionError(PfVarXError):
    stage = "tool_execution"


class ModelError(PfVarXError):
    stage = "statistical_model"


class CalibrationError(ModelError):
    stage = "calibration"


class PolicyError(PfVarXError):
    stage = "decision_policy"


class UnsupportedModeError(PfVarXError):
    stage = "mode_support"


class OutputError(PfVarXError):
    stage = "output"
