from dataclasses import dataclass

ROTOR_CUT_DEFAULT = 100.0  # wavenumbers
T_DEFAULT = 298.15  # Kelvin
p_DEFAULT = 101325  # Pascal


@dataclass(frozen=True)
class ThermoPolicy:
    """Immutable thermochemistry frequency-treatment policy (P05).

    Makes the two retained, deliberately different entry-point behaviours
    explicit instead of relying on scattered library defaults. Every value is
    passed to :func:`thermoanalysis.thermo.thermochemistry` explicitly, and the
    same policy is serialized (via :meth:`as_dict`) into the workflow JSON/YAML
    and onto the Geometry result so a payload always states which frequency
    treatment produced it.

    The field names are the stable, serialized public keys and are
    field-isomorphic between pdb2reaction and mlmm_toolkit;
    :meth:`thermochemistry_kwargs` maps them onto the ``thermochemistry``
    parameter names. Neither policy is a universal scientific default -- each is
    tied to its entry point.

    Fields
    ------
    kind
        Vibrational treatment, ``"qrrho"`` (quasi-RRHO) or ``"rrho"``.
    rotor_cutoff_cm
        Free-rotor interpolation cutoff in cm^-1 (QRRHO only).
    frequency_scale
        Multiplicative harmonic-frequency scale factor.
    zpe_scale
        Zero-point-energy scale factor, applied exactly once.
    invert_imag_from_cm
        Imaginary frequencies in ``[invert_imag_from_cm, 0]`` cm^-1 are inverted
        to real. ``0.0`` disables inversion (must be <= 0).
    positive_frequency_floor_cm
        Positive frequencies below this cm^-1 value are raised to it. ``0.0``
        disables the floor (must be >= 0).
    """

    kind: str = "qrrho"
    rotor_cutoff_cm: float = ROTOR_CUT_DEFAULT
    frequency_scale: float = 1.0
    zpe_scale: float = 1.0
    invert_imag_from_cm: float = 0.0
    positive_frequency_floor_cm: float = 0.0

    def thermochemistry_kwargs(self) -> dict:
        """Map the policy onto ``thermochemistry`` keyword arguments."""
        return {
            "kind": self.kind,
            "rotor_cutoff": self.rotor_cutoff_cm,
            "scale_factor": self.frequency_scale,
            "zpe_scale_factor": self.zpe_scale,
            "invert_imags": self.invert_imag_from_cm,
            "cutoff": self.positive_frequency_floor_cm,
        }

    def as_dict(self) -> dict:
        """Return the serialized policy with its stable public field keys."""
        return {
            "kind": self.kind,
            "rotor_cutoff_cm": self.rotor_cutoff_cm,
            "frequency_scale": self.frequency_scale,
            "zpe_scale": self.zpe_scale,
            "invert_imag_from_cm": self.invert_imag_from_cm,
            "positive_frequency_floor_cm": self.positive_frequency_floor_cm,
        }


# Standalone ``mlmm freq``: library-default QRRHO with a 100 cm^-1 rotor
# cutoff, unit frequency/ZPE scaling, NO imaginary inversion and NO
# positive-frequency floor. Reproduces the historical bare
# ``thermochemistry(qc, T, pressure=p)`` numbers exactly.
WORKFLOW_THERMO_POLICY = ThermoPolicy()

# Bundled ``Geometry.get_thermoanalysis``: the historical behaviour that inverts
# small imaginaries down to -15 cm^-1 and floors positive frequencies below
# 25 cm^-1; otherwise identical to the workflow policy.
GEOMETRY_THERMO_POLICY = ThermoPolicy(
    invert_imag_from_cm=-15.0,
    positive_frequency_floor_cm=25.0,
)
