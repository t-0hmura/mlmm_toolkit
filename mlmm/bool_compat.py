"""CLI bool-argument normalization helpers.

This module keeps backward compatibility for legacy value-style boolean flags:
`--flag True/False`.
"""

from __future__ import annotations


_BOOL_TRUE_LITERALS = {"1", "true", "t", "yes", "y", "on"}
_BOOL_FALSE_LITERALS = {"0", "false", "f", "no", "n", "off"}


def _parse_bool_literal(raw: str) -> bool | None:
    token = raw.strip().lower()
    if token in _BOOL_TRUE_LITERALS:
        return True
    if token in _BOOL_FALSE_LITERALS:
        return False
    return None


def _toggle_negative_name(
    command: str,
    positive_name: str,
    toggle_negative_aliases: dict[str, dict[str, str]],
) -> str:
    aliases = toggle_negative_aliases.get(command)
    if aliases and positive_name in aliases:
        return aliases[positive_name]
    return f"--no-{positive_name[2:]}"


def _read_bool_literal(
    args: list[str], index: int, sep: str, inline_value: str
) -> tuple[bool | None, int, bool]:
    """Return (parsed_bool, next_index, explicit_literal_used)."""
    if sep:
        return _parse_bool_literal(inline_value), index + 1, True
    if index + 1 < len(args):
        parsed_next = _parse_bool_literal(args[index + 1])
        if parsed_next is not None:
            return parsed_next, index + 2, True
    return None, index + 1, False


def _build_negative_to_positive_toggle_map(
    toggle_options: frozenset[str],
    toggle_negative_aliases: dict[str, dict[str, str]],
    command: str,
) -> dict[str, str]:
    negative_to_positive: dict[str, str] = {}
    for positive_name in toggle_options:
        if not positive_name.startswith("--"):
            continue
        canonical_negative = _toggle_negative_name(
            command, positive_name, toggle_negative_aliases
        )
        negative_to_positive[canonical_negative] = positive_name
        synthetic_no_name = f"--no-{positive_name[2:]}"
        negative_to_positive.setdefault(synthetic_no_name, positive_name)
    return negative_to_positive


def _build_negative_to_positive_single_flag_map(
    single_flag_options: frozenset[str],
) -> dict[str, str]:
    negative_to_positive: dict[str, str] = {}
    for positive_name in single_flag_options:
        if not positive_name.startswith("--"):
            continue
        negative_to_positive[f"--no-{positive_name[2:]}"] = positive_name
    return negative_to_positive


def normalize_bool_argv(
    args: list[str],
    bool_value_options_by_command: dict[str, frozenset[str]],
    bool_toggle_options_by_command: dict[str, frozenset[str]],
    toggle_negative_aliases: dict[str, dict[str, str]],
    bool_single_flag_options_by_command: dict[str, frozenset[str]] | None = None,
) -> tuple[list[str], bool]:
    """Normalize CLI argv booleans and return (normalized_args, legacy_syntax_used)."""
    if not args:
        return args, False

    command = args[0]
    bool_value_options = bool_value_options_by_command.get(command, frozenset())
    bool_toggle_options = bool_toggle_options_by_command.get(command, frozenset())
    bool_single_flag_options = (
        bool_single_flag_options_by_command or {}
    ).get(command, frozenset())

    if not bool_value_options and not bool_toggle_options and not bool_single_flag_options:
        return args, False

    toggle_negative_to_positive = _build_negative_to_positive_toggle_map(
        bool_toggle_options, toggle_negative_aliases, command
    )
    single_negative_to_positive = _build_negative_to_positive_single_flag_map(
        bool_single_flag_options
    )

    normalized: list[str] = [command]
    legacy_used = False
    i = 1
    while i < len(args):
        token = args[i]

        if token == "--":
            normalized.extend(args[i:])
            break

        if not token.startswith("--"):
            normalized.append(token)
            i += 1
            continue

        name, sep, inline_value = token.partition("=")

        if name in toggle_negative_to_positive:
            positive_name = toggle_negative_to_positive[name]
            canonical_negative = _toggle_negative_name(
                command, positive_name, toggle_negative_aliases
            )
            parsed, next_i, explicit_literal = _read_bool_literal(
                args, i, sep, inline_value
            )
            if explicit_literal:
                if parsed is None:
                    normalized.append(token)
                else:
                    legacy_used = True
                    normalized.append(canonical_negative if parsed else positive_name)
                i = next_i
                continue

            normalized.append(canonical_negative)
            i += 1
            continue

        if name in single_negative_to_positive:
            positive_name = single_negative_to_positive[name]
            parsed, next_i, explicit_literal = _read_bool_literal(
                args, i, sep, inline_value
            )
            if explicit_literal:
                if parsed is None:
                    normalized.append(token)
                else:
                    legacy_used = True
                    if not parsed:
                        normalized.append(positive_name)
                i = next_i
                continue

            i += 1
            continue

        if name.startswith("--no-"):
            positive_name = "--" + name[5:]
            if sep == "" and positive_name in bool_value_options:
                normalized.extend([positive_name, "False"])
                i += 1
                continue
            normalized.append(token)
            i += 1
            continue

        if name in bool_value_options:
            parsed, next_i, explicit_literal = _read_bool_literal(args, i, sep, inline_value)
            if explicit_literal:
                if parsed is not None:
                    legacy_used = True
                    normalized.extend([name, "True" if parsed else "False"])
                else:
                    normalized.append(token)
                i = next_i
                continue

            normalized.extend([name, "True"])
            i += 1
            continue

        if name in bool_toggle_options:
            parsed, next_i, explicit_literal = _read_bool_literal(args, i, sep, inline_value)
            if explicit_literal:
                if parsed is not None:
                    legacy_used = True
                    normalized.append(
                        name if parsed else _toggle_negative_name(
                            command, name, toggle_negative_aliases
                        )
                    )
                else:
                    normalized.append(token)
                i = next_i
                continue

            normalized.append(name)
            i += 1
            continue

        if name in bool_single_flag_options:
            parsed, next_i, explicit_literal = _read_bool_literal(args, i, sep, inline_value)
            if explicit_literal:
                if parsed is not None:
                    legacy_used = True
                    if parsed:
                        normalized.append(name)
                else:
                    normalized.append(token)
                i = next_i
                continue

            normalized.append(name)
            i += 1
            continue

        normalized.append(token)
        i += 1

    return normalized, legacy_used
