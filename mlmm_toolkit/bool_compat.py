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


def normalize_bool_argv(
    args: list[str],
    bool_value_options_by_command: dict[str, frozenset[str]],
    bool_toggle_options_by_command: dict[str, frozenset[str]],
    toggle_negative_aliases: dict[str, dict[str, str]],
) -> tuple[list[str], bool]:
    """Normalize CLI argv booleans and return (normalized_args, legacy_syntax_used)."""
    if not args:
        return args, False

    command = args[0]
    bool_value_options = bool_value_options_by_command.get(command)
    bool_toggle_options = bool_toggle_options_by_command.get(command)
    if not bool_value_options and not bool_toggle_options:
        return args, False

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
        if name.startswith("--no-"):
            positive_name = "--" + name[5:]
            if sep == "" and bool_value_options and positive_name in bool_value_options:
                normalized.extend([positive_name, "False"])
                i += 1
                continue
            if bool_toggle_options and positive_name in bool_toggle_options:
                if sep:
                    parsed_inline = _parse_bool_literal(inline_value)
                    if parsed_inline is not None:
                        legacy_used = True
                        normalized.append(name if parsed_inline else positive_name)
                    else:
                        normalized.append(token)
                    i += 1
                    continue
                if i + 1 < len(args):
                    parsed_next = _parse_bool_literal(args[i + 1])
                    if parsed_next is not None:
                        legacy_used = True
                        normalized.append(name if parsed_next else positive_name)
                        i += 2
                        continue
                normalized.append(name)
                i += 1
                continue
            normalized.append(token)
            i += 1
            continue

        if bool_value_options and name in bool_value_options:
            if sep:
                parsed_inline = _parse_bool_literal(inline_value)
                if parsed_inline is not None:
                    legacy_used = True
                    normalized.extend([name, "True" if parsed_inline else "False"])
                else:
                    normalized.append(token)
                i += 1
                continue

            if i + 1 < len(args):
                parsed_next = _parse_bool_literal(args[i + 1])
                if parsed_next is not None:
                    legacy_used = True
                    normalized.extend([name, "True" if parsed_next else "False"])
                    i += 2
                    continue

            normalized.extend([name, "True"])
            i += 1
            continue

        if bool_toggle_options and name in bool_toggle_options:
            if sep:
                parsed_inline = _parse_bool_literal(inline_value)
                if parsed_inline is not None:
                    legacy_used = True
                    normalized.append(
                        name
                        if parsed_inline
                        else _toggle_negative_name(
                            command, name, toggle_negative_aliases
                        )
                    )
                else:
                    normalized.append(token)
                i += 1
                continue

            if i + 1 < len(args):
                parsed_next = _parse_bool_literal(args[i + 1])
                if parsed_next is not None:
                    legacy_used = True
                    normalized.append(
                        name
                        if parsed_next
                        else _toggle_negative_name(
                            command, name, toggle_negative_aliases
                        )
                    )
                    i += 2
                    continue

            normalized.append(name)
            i += 1
            continue

        normalized.append(token)
        i += 1

    return normalized, legacy_used
