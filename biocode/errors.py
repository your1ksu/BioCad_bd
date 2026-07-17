"""Иерархия ошибок BIOCODE core.

Все ошибки пакета наследуются от :class:`BiocodeError`, чтобы вызывающий код мог
ловить их одним `except BiocodeError`. Никаких «тихих» падений — каждая ошибка
несёт человекочитаемое сообщение и, где уместно, подсказку по исправлению.
"""
from __future__ import annotations


class BiocodeError(Exception):
    """Базовый класс всех ошибок пакета."""


class ConfigError(BiocodeError):
    """Некорректная конфигурация прогона."""


class InputError(BiocodeError):
    """Проблема со входными данными (файл/формат/содержимое)."""


class ValidationError(InputError):
    """Данные не прошли проверки валидатора (см. validate.py)."""


class ToolNotFoundError(BiocodeError):
    """Внешний инструмент не найден.

    Сообщение включает команду установки, чтобы пользователь мог сразу исправить.
    """

    def __init__(self, tool: str, install_hint: str = "") -> None:
        msg = f"Инструмент '{tool}' не найден (ни в конфиге, ни в $PATH, ни в conda-env)."
        if install_hint:
            msg += f" Установка: {install_hint}"
        super().__init__(msg)
        self.tool = tool


class ToolRunError(BiocodeError):
    """Внешний инструмент завершился с ошибкой/таймаутом."""

    def __init__(self, tool: str, cmd: list[str], returncode: int | None,
                 log_tail: str = "") -> None:
        rc = "timeout" if returncode is None else f"код {returncode}"
        msg = f"'{tool}' завершился с ошибкой ({rc}).\n  команда: {' '.join(cmd)}"
        if log_tail:
            msg += "\n  --- хвост лога ---\n" + log_tail
        super().__init__(msg)
        self.tool = tool
        self.cmd = cmd
        self.returncode = returncode
