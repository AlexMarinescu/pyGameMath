"""Compatibility imports; use gem.bezier for supported evaluation.

Adaptive sampling remains experimental and retains known E04 defects.
"""
from gem.bezier import cubicBezierPoint, quadraticBezierPoint
from gem.experimental._bezier_legacy import BezierPath

__all__ = ["cubicBezierPoint", "quadraticBezierPoint", "BezierPath"]
