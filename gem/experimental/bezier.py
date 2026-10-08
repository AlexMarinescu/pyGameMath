"""Compatibility imports; use gem.bezier for supported Bezier operations."""
from gem.bezier import BezierPath, cubicBezierPoint, quadraticBezierPoint

__all__ = ["cubicBezierPoint", "quadraticBezierPoint", "BezierPath"]
