# Inspect graphics-facing matrix memory

**Intermediate.** Prerequisites: [matrix transforms](transformations.md) and basic
Python object ownership. Learn what the ctypes snapshot contains, why row-major
memory is separate from vector convention and how to reason about upload layout
without an OpenGL context.

## Rows, bytes and mathematical action

For Matrix4, `.matrix` holds Python numeric rows (floats are binary64 on the
verified environment), while `.c_matrix` is a nested ctypes.c_float (binary32)
array snapshot. Offset 4*row+column indexes a
flattened row-major buffer. XYZ translation appears at offsets 12,13,14 because
gem uses row-vector transforms with translation in the final row.

Memory order and mathematical convention are separate choices. If a consumer
interprets the same flattened bytes as a column-major matrix, it reads M^T.
Applying that with column vectors gives M^T*p_column, the transpose of gem's
p_row*M. This can be mathematically compatible without another transpose,
but only when the consumer's full upload/shader conventions agree. A caller
that explicitly transposes bytes or selects an upload transpose flag must account
for both steps; no universal flag setting is prescribed here.

## Verify the snapshot and its lifetime

```python
import ctypes
from gem.matrix import Matrix
from gem.vector import Vector

matrix = Matrix(4).translate(Vector(3, [2, -3, 4]))
owner = matrix.c_matrix
assert ctypes.sizeof(ctypes.c_float) == 4  # executed environment assumption
assert ctypes.sizeof(owner) == 16 * 4
pointer = ctypes.cast(owner, ctypes.POINTER(ctypes.c_float))
flat = [pointer[i] for i in range(16)]
assert flat[12:16] == [2.0, -3.0, 4.0, 1.0]
assert (matrix * Vector(4, [1, 2, 3, 1])).vector == [3.0, -1.0, 7.0, 1.0]
assert matrix.i_translate(Vector(3, [1, 0, 0])) is matrix
assert matrix.c_matrix is not owner
new_pointer = ctypes.cast(matrix.c_matrix, ctypes.POINTER(ctypes.c_float))
assert [new_pointer[i] for i in range(12, 16)] == [3.0, -3.0, 4.0, 1.0]
assert [pointer[i] for i in range(12, 16)] == [2.0, -3.0, 4.0, 1.0]  # owner still alive
matrix.matrix[3][0] = 99
assert matrix.c_matrix[3][0] == 3.0   # direct row edits do not refresh snapshots
print(flat[12:16])
```

Output: `[2.0, -3.0, 4.0, 1.0]`. The point moved by the final-row translation.
A supported in-place method then produced a new snapshot, so an old pointer still
refers to the old array. Holding `owner` keeps that buffer alive; reacquire the
pointer when new data is required. Direct row edits deliberately demonstrate that
the snapshot is not a live view, not a recommended synchronization mechanism.

## External graphics applications and limits

An OpenGL-facing application can pass an appropriately owned float buffer to its
binding, provided its expected ordering, transpose argument, shader multiplication
and context requirements are confirmed. This guide makes no graphics call and
verifies no driver/binding; gem does not render. Keep arrays alive until the
consumer has finished accessing them, especially for asynchronous foreign code.

Python float precision/range is not preserved in float32 export: values may round
or overflow while still representable in Python rows. Native ctypes memory is
not a portable endian-tagged file format. The 4-byte c_float assumption is checked
on the executed environment, not advertised as a newly verified platform matrix.
For fresh explicit arrays, use common.conv_list/conv_list_2d; these require no GPU.

See [Matrix ctypes contract](../api/matrix.md#numerical-domain-and-ctypes-snapshots),
[common conversion helpers](../api/common.md#ctypes-ownership),
[numerical accuracy](numerical.md), [camera conventions](camera.md) and [index](index.md).
