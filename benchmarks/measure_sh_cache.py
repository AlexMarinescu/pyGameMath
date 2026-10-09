"""Measure bounded SH-layout cache cost separately from warm timings."""
import json
from pathlib import Path
import sys
import tracemalloc

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from gem import spherical_harmonics as sh


def main():
    sh._basis_layout.cache_clear();tracemalloc.start()
    sh._basis(3,.73,1.27)
    retained3,peak3=tracemalloc.get_traced_memory()
    tracemalloc.reset_peak()
    for bands in range(1,17):sh._basis(bands,.73,1.27)
    retained16,peak16=tracemalloc.get_traced_memory()
    info=sh._basis_layout.cache_info();tracemalloc.stop()
    result={'cold_l2_layout_and_basis_peak_bytes':peak3,'l2_layout_retained_traced_bytes':retained3,
            'layouts_1_through_16_retained_traced_bytes':retained16,'fill_16_layouts_peak_bytes':peak16,
            'cache':{'maxsize':info.maxsize,'currsize':info.currsize,'hits':info.hits,'misses':info.misses},
            'method':'clear cache before tracing; discard each result; peak includes recurrence/output temporaries',
            'limitations':'tracemalloc Python allocations only, not RSS/native allocations; at most 16 layouts, each O(bands^2)'}
    Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n');print(result)


if __name__=='__main__':main()
