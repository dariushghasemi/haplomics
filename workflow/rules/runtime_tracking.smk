# Include once, after ws_path is available. Merge these calls into any existing
# workflow hooks instead of registering a second set of handlers.

import runpy

_haplotype_runtime = runpy.run_path("workflow/scripts/runtime_tracking.py")["RuntimeTracker"](
    ws_path("benchmarks"), ws_path("logs/runtime")
    )

onstart:
    _haplotype_runtime.start()
    
onsuccess:
    _haplotype_runtime.finish("success")

onerror:
    _haplotype_runtime.finish("failed")
