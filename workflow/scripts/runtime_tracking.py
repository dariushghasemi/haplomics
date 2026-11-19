"""Collect changed Snakemake benchmarks and workflow elapsed time using hooks.

No third-party dependencies. Benchmark paths must be benchmarks/<rule>/<job>.tsv.
One workflow invocation at a time per output workspace, with shared filesystem
visibility for benchmarks. Skipped jobs retain their previous benchmark files
but are excluded from the current run's report.
"""
from collections import defaultdict
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import time
from uuid import uuid4
import warnings


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class RuntimeTracker:
    def __init__(self, benchmark_dir, report_dir):
        self.benchmark_dir = Path(benchmark_dir)
        self.report_dir = Path(report_dir)
        self.started = None

    def snapshot(self):
        files = {}
        for path in self.benchmark_dir.glob("*/*.tsv"):
            stat = path.stat()
            files[path] = (stat.st_mtime_ns, stat.st_ctime_ns, stat.st_size)
        return files

    def start(self):
        self.started = time.perf_counter()
        self.started_at = timestamp()
        self.before = self.snapshot()
        run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ") + "_" + uuid4().hex[:8]
        self.run_dir = self.report_dir / run_id
        self.run_dir.mkdir(parents=True, exist_ok=False)
        # Preserve evidence of an interrupted run even if no finish hook executes.
        self.write_run("running", None, [], [])

    def write_run(self, status, elapsed, jobs, errors):
        report = {
            "status": status,
            "started_at_utc": self.started_at,
            "finished_at_utc": timestamp() if elapsed is not None else None,
            "total_elapsed_seconds": elapsed,
            "benchmarked_job_seconds": sum(job["seconds"] for job in jobs),
            "benchmark_records": len(jobs),
            "benchmark_read_errors": errors,
            "scope": "Elapsed time from onstart to finish hook; job times only for changed benchmarks in this invocation.",
        }
        (self.run_dir / "run.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    def finish(self, status):
        if self.started is None:
            return  # e.g. a failure before onstart
        elapsed = time.perf_counter() - self.started
        jobs, errors = [], []
        for path, signature in sorted(self.snapshot().items()):
            if self.before.get(path) == signature:
                continue
            try:
                with path.open(encoding="utf-8", newline="") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                if not rows:
                    raise ValueError("empty benchmark")
                measured = []
                for repeat, row in enumerate(rows, 1):
                    seconds = float(row["s"])
                    if not math.isfinite(seconds) or seconds < 0:
                        raise ValueError("invalid benchmark seconds")
                    measured.append(dict(rule=path.parent.name, job=path.stem,
                                         repeat=repeat, seconds=seconds,
                                         benchmark=str(path)))
                jobs.extend(measured)
            except (OSError, KeyError, TypeError, ValueError, csv.Error) as error:
                errors.append(f"{path}: {error}")

        grouped = defaultdict(list)
        for job in jobs:
            grouped[job["rule"]].append(job["seconds"])
        summary = [dict(rule=rule, measurements=len(values), total_seconds=sum(values),
                        mean_seconds=sum(values)/len(values), min_seconds=min(values),
                        max_seconds=max(values)) for rule, values in sorted(grouped.items())]
        write_tsv(self.run_dir / "jobs.tsv", ["rule", "job", "repeat", "seconds", "benchmark"], jobs)
        write_tsv(self.run_dir / "rules.tsv",
                  ["rule", "measurements", "total_seconds", "mean_seconds", "min_seconds", "max_seconds"], summary)
        self.write_run(status, elapsed, jobs, errors)
        for error in errors:
            warnings.warn(f"Runtime report: {error}")
        print(f"Runtime report: {self.run_dir} (workflow elapsed: {elapsed:.3f} s; {status})")
