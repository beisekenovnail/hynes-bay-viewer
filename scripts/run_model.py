"""Run the yearly runs in order; the end state of each year starts the next (hotstart).
Usage: python3 run_model.py [runs folder]. Env: YEARS (default: every year folder)."""
import os, shutil, subprocess, sys
from config import years, RUNS, BIN


def run_chain(runs, yrs, bin=BIN):
    bin = os.path.abspath(bin)
    for y in yrs:
        d = f"{runs}/{y}"
        if y != yrs[0]: shutil.copy(f"{runs}/{y - 1}/hotstart_out.dat", f"{d}/hotstart_in.dat")
        with open(f"{d}/run_stdout.txt", "w") as fo:
            rc = subprocess.run([bin], cwd=d, stdout=fo, stderr=subprocess.STDOUT).returncode
        n = sum(1 for _ in open(f"{d}/SAL.out")) if os.path.exists(f"{d}/SAL.out") else 0
        want = sum(1 for _ in open(f"{d}/Precip.csv")) - 1
        print(f"  {y}: exit {rc}, {n} of {want} days", flush=True)
        if rc != 0 or n < want: return False
    return True


if __name__ == "__main__":
    runs = sys.argv[1] if len(sys.argv) > 1 else RUNS
    yrs = years("") or sorted(int(d) for d in os.listdir(runs) if d.isdigit())
    print("RUN", "OK" if run_chain(runs, yrs) else "FAILED")
