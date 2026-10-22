import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import generic_slice_C as GC
import density_loop_D4_factors as F4
GC.NT = 1100
F4.NT = 1100
import generic_slice_D as GD
GD.main(Bp=6, a=8, Dg=60)
