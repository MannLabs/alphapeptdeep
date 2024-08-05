# Starting the gui will raise by OpenMP in LLVM package
# `OMP: Error #15: Initializing libomp.dylib, but found libomp.dylib already initialized.`
# This is a quick fix, and it will only affect the GUI rather than the kernel.
import os

os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
