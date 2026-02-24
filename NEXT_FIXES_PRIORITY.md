# Next Fixes - Priority List

Based on the codebase analysis, here are the recommended next improvements, ordered by priority:

---

## 🔴 **HIGH PRIORITY - Security & Critical Issues**

### 1. **Replace `eval()` with `ast.literal_eval()` or `json.loads()`** ⚠️ SECURITY RISK
**Files:** `pyext/src/analysis_trajectories.py` (lines 397, 482)

**Issue:** Using `eval()` on file input is a **security vulnerability**. If an attacker can modify the input file, they can execute arbitrary code.

**Current code:**
```python
d = eval(line)  # DANGEROUS!
```

**Should be:**
```python
import ast
d = ast.literal_eval(line)  # Safe - only evaluates literals
# OR if the file is JSON:
import json
d = json.loads(line)
```

**Impact:** Critical security fix - prevents code injection attacks.

---

### 2. **Replace deprecated pandas `.append()` method**
**Files:** `pyext/src/validation.py` (line 415)

**Issue:** `DataFrame.append()` is deprecated in pandas 2.0+ and will be removed. It's also inefficient.

**Current code:**
```python
sPEMAP_all = sPEMAP_all.append(
    {'cluster': cl, 'pEMAP_satisfaction': percent_satif},
    ignore_index=True)
```

**Should be:**
```python
# Option 1: Use pd.concat (recommended)
new_row = pd.DataFrame([{'cluster': cl, 'pEMAP_satisfaction': percent_satif}])
sPEMAP_all = pd.concat([sPEMAP_all, new_row], ignore_index=True)

# Option 2: Build list and create DataFrame once (more efficient)
# Collect rows in a list, then create DataFrame at the end
```

**Impact:** Prevents future pandas compatibility issues, improves performance.

---

## 🟡 **MEDIUM PRIORITY - Code Quality & Maintainability**

### 3. **Replace bare `except:` clauses with specific exceptions**
**Files:** Multiple files (13 instances found)

**Issue:** Bare `except:` clauses catch all exceptions including `KeyboardInterrupt` and `SystemExit`, making debugging difficult.

**Current code:**
```python
except:  # noqa: E722
    print("No S_info")
```

**Should be:**
```python
except (KeyError, AttributeError, IndexError) as e:
    print(f"No S_info: {e}")
    # Or log it properly
```

**Locations:**
- `analysis_trajectories.py`: lines 483, 759, 1081, 1704, 1715
- `tools.py`: lines 49, 57

**Impact:** Better error handling, easier debugging, prevents catching unintended exceptions.

---

### 4. **Fix hardcoded system paths in example files**
**Files:** All example files

**Issue:** Hardcoded absolute paths prevent code from working on other systems.

**Current code:**
```python
sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
```

**Should be:**
```python
import os
from pathlib import Path

# Get the directory of this script
script_dir = Path(__file__).parent.parent.parent
sys.path.insert(0, str(script_dir / 'pyext' / 'src'))
# OR use relative imports if possible
```

**Files affected:**
- `example/rerun_clustering.py`
- `example/run_extract_models.py`
- `example/run_analysis_trajectories.py`
- `example/get_accuracy.py`
- `example/get_CMs.py`

**Impact:** Makes code portable and usable by others.

---

### 5. **Standardize `__future__` imports**
**Issue:** Some files have `from __future__ import print_function`, others don't. Since Python 3 is required, these are mostly unnecessary, but consistency matters.

**Files with `__future__`:**
- `contact_maps.py`
- `validation.py`
- `accuracy.py`
- `tools.py`
- `equilibration.py`

**Files without:**
- `analysis_trajectories.py` (main file)
- `align_rmf.py`
- `compute_distance_metrics.py`

**Recommendation:** 
- Remove `from __future__ import print_function` (not needed in Python 3)
- Keep `from __future__ import division` if needed for integer division behavior
- Or add consistently if you want to maintain Python 2 compatibility (unlikely)

**Impact:** Code consistency, cleaner imports.

---

## 🟢 **LOW PRIORITY - Code Quality Improvements**

### 6. **Replace `print()` statements with proper logging**
**Issue:** 70+ print statements throughout the codebase. Should use Python's `logging` module for better control.

**Current code:**
```python
print("Trajectory, ts_eqs: ", traj, ts_eq)
```

**Should be:**
```python
import logging
logger = logging.getLogger(__name__)
logger.info("Trajectory, ts_eqs: %s, %s", traj, ts_eq)
```

**Impact:** Better control over output, can redirect to files, different log levels, etc.

---

### 7. **Add type hints**
**Issue:** No type hints in the codebase. Would improve IDE support and code documentation.

**Example:**
```python
# Before:
def get_keys(self, stat_file):

# After:
from typing import Dict, Any
def get_keys(self, stat_file: str) -> Dict[str, Any]:
```

**Impact:** Better IDE support, self-documenting code, catch errors earlier.

---

### 8. **Fix string concatenation in validation.py**
**Issue:** Still some string concatenation for filenames.

**Location:** `validation.py` line 57-59
```python
# Current:
'h1_'+row.traj+'_'+str(int(row.MC_frame))+'.rmf3'

# Should be:
f'h1_{row.traj}_{int(row.MC_frame)}.rmf3'
```

**Impact:** Consistency with other fixes, readability.

---

### 9. **Remove duplicate imports**
**Issue:** Some files import the same module twice.

**Example:** `validation.py` imports `IMP` twice (lines 2 and 16)

**Impact:** Code cleanliness.

---

### 10. **Add docstrings to missing functions**
**Issue:** Some functions lack proper docstrings.

**Impact:** Better code documentation and IDE help.

---

## Recommended Order of Fixes

1. **First:** Fix `eval()` security issue (#1) - **CRITICAL**
2. **Second:** Fix deprecated pandas `.append()` (#2) - **HIGH**
3. **Third:** Fix hardcoded paths in examples (#4) - **HIGH** (affects usability)
4. **Fourth:** Replace bare except clauses (#3) - **MEDIUM**
5. **Fifth:** Standardize imports (#5) - **MEDIUM**
6. **Later:** Logging, type hints, etc. (#6-10) - **LOW** (nice to have)

---

## Quick Win: Start with #1 and #2

These two fixes are:
- **Critical for security** (#1)
- **Prevents future breakage** (#2)
- **Relatively quick to implement**
- **High impact**

Would you like me to start with these?



