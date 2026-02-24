# Path and File Name Handling Issues Report

This document identifies issues with path and file name handling across the codebase that should be updated for better cross-platform compatibility, maintainability, and modern Python practices.

## Summary of Issues

### 1. **Hardcoded Path Separators** (Critical - Cross-platform issue)
Using hardcoded `/` or `\` instead of `os.path.join()` or `pathlib.Path`

### 2. **String Concatenation for Paths** (Critical - Cross-platform issue)
Using `+` operator to build file paths instead of proper path joining

### 3. **String Splitting for Path Extraction** (Moderate - Fragile)
Manually splitting paths with `.split("/")` instead of using pathlib or os.path functions

### 4. **Nested os.path.join()** (Minor - Redundant)
Unnecessary nesting of `os.path.join()` calls

### 5. **os.mkdir() instead of os.makedirs()** (Moderate - Error-prone)
Using `os.mkdir()` which fails if parent directories don't exist

### 6. **F-string Path Construction** (Minor - Inconsistent)
Mixing f-strings with path construction instead of consistent use of os.path.join

### 7. **Missing Path Validation** (Moderate - Error-prone)
No validation that paths exist or are valid before use

---

## Detailed Issues by File

### `pyext/src/analysis_trajectories.py`

#### Line 128: `os.mkdir()` should be `os.makedirs()`
```python
# Current:
os.mkdir(self.analysis_dir)

# Should be:
os.makedirs(self.analysis_dir, exist_ok=True)
```
**Issue**: Will fail if parent directories don't exist. Should use `exist_ok=True` to avoid errors if directory already exists.

#### Line 597: Hardcoded path separator
```python
# Current:
traj = [x for x in out.split("/") if self.dir_name in x][0]

# Should use:
from pathlib import Path
traj = [x for x in Path(out).parts if self.dir_name in x][0]
# Or:
traj = os.path.basename(out) if self.dir_name in out else None
```

#### Lines 822, 831, 840: Hardcoded path separator and string splitting
```python
# Current:
kk = k.split(self.dir_name)[-1].split("/")[0]

# Should use:
from pathlib import Path
kk = Path(k).name.split(self.dir_name)[-1] if self.dir_name in k else Path(k).name
```

#### Lines 860, 871: String splitting for filename extraction
```python
# Current:
k = f.split("scores_info_")[-1].split(".csv")[0]

# Should use:
from pathlib import Path
k = Path(f).stem.replace("scores_info_", "")
```

#### Lines 1000, 1004: Hardcoded path separator
```python
# Current:
runs_A.append([x for x in t.split("/") if self.dir_name in x][0])

# Should use:
from pathlib import Path
runs_A.append([x for x in Path(t).parts if self.dir_name in x][0])
```

#### Lines 824-826, 833-835, 842-844: Nested os.path.join()
```python
# Current:
T.to_csv(
    os.path.join(
        os.path.join(self.analysis_dir, f"scores_info_{kk}.csv")
    )
)

# Should be:
T.to_csv(os.path.join(self.analysis_dir, f"scores_info_{kk}.csv"))
```

#### Line 1097-1099: String concatenation for filename
```python
# Current:
"h1_" + row.traj + "_" + str(int(row.MC_frame)) + ".rmf3"

# Should use:
f"h1_{row.traj}_{int(row.MC_frame)}.rmf3"
# Or better, use os.path.join for directory + filename separately
```

#### Line 1222: String concatenation for filename
```python
# Current:
os.path.join(gsms_dir, filename + ".txt")

# Should be:
os.path.join(gsms_dir, f"{filename}.txt")
```

#### Line 1333: String concatenation for filename
```python
# Current:
os.path.join(analysis_dir, scores_prefix + ".txt")

# Should be:
os.path.join(analysis_dir, f"{scores_prefix}.txt")
```

#### Line 1427: `os.mkdir()` should be `os.makedirs()`
```python
# Current:
os.mkdir(d)

# Should be:
os.makedirs(d, exist_ok=True)
```

#### Line 1994: Typo in variable name
```python
# Current:
file_out_hist = f"plot_XLs_histogram_cluster{cluster}_{types_XLs}.{self.plot_fmt}"
# Note: types_XLs should be type_XLs
```

---

### `pyext/src/align_rmf.py`

#### Lines 112-115: Hardcoded path separator
```python
# Current:
if "/" in rmf_in:
    name_in = rmf_in.split("/")[-1].split(".")[0]
else:
    name_in = rmf_in.split(".")[0]

# Should use:
from pathlib import Path
name_in = Path(rmf_in).stem
```

#### Line 116: String concatenation for filename
```python
# Current:
out_RMSD = f"{RMSD}_{name_in}.txt"

# Should use os.path.join if this is meant to be a full path:
# Or keep as f-string if it's just a filename
```

---

### `pyext/src/contact_maps.py`

#### Line 152: F-string with hardcoded separator
```python
# Current:
file_name = f"{self.clustering_dir}/cluster.{self.cluster}.sample_{half}.txt"

# Should be:
file_name = os.path.join(
    self.clustering_dir, 
    f"cluster.{self.cluster}.sample_{half}.txt"
)
```

#### Line 639: String splitting for filename extraction
```python
# Current:
prot = f.split("ContMap_")[-1].split(".dat")[0]

# Should use:
from pathlib import Path
prot = Path(f).stem.replace("ContMap_", "")
```

---

### `pyext/src/validation.py`

#### Line 199: String concatenation for path
```python
# Current:
dist = pd.read_csv(self.analysis_dir + '/XLs_dist_info_'
                   + str(t) + '.csv')

# Should be:
dist = pd.read_csv(
    os.path.join(self.analysis_dir, f"XLs_dist_info_{t}.csv")
)
```

#### Lines 286, 289: String concatenation with hardcoded separator
```python
# Current:
stats_XLs.to_csv(
    os.path.join(self.clustering_dir,
                 '/XLs_satisfaction_cluster_'+str(cluster)+'.csv'))
# Note: Leading '/' makes this an absolute path!

# Should be:
stats_XLs.to_csv(
    os.path.join(self.clustering_dir,
                 f"XLs_satisfaction_cluster_{cluster}.csv"))
```

---

### `pyext/src/compute_distance_metrics.py`

#### Lines 137-138: String concatenation for path
```python
# Current:
rmf3 = self.clustering_dir + '/cluster.' + \
    str(self.cluster) + '/cluster_center_model.rmf3'

# Should be:
rmf3 = os.path.join(
    self.clustering_dir,
    f"cluster.{self.cluster}",
    "cluster_center_model.rmf3"
)
```

#### Line 59: String concatenation in loop
```python
# Current:
name += str(vv[0])+' '+str(vv[1])+'/'

# Should use:
name_parts = [f"{vv[0]} {vv[1]}" for vv in v]
name = "/".join(name_parts)
```

---

### `pyext/src/accuracy.py`

#### Lines 154, 164: String concatenation for filenames
```python
# Current:
"accuracy_" + self.out_header + "_clusters.dat"
"accuracy_" + self.out_header + "_cl" + str(k) + ".dat"

# Should be:
f"accuracy_{self.out_header}_clusters.dat"
f"accuracy_{self.out_header}_cl{k}.dat"
```

---

### `example/run_extract_models.py`

#### Lines 36-38: F-string with hardcoded separator
```python
# Current:
f'{analysis_dir}/selected_models_A_cluster{cluster}_detailed.csv'

# Should be:
os.path.join(analysis_dir, f'selected_models_A_cluster{cluster}_detailed.csv')
```

---

### `example/rerun_clustering.py`

#### Line 16: String concatenation for path
```python
# Current:
analys_dir = top_dir+'/analys/'

# Should be:
analys_dir = os.path.join(top_dir, 'analys')
```

#### Line 20: String concatenation for path
```python
# Current:
out_dirs = glob.glob(top_dir+'/'+dir_head+'*/output/')

# Should be:
out_dirs = glob.glob(os.path.join(top_dir, f"{dir_head}*/output"))
```

---

## Recommendations

### 1. **Use `pathlib.Path` for Modern Python (Python 3.4+)**
   - More intuitive and object-oriented
   - Better cross-platform support
   - Cleaner code

### 2. **Use `os.makedirs(..., exist_ok=True)` instead of `os.mkdir()`**
   - Handles missing parent directories
   - Prevents errors if directory exists

### 3. **Always use `os.path.join()` or `pathlib.Path` for path construction**
   - Never use string concatenation with `+`
   - Never use hardcoded separators

### 4. **Use `Path.stem`, `Path.name`, `Path.parent` for path manipulation**
   - Instead of string splitting
   - More robust and readable

### 5. **Add path validation**
   - Check if paths exist before use
   - Use `os.path.exists()` or `Path.exists()`
   - Provide clear error messages

### 6. **Consistent filename construction**
   - Use f-strings for filenames: `f"{prefix}_{suffix}.ext"`
   - Use `os.path.join()` for combining directory + filename

---

## Priority Order for Fixes

1. **High Priority**: Hardcoded separators and string concatenation (cross-platform issues)
2. **Medium Priority**: `os.mkdir()` → `os.makedirs()` (error prevention)
3. **Medium Priority**: String splitting for path extraction (fragility)
4. **Low Priority**: Nested `os.path.join()` (code cleanup)
5. **Low Priority**: F-string path construction (consistency)

---

## Example Modernization

### Before:
```python
traj = [x for x in out.split("/") if self.dir_name in x][0]
kk = k.split(self.dir_name)[-1].split("/")[0]
file_path = self.analysis_dir + "/" + filename + ".csv"
os.mkdir(self.analysis_dir)
```

### After:
```python
from pathlib import Path

traj = next((x for x in Path(out).parts if self.dir_name in x), None)
kk = Path(k).name.split(self.dir_name)[-1] if self.dir_name in k else Path(k).name
file_path = Path(self.analysis_dir) / f"{filename}.csv"
Path(self.analysis_dir).mkdir(parents=True, exist_ok=True)
```

Or using `os.path`:
```python
import os

traj = os.path.basename(out) if self.dir_name in out else None
kk = os.path.basename(k).split(self.dir_name)[-1] if self.dir_name in k else os.path.basename(k)
file_path = os.path.join(self.analysis_dir, f"{filename}.csv")
os.makedirs(self.analysis_dir, exist_ok=True)
```

