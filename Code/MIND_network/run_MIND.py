import os
import sys
import shutil
import pandas as pd
from pathlib import Path

# ================= USER CONFIGURATION =================
MIND_CODE_DIR = Path("/path/to/function/MIND/code")

DATA_ROOT = Path("/path/to/data/morphology/results")

OUTPUT_DIR = Path("/path/to/results/MIND_output")

MIND_FEATURES = ['CT', 'SA', 'Vol', 'SD', 'MC']

PARCELLATION_NAME = 'DK318'
# =======================================================

def setup_environment():
    if str(MIND_CODE_DIR) not in sys.path:
        sys.path.insert(0, str(MIND_CODE_DIR))
    try:
        from MIND import compute_MIND
        return compute_MIND
    except ImportError as e:
        print("Error importing MIND: {}".format(e))
        sys.exit(1)

def run_single_subject(subject_id):
    compute_MIND_func = setup_environment()
    
    pid = os.getpid()
    
    temp_root = OUTPUT_DIR / "temp_processing" / "{}_{}".format(subject_id, pid)
    temp_surf = temp_root / "surf"
    temp_label = temp_root / "label"
    
    src_fs_dir = DATA_ROOT / subject_id / "surf"
    src_label_dir = DATA_ROOT / subject_id / "label"

    if not src_fs_dir.exists() or not src_label_dir.exists():
        print("[SKIP] Data missing for {}".format(subject_id))
        return

    final_output_csv = OUTPUT_DIR / "{}.MIND_network.csv".format(subject_id)
    if final_output_csv.exists():
        print("[SKIP] Result already exists for {}".format(subject_id))
        return

    print("--- Processing {} (PID: {}) ---".format(subject_id, pid))
    
    try:
        if not temp_surf.exists(): os.makedirs(str(temp_surf))
        if not temp_label.exists(): os.makedirs(str(temp_label))

        for file_path in src_fs_dir.glob("*"):
            if file_path.is_file():
                try:
                    os.symlink(str(file_path), str(temp_surf / file_path.name))
                except OSError:
                    shutil.copy2(str(file_path), str(temp_surf / file_path.name))

        
        lh_name = "{}.L.{}.native.annot".format(subject_id, PARCELLATION_NAME)
        rh_name = "{}.R.{}.native.annot".format(subject_id, PARCELLATION_NAME)
        
        src_lh = src_label_dir / lh_name
        src_rh = src_label_dir / rh_name
        
        if src_lh.exists() and src_rh.exists():
            shutil.copy2(str(src_lh), str(temp_label / "lh.{}.annot".format(PARCELLATION_NAME)))
            shutil.copy2(str(src_rh), str(temp_label / "rh.{}.annot".format(PARCELLATION_NAME)))
        else:
            print("[ERROR] Label files missing for {}".format(subject_id))
            print("Expected: {} and {}".format(lh_name, rh_name))
            return


        print("Computing MIND network...")
        sys.stdout.flush()
        MIND_network = compute_MIND_func(str(temp_root), MIND_FEATURES, PARCELLATION_NAME, True)

        if not OUTPUT_DIR.exists():
            os.makedirs(str(OUTPUT_DIR))
            
        MIND_network.to_csv(str(final_output_csv), index=False)
        print("[SUCCESS] Saved to {}".format(final_output_csv.name))

    except Exception as e:
        print("[ERROR] Failed processing {}: {}".format(subject_id, e))
        if final_output_csv.exists():
            os.remove(str(final_output_csv))
    finally:
        if temp_root.exists():
            shutil.rmtree(str(temp_root))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 run_MIND.py <subject_id>")
        sys.exit(1)
    
    subj_arg = sys.argv[1]
    run_single_subject(subj_arg)