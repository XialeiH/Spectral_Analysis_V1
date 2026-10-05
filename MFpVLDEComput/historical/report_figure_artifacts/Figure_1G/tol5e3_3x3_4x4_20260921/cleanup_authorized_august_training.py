import datetime, json, os, pathlib, shutil
root = pathlib.Path("/scratch/xh2906/librarySCI_runs")
parents = ["grating50_angle0_global_residual_20260817_140505","grating50_wide_k8_angle0_pcgrad_20260817_120942","grating_ann_architecture_benchmark_20260813_205219","grating_capacity_symmetry_study_20260813_204040","grating_image_cnn_convrnn_SCI_20260812_233101","grating_image_cnn_convrnn_SCI_clean_20260812_233702","grating_image_cnn_convrnn_SCI_realcontrast_20260812_234033","grating_resolution50_coordinatefree_cnn_pinn_20260817_094000","grating_resolution50_mlp_encoder_convrnn_20260817_160615","grating_resolution50_positionfree_architecture_benchmark_20260816_222900"]
cutoff = datetime.datetime(2026,9,1,tzinfo=datetime.timezone.utc).timestamp()
plan = []
for name in parents:
    target = root / name / "training"
    assert not target.is_symlink() and target.is_dir(), str(target)
    assert target.resolve() == target
    latest = target.stat().st_mtime
    allocated = 0
    entries = 0
    for base, dirs, files in os.walk(target, followlinks=False):
        for child in dirs + files:
            st = os.lstat(os.path.join(base,child))
            latest = max(latest,st.st_mtime)
            allocated += st.st_blocks*512
            entries += 1
    assert latest < cutoff, "Recent changes: "+str(target)
    plan.append({"path":str(target),"bytes_before":allocated,"entries":entries,"latest_mtime":latest})
for item in plan:
    shutil.rmtree(item["path"])
    item["deleted"] = not pathlib.Path(item["path"]).exists()
    print(json.dumps(item), flush=True)
print(json.dumps({"total_bytes_before":sum(p["bytes_before"] for p in plan),"deleted_directories":len(plan)}))

