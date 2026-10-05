"""Collect lightweight Torch tables without polling Slurm or transferring state MATs."""
import io
import subprocess
import tarfile
from pathlib import Path
from analyze_impl import analyze

ROOT = Path(__file__).resolve().parents[1]
REMOTE = '/scratch/xh2906/librarySCI_runs/figure1b_dnn3_impl_optimization_20260921'
SCRIPT = r'''
import io, pathlib, sys, tarfile
root = pathlib.Path(REMOTE_ROOT)
with tarfile.open(fileobj=sys.stdout.buffer, mode='w|') as archive:
    for phase in ['probe', 'validation', 'control']:
        for path in sorted((root / phase).glob('*.tsv')):
            archive.add(path, arcname=str(path.relative_to(root)))
    for name in ['selected_variant.txt', 'submission.txt']:
        path = root / name
        if path.exists():
            archive.add(path, arcname=name)
'''.replace('REMOTE_ROOT', repr(REMOTE))

result = subprocess.run(['ssh', '-S', '/Users/xialeihuang/.ssh/torch-codex-cm',
                         '-o', 'BatchMode=yes', 'xh2906@login.torch.hpc.nyu.edu', 'python3 -'],
                        input=SCRIPT.encode(), capture_output=True, check=True)
with tarfile.open(fileobj=io.BytesIO(result.stdout)) as archive:
    for item in archive:
        relative = Path(item.name)
        assert item.isfile() and not relative.is_absolute() and '..' not in relative.parts
        target = ROOT / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(archive.extractfile(item).read())
analyze(ROOT)
