"""Read only saved spike fields; leave voltage records and simulation timing untouched."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from scipy.io import loadmat

p = argparse.ArgumentParser()
p.add_argument('--trial', required=True)
p.add_argument('--duration', type=int, required=True)
p.add_argument('--start', type=int, default=0)
p.add_argument('--out', required=True)
p.add_argument('--reference')
a = p.parse_args()
trial = Path(a.trial)
out = Path(a.out)
out.mkdir(parents=True, exist_ok=True)
inputs = '/scratch/xh2906/librarySCI_runs/figure1g2_dense_cold_equilibrium_20260829/snn_paper2_original_exact_20260901_160003/inputs/AllMFPixPara_Paper2TuneFig1V4D2_torch_staged.mat'
meta = loadmat(inputs, variable_names=['NnS','NnC','NnI','EcplxInd','N_E','N_I','N_HC','NPixX','NPixY'], simplify_cells=True)
assert meta['N_HC']==3 and meta['NPixX']==10 and meta['NPixY']==10
ecplx = meta['EcplxInd'].astype(bool)
ids = [np.flatnonzero(~ecplx),np.flatnonzero(ecplx),np.arange(meta['N_I'])]
pixel = []
for name in ['NnS','NnC','NnI']:
    n=meta[name]
    x=np.ceil(n['X']/(n['X'].max()/30)).astype(int)
    y=np.ceil(n['Y']/(n['Y'].max()/30)).astype(int)
    pixel.append(y-1+(x-1)*30)
neurons=np.stack([np.bincount(v,minlength=900) for v in pixel])
assert np.all(neurons>0)
counts=np.zeros((3,900,a.duration-a.start),dtype=np.int64)
last9=np.zeros((3,900),dtype=np.int64)
data=trial/'Data/Paper2_NetworkTuning/Fig1V4/Paper2NWSimulationData'
for second in range(a.start+1,a.duration+1):
    path=data/f'DriveWkSp_SCSepa_Cconst_4Hz_Deg0.0_{second}s_NewSmear.mat'
    totals=[]; final=[]
    with h5py.File(path,'r') as f:
        for field,n in [('SpEs',meta['N_E']),('SpIs',meta['N_I'])]:
            refs=f['NWTrace/'+field][()].ravel()
            assert len(refs)==a.duration*20
            # Original SampleInd increases across seconds; only these 20 cells are populated.
            selected=refs[(second-1)*20:second*20]
            totals.append(np.zeros(n,dtype=np.int64)); final.append(np.zeros(n,dtype=np.int64))
            for ref in selected:
                d=f[ref]
                if d.attrs.get('MATLAB_empty',0): continue
                spikes=d[()]
                assert spikes.shape[0]==2
                ix=spikes[0].astype(np.int64)-1
                assert np.all((ix>=0)&(ix<n)) and np.all((spikes[1]>=0)&(spikes[1]<=1000))
                totals[-1]+=np.bincount(ix,minlength=n)
                if second>=a.duration-9:
                    times=spikes[1]+np.float32((second-1)*1000)
                    use=(times>=(a.duration*1000-9000))&(times<=a.duration*1000)
                    final[-1]+=np.bincount(ix[use],minlength=n)
    for k,j in enumerate([0,0,1]):
        counts[k,:,second-a.start-1]=np.bincount(pixel[k],weights=totals[j][ids[k]],minlength=900).astype(np.int64)
        last9[k]+=np.bincount(pixel[k],weights=final[j][ids[k]],minlength=900).astype(np.int64)
    if second%25==0: print(f'Read {second}/{a.duration}',flush=True)
if a.reference:
    previous=loadmat(a.reference,simplify_cells=True)
    old=np.stack([np.asarray(v) for v in previous['allPixelCounts']])
    assert old.shape==counts.shape
    assert np.array_equal(old,counts), 'Spike-field reader differs from verified MATLAB extraction'
    print('PASS exact equality to all MATLAB pixel counts',flush=True)
saved=loadmat(data/'NWSimulationPix_0.0deg_NewSmear.mat',simplify_cells=True)['NWSmlt']
errors={name:float(np.max(np.abs(last9[k]/neurons[k]/9-saved[name]))) for k,name in enumerate(['FS','FC','FI'])}
assert max(errors.values())<1e-10,errors
np.savez_compressed(out/'counts.npz',counts=counts,neurons=neurons,duration=a.duration,start=a.start)
(out/'summary.json').write_text(json.dumps(dict(duration=a.duration,start=a.start,last9_max_abs=errors,shape=list(counts.shape),reference_verified=bool(a.reference)),indent=2)+'\n')
print(json.dumps(errors),flush=True)
