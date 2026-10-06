"""Joint interference-order S,D fit, SM Asimov data; no detector likelihood is inferred.
Requires NumPy. Uses the included order-96 and order-128 response integrals.
Parameters x=n7/1e-4, y=C_Hq3/(0.01 TeV^-2). Channels ordered Z1,h1,Z2,h2.
Covariance models are illustrative fixed-error assumptions, not measured covariances.
"""
from pathlib import Path
import json,csv
import numpy as np
P=Path(__file__).resolve().parent
r=json.loads((P/'response_integrals_order128.json').read_text())
r96=json.loads((P/'response_integrals_order96.json').read_text())
convergence=max(abs(a/b-1) for k in r for aa,bb in zip(r96[k],r[k]) for a,b in zip(aa,bb) if abs(b)>1e-28)
assert convergence<2e-5
T=np.kron(np.eye(2),np.array([[.5,.5],[.5,-.5]]))
e=np.array([.088,.15,.091,.15]);scales=np.array([1e-4,.01])
def jacobian(p,cut):
 J=[]
 for i in range(2):
  for ch in ['WZ','Wh']:
   a=np.array(r[f'{ch}_long_{p}_{cut}'][i]);J.append(a[1:]/a[0])
 return np.array(J)
def covariance(rho=0,common=0):
 V=np.diag(e**2)
 if common:
  assert common<min(e)
  # Preserve total diagonal errors: allocate a common component, not add it twice.
  V=np.diag(e**2-common**2)+common**2*np.ones((4,4))
 else:
  for i in range(2):V[2*i,2*i+1]=V[2*i+1,2*i]=rho*e[2*i]*e[2*i+1]
 assert np.linalg.eigvalsh(V).min()>0
 return V
def fit(J,V):
 A=T@J;W=T@V@T.T
 F=J.T@np.linalg.solve(V,J)
 Fsd=A.T@np.linalg.solve(W,A)
 assert np.allclose(F,Fsd,rtol=1e-11,atol=1e-12)
 C=np.linalg.inv(F)
 sd=np.sqrt(np.diag(C));fixed=1/np.sqrt(np.diag(F))
 for k in range(2):
  theta=np.zeros(2);theta[k]=sd[k];j=1-k;theta[j]=-F[j,k]/F[j,j]*theta[k]
  assert abs(theta@F@theta-1)<1e-12
 ang=np.linspace(0,2*np.pi,1000)
 ellipse=np.sqrt(5.991)*np.linalg.cholesky(C)@np.vstack([np.cos(ang),np.sin(ang)])
 assert np.min(1+J@ellipse)>0
 return dict(F=F.tolist(),C=C.tolist(),A=A.tolist(),Vsd=W.tolist(),
  sigma_profile=(sd*scales).tolist(),sigma_fixed=(fixed*scales).tolist(),
  parameter_correlation=float(C[0,1]/np.prod(sd)))
models=[('independent',0,0),('rho=-0.5',-.5,0),('rho=+0.5',.5,0),('common=3%',0,.03),('common=5%',0,.05)]
results={};rows=[];lines=['Joint n7 / C_Hq3 fit through S_i and D_i','',
 'SM Asimov, interference-only rates; kappa7=1; all other deformations fixed.',
 'pT bins 100-200 and 200-1500 GeV; same upper mass cutoff in both channels.',
 'C_Hq3 in TeV^-2. Normalized response coordinates x=n7/1e-4, y=C_Hq3/0.01.',
 'Fixed fractional channel errors: WZ=(0.088,0.091), Wh=(0.15,0.15).',
 'Baseline has no cross-channel or cross-bin experimental covariance.',
 'Reported 1D intervals profile over the other coefficient (Delta chi2=1).',
 '2D contours use Delta chi2=2.30 and 5.991 (nominal 68%,95%).',
 'This is a two-deformation benchmark, not a general SMEFT-versus-HEFT likelihood test.',
 'Source flavour/charge content is preserved (source W+ luminosity convention); no detector or new charge-summed model is generated.',
 'At linear order the two derivative templates may be added without double-counting SM.',
 'No joint quadratic/full-amplitude model is claimed: it requires consistent mixed coefficient terms.',
 'Total errors remain fixed after changing cuts and extending pT acceptance.',
 'Correlation sweeps below are assumptions, not extracted nuisance estimates.','']
for p in [1,2]:
 for cut in [500,1000,2000,3000]:
  for name,rho,common in models:
   J=jacobian(p,cut);V=covariance(rho,common);f=fit(J,V)
   if common:
    # Profile one N(0,1) common nuisance and verify equivalence to marginal covariance.
    D=np.diag(e**2-common**2);Di=np.linalg.inv(D);b=np.full(4,common)
    prof=Di-np.outer(Di@b,Di@b)/(1+b@Di@b)
    assert np.allclose(J.T@prof@J,np.array(f['F']),rtol=1e-11,atol=1e-12)
   key=f'{p}_{cut}_{name}';results[key]=f
   rows.append([p,cut,name,*f['sigma_profile'],*f['sigma_fixed'],f['parameter_correlation']])
  f=results[f'{p}_{cut}_independent']
  lines.append(f'PDF {p}, Mmax={cut}: profiled +/- {f["sigma_profile"]}; other coefficient fixed +/- {f["sigma_fixed"]}; parameter correlation={f["parameter_correlation"]:.8f}')
lines+=['','3 TeV nCTEQ15 response matrix (rows S1,D1,S2,D2; columns x,y):',str(np.array(results['2_3000_independent']['A'])),
 'S,D covariance (rows S1,D1,S2,D2):',str(np.array(results['2_3000_independent']['Vsd'])), '',
 '3 TeV nCTEQ15 correlation sensitivity (profiled sigmas n7,C_Hq3):']
for name,_,_ in models:
 f=results[f'2_3000_{name}'];lines.append(f'{name}: {f["sigma_profile"]}, coefficient correlation {f["parameter_correlation"]:.6f}')
# S-only and D-only Fisher matrices: informative diagnostics, not extra independent data.
J=jacobian(2,3000);V=covariance();A=T@J;W=T@V@T.T
for name,inds in [('S',[0,2]),('D',[1,3])]:
 a=A[inds,:];v=W[np.ix_(inds,inds)];F=a.T@np.linalg.solve(v,a)
 results[f'{name}_only_2_3000']={'F':F.tolist(),'condition_number':float(np.linalg.cond(F))}
 # Their very weak profiled directions need not remain physically meaningful far from SM.
 lines.append(f'{name}-only Fisher condition number: {np.linalg.cond(F):.6g}; not used for reported joint limits.')
lines+=['',f'Quadrature convergence max relative coefficient change: {convergence:.6g}',
 'Checks passed: channel-space and S,D Fisher matrices identical; profiled endpoints chi2=1;',
 '95% joint ellipses predict positive integrated rates; illustrative common-nuisance profiling equals covariance treatment.',
 'S_i and D_i must be fitted with their covariance; adding a separate ratio or WZ-fraction likelihood would require their shared-data covariance.']
(P/'results.txt').write_text('\n'.join(lines)+'\n')
(P/'fit_results.json').write_text(json.dumps(results,indent=2)+'\n')
with (P/'sensitivity.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['PDF set','Mmax GeV','covariance model','profile sigma n7','profile sigma C TeV^-2','fixed sigma n7','fixed sigma C TeV^-2','fit correlation']);w.writerows(rows)
# Export exact numeric matrices for a native Mathematica scientific figure.
fmt=lambda a:'{'+','.join(fmt(x) if isinstance(x,list) else format(x,'.17g') for x in a)+'}'
(P/'plot_matrices.wl').write_text('\n'.join(f'{n}={fmt(results[k]["F"])};' for n,k in [('fJoint','2_3000_independent'),('fS','S_only_2_3000'),('fD','D_only_2_3000')])+'\n')
print('\n'.join(lines))
